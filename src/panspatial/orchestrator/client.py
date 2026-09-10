"""Provider-independent execution layer: runs a rendered prompt and records the run.

Everything that makes this pipeline defensible lives here or above, not in the provider:

* The system anchor is byte-identical on every call and sits first in the request, so the
  statistical constraints are always in force -- and so both providers' prompt caching
  (explicit on Anthropic, automatic prefix-based on Azure) actually hits.
* Every call writes a run record under ``runs/`` containing the prompt fingerprint, the
  bound variable values, token usage, and the response, so any generated sentence can be
  traced to the inputs that produced it.
* Responses are never silently truncated or silently filtered; the backend raises.

Provider selection is in :mod:`panspatial.orchestrator.backends`.
"""

from __future__ import annotations

import json
import logging
import os
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from panspatial.orchestrator.backends import (
    Backend,
    BackendResponse,
    LLMError,
    resolve_backend,
)
from panspatial.orchestrator.registry import RenderedPrompt

log = logging.getLogger("panspatial.orchestrator")

__all__ = ["LLMError", "LLMResult", "OrchestratorClient"]


@dataclass
class LLMResult:
    fingerprint: str
    template_id: str
    template_version: str
    provider: str
    model: str
    effort: str
    text: str
    parsed: Any | None
    stop_reason: str | None
    input_tokens: int
    output_tokens: int
    cache_read_tokens: int
    cache_creation_tokens: int
    cost_usd: float | None
    duration_s: float
    request_id: str | None
    created_at: str
    extra: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @property
    def cost_display(self) -> str:
        """Cost is only shown where list pricing is actually known.

        Azure pricing depends on region, tier, and commercial agreement, so the run record
        reports tokens and says so rather than printing a number that is probably wrong.
        """
        return f"${self.cost_usd:.4f}" if self.cost_usd is not None else "unpriced"


class OrchestratorClient:
    """Executes a :class:`RenderedPrompt` against a backend and records the run."""

    def __init__(
        self,
        run_dir: str | Path = "runs",
        *,
        backend: Backend | str | None = None,
        dry_run: bool = False,
    ) -> None:
        self.run_dir = Path(run_dir)
        self.run_dir.mkdir(parents=True, exist_ok=True)
        self.dry_run = dry_run
        if dry_run:
            self.backend: Backend | None = backend if isinstance(backend, Backend) else None
        elif isinstance(backend, Backend):
            self.backend = backend
        else:
            self.backend = resolve_backend(backend)
        if self.backend is not None and not dry_run:
            log.info("backend: %s", self.backend.describe())

    # ------------------------------------------------------------------ requests

    def run(
        self,
        prompt: RenderedPrompt,
        *,
        schema: dict | None = None,
        write_record: bool = True,
    ) -> LLMResult:
        if self.dry_run:
            return self._dry_run_result(prompt)
        assert self.backend is not None

        target = self.backend.resolve_model(prompt.model)
        log.info(
            "calling %s [%s v%s, model=%s, effort=%s, fp=%s]",
            self.backend.name,
            prompt.template_id,
            prompt.template_version,
            target,
            prompt.effort,
            prompt.fingerprint,
        )
        started = time.monotonic()
        response = self.backend.complete(prompt, schema=schema)
        duration = time.monotonic() - started

        result = self._to_result(prompt, response, duration)
        log.info(
            "done fp=%s in %.1fs, %d out tok, cached_in=%d, %s",
            result.fingerprint,
            result.duration_s,
            result.output_tokens,
            result.cache_read_tokens,
            result.cost_display,
        )
        if result.cache_read_tokens == 0 and result.cache_creation_tokens == 0:
            log.warning(
                "no cache activity on fp=%s -- the system anchor should be a cached prefix; "
                "check that nothing per-run leaked into the system message",
                result.fingerprint,
            )
        if write_record:
            self._write_record(prompt, result)
        return result

    def _to_result(
        self, prompt: RenderedPrompt, response: BackendResponse, duration: float
    ) -> LLMResult:
        return LLMResult(
            fingerprint=prompt.fingerprint,
            template_id=prompt.template_id,
            template_version=prompt.template_version,
            provider=response.provider,
            model=response.model,
            effort=prompt.effort,
            text=response.text,
            parsed=response.parsed,
            stop_reason=response.stop_reason,
            input_tokens=response.input_tokens,
            output_tokens=response.output_tokens,
            cache_read_tokens=response.cached_input_tokens,
            cache_creation_tokens=response.cache_write_tokens,
            cost_usd=response.cost_usd,
            duration_s=round(duration, 2),
            request_id=response.request_id,
            created_at=datetime.now(timezone.utc).isoformat(),
            extra=dict(response.extra),
        )

    def _dry_run_result(self, prompt: RenderedPrompt) -> LLMResult:
        target = self.backend.resolve_model(prompt.model) if self.backend else prompt.model
        log.info(
            "dry run: %s fp=%s -> %s (%d chars)",
            prompt.template_id,
            prompt.fingerprint,
            target,
            len(prompt.user),
        )
        return LLMResult(
            fingerprint=prompt.fingerprint,
            template_id=prompt.template_id,
            template_version=prompt.template_version,
            provider=self.backend.name if self.backend else "dry-run",
            model=target,
            effort=prompt.effort,
            text="",
            parsed=None,
            stop_reason="dry_run",
            input_tokens=0,
            output_tokens=0,
            cache_read_tokens=0,
            cache_creation_tokens=0,
            cost_usd=0.0,
            duration_s=0.0,
            request_id=None,
            created_at=datetime.now(timezone.utc).isoformat(),
        )

    # ------------------------------------------------------------------ ledger

    def _write_record(self, prompt: RenderedPrompt, result: LLMResult) -> Path:
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
        slug = prompt.template_id.replace(".", "_")
        path = self.run_dir / f"{stamp}_{slug}_{result.fingerprint}.json"
        record = {
            "prompt": {
                "id": prompt.template_id,
                "version": prompt.template_version,
                "fingerprint": prompt.fingerprint,
                "requested_model": prompt.model,
                "resolved_model": result.model,
                "provider": result.provider,
                "effort": prompt.effort,
                "bound_values": prompt.bound_values,
                "user_prompt": prompt.user,
            },
            "result": result.to_dict(),
            "env": {
                "git_commit": os.environ.get("PANSPATIAL_GIT_COMMIT", "unrecorded"),
                "host": os.uname().nodename,
            },
        }
        path.write_text(json.dumps(record, indent=2, sort_keys=True), encoding="utf-8")
        log.info("run record -> %s", path)
        return path
