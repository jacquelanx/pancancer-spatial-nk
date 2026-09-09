"""Anthropic client wrapper: cached system anchor, streaming, and a provenance ledger.

Design points that matter for a study that has to be defensible:

* The system anchor (``prompts/system/00_meta.md``) is byte-identical on every call and is
  the cached prefix, so the statistical constraints are always in force and cheap.
* Every call writes a run record under ``runs/`` containing the prompt fingerprint, the
  bound variable values, token usage, and the response. Any generated sentence can be
  traced to the inputs that produced it.
* Responses are never silently truncated: ``stop_reason == "max_tokens"`` raises.
"""

from __future__ import annotations

import json
import logging
import os
import time
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from panspatial.orchestrator.registry import RenderedPrompt

log = logging.getLogger("panspatial.orchestrator")

# USD per million tokens, Claude API list price. Used only for the cost line in run records.
_PRICING = {
    "claude-opus-5": (5.0, 25.0),
    "claude-sonnet-5": (2.0, 10.0),
    "claude-haiku-4-5": (1.0, 5.0),
}
# Above this, the SDK wants streaming or the request risks an HTTP timeout.
_STREAM_THRESHOLD = 16000


class LLMError(RuntimeError):
    """Raised when a call fails or returns output that must not be used downstream."""


@dataclass
class LLMResult:
    fingerprint: str
    template_id: str
    template_version: str
    model: str
    effort: str
    text: str
    parsed: Any | None
    stop_reason: str | None
    input_tokens: int
    output_tokens: int
    cache_read_tokens: int
    cache_creation_tokens: int
    cost_usd: float
    duration_s: float
    request_id: str | None
    created_at: str

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _estimate_cost(model: str, usage: Any) -> float:
    price_in, price_out = _PRICING.get(model, (0.0, 0.0))
    fresh = getattr(usage, "input_tokens", 0) or 0
    cached_read = getattr(usage, "cache_read_input_tokens", 0) or 0
    cached_write = getattr(usage, "cache_creation_input_tokens", 0) or 0
    out = getattr(usage, "output_tokens", 0) or 0
    # Cache reads bill at ~0.1x, cache writes at ~1.25x of the input rate.
    return (
        fresh * price_in
        + cached_read * price_in * 0.1
        + cached_write * price_in * 1.25
        + out * price_out
    ) / 1_000_000


class OrchestratorClient:
    """Executes a :class:`RenderedPrompt` against the Claude API and records the run."""

    def __init__(
        self,
        run_dir: str | Path = "runs",
        *,
        client: Any | None = None,
        dry_run: bool = False,
    ) -> None:
        self.run_dir = Path(run_dir)
        self.run_dir.mkdir(parents=True, exist_ok=True)
        self.dry_run = dry_run
        self._client = client
        if client is None and not dry_run:
            try:
                import anthropic
            except ImportError as exc:  # pragma: no cover - environment dependent
                raise LLMError(
                    "the `anthropic` package is required for live calls; "
                    "`pip install anthropic`, or use dry_run=True to render prompts only"
                ) from exc
            # Zero-arg construction resolves ANTHROPIC_API_KEY, ANTHROPIC_AUTH_TOKEN,
            # or an `ant auth login` profile, in that order.
            self._client = anthropic.Anthropic(timeout=1800.0, max_retries=4)

    # ------------------------------------------------------------------ requests

    def _request_kwargs(self, prompt: RenderedPrompt, schema: dict | None) -> dict[str, Any]:
        output_config: dict[str, Any] = {"effort": prompt.effort}
        if schema is not None:
            output_config["format"] = {"type": "json_schema", "schema": schema}
        return {
            "model": prompt.model,
            "max_tokens": prompt.max_tokens,
            # Frozen prefix -> cached. Never interpolate anything per-run into `system`.
            "system": [
                {
                    "type": "text",
                    "text": prompt.system,
                    "cache_control": {"type": "ephemeral", "ttl": "1h"},
                }
            ],
            "messages": [{"role": "user", "content": prompt.user}],
            "thinking": {"type": "adaptive"},
            "output_config": output_config,
        }

    def run(
        self,
        prompt: RenderedPrompt,
        *,
        schema: dict | None = None,
        write_record: bool = True,
    ) -> LLMResult:
        if self.dry_run:
            return self._dry_run_result(prompt)

        kwargs = self._request_kwargs(prompt, schema)
        started = time.monotonic()
        log.info(
            "calling %s [%s v%s, effort=%s, fp=%s]",
            prompt.model,
            prompt.template_id,
            prompt.template_version,
            prompt.effort,
            prompt.fingerprint,
        )
        response = self._call(kwargs)
        duration = time.monotonic() - started

        if response.stop_reason == "refusal":
            details = getattr(response, "stop_details", None)
            raise LLMError(
                f"{prompt.template_id}: request declined "
                f"({getattr(details, 'category', 'unknown')}): "
                f"{getattr(details, 'explanation', '')}"
            )
        if response.stop_reason == "max_tokens":
            raise LLMError(
                f"{prompt.template_id}: response hit max_tokens ({prompt.max_tokens}) and is "
                "truncated. Raise max_tokens in the template rather than using partial output."
            )

        text = "".join(b.text for b in response.content if b.type == "text")
        parsed = None
        if schema is not None:
            try:
                parsed = json.loads(text)
            except json.JSONDecodeError as exc:
                raise LLMError(f"{prompt.template_id}: structured output is not valid JSON") from exc

        usage = response.usage
        result = LLMResult(
            fingerprint=prompt.fingerprint,
            template_id=prompt.template_id,
            template_version=prompt.template_version,
            model=prompt.model,
            effort=prompt.effort,
            text=text,
            parsed=parsed,
            stop_reason=response.stop_reason,
            input_tokens=getattr(usage, "input_tokens", 0) or 0,
            output_tokens=getattr(usage, "output_tokens", 0) or 0,
            cache_read_tokens=getattr(usage, "cache_read_input_tokens", 0) or 0,
            cache_creation_tokens=getattr(usage, "cache_creation_input_tokens", 0) or 0,
            cost_usd=_estimate_cost(prompt.model, usage),
            duration_s=round(duration, 2),
            request_id=getattr(response, "_request_id", None),
            created_at=datetime.now(timezone.utc).isoformat(),
        )
        log.info(
            "done fp=%s in %.1fs, %d out tok, cache_read=%d, $%.4f",
            result.fingerprint,
            result.duration_s,
            result.output_tokens,
            result.cache_read_tokens,
            result.cost_usd,
        )
        if result.cache_read_tokens == 0 and result.cache_creation_tokens == 0:
            log.warning(
                "no cache activity on fp=%s -- the system anchor should be a cached prefix; "
                "check that nothing per-run leaked into `system`",
                result.fingerprint,
            )
        if write_record:
            self._write_record(prompt, result)
        return result

    def _call(self, kwargs: dict[str, Any]) -> Any:
        import anthropic

        try:
            if kwargs["max_tokens"] > _STREAM_THRESHOLD:
                with self._client.messages.stream(**kwargs) as stream:
                    return stream.get_final_message()
            return self._client.messages.create(**kwargs)
        except anthropic.NotFoundError as exc:
            raise LLMError(f"unknown model {kwargs['model']!r}: {exc}") from exc
        except anthropic.AuthenticationError as exc:
            raise LLMError(
                "authentication failed: set ANTHROPIC_API_KEY or run `ant auth login`"
            ) from exc
        except anthropic.RateLimitError as exc:
            raise LLMError(f"rate limited after SDK retries: {exc}") from exc
        except anthropic.BadRequestError as exc:
            raise LLMError(f"request rejected: {exc}") from exc
        except anthropic.APIConnectionError as exc:
            raise LLMError(f"network failure reaching the API: {exc}") from exc

    def _dry_run_result(self, prompt: RenderedPrompt) -> LLMResult:
        log.info("dry run: %s fp=%s (%d chars)", prompt.template_id, prompt.fingerprint, len(prompt.user))
        return LLMResult(
            fingerprint=prompt.fingerprint,
            template_id=prompt.template_id,
            template_version=prompt.template_version,
            model=prompt.model,
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
                "model": prompt.model,
                "effort": prompt.effort,
                "bound_values": prompt.bound_values,
                "user_prompt": prompt.user,
                "system_sha256_prefix": prompt.fingerprint,
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
