"""Provider backends for the orchestrator: Anthropic and Azure OpenAI.

The orchestration guarantees -- strict prompt binding, the frozen system anchor, the
verified-results gate, run records -- live above this layer and are provider-independent.
A backend's only job is to send one rendered prompt and normalize what comes back into a
:class:`BackendResponse`, including the two failure modes that must never pass silently:

* **Truncation.** ``finish_reason == "length"`` (Azure) or ``stop_reason == "max_tokens"``
  (Anthropic) produces a response that reads complete and is not. Both raise.
* **Policy refusal.** Azure's content filter and Anthropic's refusal stop reason are
  normalized to the same error, so a blocked call cannot be mistaken for an empty result.

Azure OpenAI has three quirks this module absorbs so callers never see them:

1. ``model`` is a *deployment name*, not a model name, and deployment names are arbitrary.
   ``AZURE_OPENAI_DEPLOYMENTS`` maps the logical model in a prompt template to whatever the
   deployment is actually called.
2. Reasoning deployments reject ``max_tokens`` (they want ``max_completion_tokens``) and
   accept ``reasoning_effort``; chat deployments are the reverse. Which is which cannot be
   inferred from a deployment name, so capabilities are probed once per deployment from the
   API's own 400 message and cached for the process.
3. Prompt caching is automatic and prefix-based above ~1024 tokens -- there is no
   ``cache_control``. The frozen system anchor is therefore still the right design; it is
   just cached implicitly rather than explicitly.
"""

from __future__ import annotations

import json
import logging
import os
import re
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from panspatial.orchestrator.registry import RenderedPrompt

log = logging.getLogger("panspatial.orchestrator.backend")

# USD per million tokens for Anthropic list pricing. Azure pricing varies by region and
# agreement, so Azure calls report tokens and leave cost unpriced rather than guessing.
ANTHROPIC_PRICING = {
    "claude-opus-5": (5.0, 25.0),
    "claude-sonnet-5": (2.0, 10.0),
    "claude-haiku-4-5": (1.0, 5.0),
}

# Above this, the Anthropic SDK wants streaming or the request risks an HTTP timeout.
_STREAM_THRESHOLD = 16000

# Our templates use Anthropic's five levels; Azure reasoning models take three.
_EFFORT_TO_AZURE = {
    "low": "low",
    "medium": "medium",
    "high": "high",
    "xhigh": "high",
    "max": "high",
}


class LLMError(RuntimeError):
    """A call failed, or returned output that must not be used downstream."""


@dataclass
class BackendResponse:
    """Normalized result of one completion, whatever provider produced it."""

    text: str
    parsed: Any | None
    model: str
    stop_reason: str | None
    input_tokens: int
    output_tokens: int
    cached_input_tokens: int
    cache_write_tokens: int
    cost_usd: float | None
    request_id: str | None
    provider: str
    extra: dict[str, Any] = field(default_factory=dict)


# --------------------------------------------------------------------------- env


def load_dotenv(path: str | Path = ".env", *, override: bool = False) -> dict[str, str]:
    """Load ``KEY=value`` pairs from a dotenv file into ``os.environ``.

    Deliberately dependency-free and deliberately conservative: existing environment
    variables win unless ``override`` is set, so a shell export always beats a stale file.
    Values may be quoted; ``export`` prefixes and ``#`` comments are tolerated.
    """
    file = Path(path)
    if not file.is_file():
        return {}
    loaded: dict[str, str] = {}
    for raw in file.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("export "):
            line = line[len("export ") :].lstrip()
        key, sep, value = line.partition("=")
        if not sep:
            continue
        key = key.strip()
        value = value.strip().strip('"').strip("'")
        if not key or (key in os.environ and not override):
            continue
        os.environ[key] = value
        loaded[key] = value
    if loaded:
        log.info("loaded %d variable(s) from %s: %s", len(loaded), file, ", ".join(sorted(loaded)))
    return loaded


def parse_deployment_map(spec: str | None) -> dict[str, str]:
    """Parse ``AZURE_OPENAI_DEPLOYMENTS`` into a logical-model -> deployment mapping.

    Two accepted forms, mixable::

        AZURE_OPENAI_DEPLOYMENTS=gpt-5-prod
        AZURE_OPENAI_DEPLOYMENTS=claude-opus-5=gpt-5-prod,fast=gpt-4o-mini-dev

    A bare entry maps to itself. The first entry is the default deployment, used for any
    model a prompt template requests that has no explicit mapping -- which is what lets the
    same templates run on either provider without edits.
    """
    mapping: dict[str, str] = {}
    if not spec:
        return mapping
    for chunk in spec.split(","):
        item = chunk.strip()
        if not item:
            continue
        logical, sep, deployment = item.partition("=")
        logical, deployment = logical.strip(), deployment.strip()
        if sep and deployment:
            mapping[logical] = deployment
        else:
            mapping[logical] = logical
    return mapping


# ----------------------------------------------------------------------- interface


class Backend(ABC):
    """One provider. Implementations normalize errors into :class:`LLMError`."""

    name: str

    @abstractmethod
    def resolve_model(self, requested: str) -> str:
        """Map the model named in a prompt template to what this provider expects."""

    @abstractmethod
    def complete(self, prompt: RenderedPrompt, schema: dict | None = None) -> BackendResponse:
        """Send one rendered prompt and return the normalized response."""

    def describe(self) -> str:
        return self.name


# ----------------------------------------------------------------------- anthropic


class AnthropicBackend(Backend):
    """Claude via the Anthropic API, with the system anchor as an explicit cached prefix."""

    name = "anthropic"

    def __init__(self, client: Any | None = None, *, timeout: float = 1800.0) -> None:
        if client is not None:
            self._client = client
            return
        try:
            import anthropic
        except ImportError as exc:  # pragma: no cover - environment dependent
            raise LLMError(
                "the `anthropic` package is required for this backend: pip install anthropic"
            ) from exc
        # Zero-arg construction resolves ANTHROPIC_API_KEY, ANTHROPIC_AUTH_TOKEN, or an
        # `ant auth login` profile, in that order.
        self._client = anthropic.Anthropic(timeout=timeout, max_retries=4)

    def resolve_model(self, requested: str) -> str:
        return requested

    def complete(self, prompt: RenderedPrompt, schema: dict | None = None) -> BackendResponse:
        import anthropic

        output_config: dict[str, Any] = {"effort": prompt.effort}
        if schema is not None:
            output_config["format"] = {"type": "json_schema", "schema": schema}
        kwargs: dict[str, Any] = {
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

        try:
            if prompt.max_tokens > _STREAM_THRESHOLD:
                with self._client.messages.stream(**kwargs) as stream:
                    response = stream.get_final_message()
            else:
                response = self._client.messages.create(**kwargs)
        except anthropic.NotFoundError as exc:
            raise LLMError(f"unknown model {prompt.model!r}: {exc}") from exc
        except anthropic.AuthenticationError as exc:
            raise LLMError(
                "authentication failed: set ANTHROPIC_API_KEY or run `ant auth login`"
            ) from exc
        except anthropic.RateLimitError as exc:
            raise LLMError(f"rate limited after SDK retries: {exc}") from exc
        except anthropic.BadRequestError as exc:
            raise LLMError(f"request rejected: {exc}") from exc
        except anthropic.APIConnectionError as exc:
            raise LLMError(f"network failure reaching the Anthropic API: {exc}") from exc

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
        usage = response.usage
        return BackendResponse(
            text=text,
            parsed=_parse_json(text, prompt.template_id) if schema is not None else None,
            model=prompt.model,
            stop_reason=response.stop_reason,
            input_tokens=getattr(usage, "input_tokens", 0) or 0,
            output_tokens=getattr(usage, "output_tokens", 0) or 0,
            cached_input_tokens=getattr(usage, "cache_read_input_tokens", 0) or 0,
            cache_write_tokens=getattr(usage, "cache_creation_input_tokens", 0) or 0,
            cost_usd=_anthropic_cost(prompt.model, usage),
            request_id=getattr(response, "_request_id", None),
            provider=self.name,
        )


def _anthropic_cost(model: str, usage: Any) -> float:
    price_in, price_out = ANTHROPIC_PRICING.get(model, (0.0, 0.0))
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


# -------------------------------------------------------------------- azure openai


@dataclass
class DeploymentCapabilities:
    """What one Azure deployment actually accepts, learned from the API's own errors."""

    token_param: str = "max_completion_tokens"
    send_reasoning_effort: bool = True
    structured_output: str = "json_schema"  # -> "json_object" -> "prompt_only"

    def describe(self) -> str:
        return (
            f"{self.token_param}"
            f"{', reasoning_effort' if self.send_reasoning_effort else ''}"
            f", {self.structured_output}"
        )


_UNSUPPORTED = re.compile(
    r"unsupported|unrecognized|not supported|invalid[_ ]?(?:request[_ ]?)?(?:argument|parameter)"
    r"|unknown parameter|extra inputs are not permitted",
    re.IGNORECASE,
)


class AzureOpenAIBackend(Backend):
    """Azure OpenAI via the ``openai`` SDK's ``AzureOpenAI`` client.

    Requires ``AZURE_OPENAI_API_KEY``, ``AZURE_OPENAI_ENDPOINT``,
    ``AZURE_OPENAI_API_VERSION``, and ``AZURE_OPENAI_DEPLOYMENTS``.
    """

    name = "azure-openai"

    def __init__(
        self,
        client: Any | None = None,
        *,
        endpoint: str | None = None,
        api_key: str | None = None,
        api_version: str | None = None,
        deployments: str | dict[str, str] | None = None,
        timeout: float = 1800.0,
    ) -> None:
        if isinstance(deployments, dict):
            self.deployments = dict(deployments)
        else:
            self.deployments = parse_deployment_map(
                deployments if deployments is not None else os.environ.get("AZURE_OPENAI_DEPLOYMENTS")
            )
        if not self.deployments:
            raise LLMError(
                "AZURE_OPENAI_DEPLOYMENTS is empty. Azure routes on deployment name, not "
                "model name, so at least one is required. Example:\n"
                "  AZURE_OPENAI_DEPLOYMENTS=claude-opus-5=my-gpt5-deployment\n"
                "A bare name (AZURE_OPENAI_DEPLOYMENTS=my-deployment) maps to itself and "
                "becomes the default for every prompt template."
            )
        self.default_deployment = next(iter(self.deployments.values()))
        self._caps: dict[str, DeploymentCapabilities] = {}
        self._seed_capability_overrides()

        if client is not None:
            self._client = client
            return
        try:
            from openai import AzureOpenAI
        except ImportError as exc:  # pragma: no cover - environment dependent
            raise LLMError(
                "the `openai` package is required for the Azure backend: pip install openai"
            ) from exc

        endpoint = endpoint or os.environ.get("AZURE_OPENAI_ENDPOINT")
        api_key = api_key or os.environ.get("AZURE_OPENAI_API_KEY")
        api_version = api_version or os.environ.get("AZURE_OPENAI_API_VERSION")
        missing = [
            n
            for n, v in (
                ("AZURE_OPENAI_ENDPOINT", endpoint),
                ("AZURE_OPENAI_API_KEY", api_key),
                ("AZURE_OPENAI_API_VERSION", api_version),
            )
            if not v
        ]
        if missing:
            raise LLMError(
                f"missing Azure configuration: {', '.join(missing)}. Set them in the "
                "environment or in a .env file at the repository root."
            )
        self._client = AzureOpenAI(
            azure_endpoint=endpoint,
            api_key=api_key,
            api_version=api_version,
            timeout=timeout,
            max_retries=4,
        )
        log.info(
            "Azure OpenAI ready | endpoint=%s api-version=%s deployments=%s",
            endpoint,
            api_version,
            ", ".join(f"{k}->{v}" for k, v in self.deployments.items()),
        )

    def _seed_capability_overrides(self) -> None:
        """Skip probing for deployments the operator has already classified."""
        reasoning = {
            d.strip()
            for d in os.environ.get("AZURE_OPENAI_REASONING_DEPLOYMENTS", "").split(",")
            if d.strip()
        }
        chat = {
            d.strip()
            for d in os.environ.get("AZURE_OPENAI_CHAT_DEPLOYMENTS", "").split(",")
            if d.strip()
        }
        for d in reasoning:
            self._caps[d] = DeploymentCapabilities(
                token_param="max_completion_tokens", send_reasoning_effort=True
            )
        for d in chat:
            self._caps[d] = DeploymentCapabilities(
                token_param="max_tokens", send_reasoning_effort=False
            )

    def describe(self) -> str:
        return f"{self.name} ({len(self.deployments)} deployment(s))"

    def resolve_model(self, requested: str) -> str:
        """Map a template's logical model to a deployment name.

        An unmapped model falls back to the default deployment rather than failing: prompt
        templates name Anthropic models, and requiring an edit to every template to run on
        Azure would defeat the point of a shared registry.
        """
        deployment = self.deployments.get(requested)
        if deployment is None:
            deployment = self.default_deployment
            log.info(
                "model %r has no entry in AZURE_OPENAI_DEPLOYMENTS; using default deployment "
                "%r. Add '%s=<deployment>' to pin it.",
                requested,
                deployment,
                requested,
            )
        return deployment

    def complete(self, prompt: RenderedPrompt, schema: dict | None = None) -> BackendResponse:
        import openai

        deployment = self.resolve_model(prompt.model)
        caps = self._caps.setdefault(deployment, DeploymentCapabilities())

        last_error: Exception | None = None
        for attempt in range(4):
            kwargs = self._build_kwargs(prompt, deployment, caps, schema)
            try:
                response = self._client.chat.completions.create(**kwargs)
                break
            except openai.BadRequestError as exc:
                if self._is_content_filter(exc):
                    raise LLMError(
                        f"{prompt.template_id}: blocked by the Azure content filter before "
                        f"generation. This is a policy decision by the deployment's filter "
                        f"configuration, not an empty result -- do not treat it as one. "
                        f"Detail: {exc}"
                    ) from exc
                adjusted = self._adapt(caps, str(exc))
                if not adjusted:
                    raise LLMError(f"{prompt.template_id}: Azure rejected the request: {exc}") from exc
                log.info(
                    "deployment %s: %s; retrying with %s", deployment, adjusted, caps.describe()
                )
                last_error = exc
            except openai.NotFoundError as exc:
                raise LLMError(
                    f"deployment {deployment!r} not found at this endpoint. Azure routes on "
                    f"deployment name, not model name -- check AZURE_OPENAI_DEPLOYMENTS "
                    f"against the deployments in your Azure OpenAI resource. Detail: {exc}"
                ) from exc
            except openai.AuthenticationError as exc:
                raise LLMError(f"Azure authentication failed: check AZURE_OPENAI_API_KEY: {exc}") from exc
            except openai.PermissionDeniedError as exc:
                raise LLMError(f"Azure denied access to {deployment!r}: {exc}") from exc
            except openai.RateLimitError as exc:
                raise LLMError(
                    f"Azure rate limited after SDK retries (deployment {deployment!r} may be "
                    f"under-provisioned in TPM): {exc}"
                ) from exc
            except openai.APIConnectionError as exc:
                raise LLMError(f"network failure reaching {deployment!r}: {exc}") from exc
        else:
            raise LLMError(
                f"{prompt.template_id}: could not find a parameter set {deployment!r} accepts "
                f"after 4 attempts. Last error: {last_error}"
            )

        return self._normalize(prompt, deployment, caps, response, schema)

    # -- request construction -------------------------------------------------

    def _build_kwargs(
        self,
        prompt: RenderedPrompt,
        deployment: str,
        caps: DeploymentCapabilities,
        schema: dict | None,
    ) -> dict[str, Any]:
        # System anchor first and byte-identical every call: Azure's prompt caching is
        # automatic and prefix-based, so a stable prefix is what earns the cache hit.
        kwargs: dict[str, Any] = {
            "model": deployment,
            "messages": [
                {"role": "system", "content": prompt.system},
                {"role": "user", "content": prompt.user},
            ],
            caps.token_param: prompt.max_tokens,
        }
        if caps.send_reasoning_effort:
            kwargs["reasoning_effort"] = _EFFORT_TO_AZURE.get(prompt.effort, "high")
        if schema is not None:
            if caps.structured_output == "json_schema":
                kwargs["response_format"] = {
                    "type": "json_schema",
                    "json_schema": {
                        "name": prompt.template_id.replace(".", "_"),
                        "strict": True,
                        "schema": schema,
                    },
                }
            elif caps.structured_output == "json_object":
                kwargs["response_format"] = {"type": "json_object"}
                kwargs["messages"][1] = {
                    "role": "user",
                    "content": (
                        f"{prompt.user}\n\nReturn a single JSON object conforming to this "
                        f"schema:\n{json.dumps(schema, indent=2)}"
                    ),
                }
        return kwargs

    def _adapt(self, caps: DeploymentCapabilities, message: str) -> str | None:
        """Adjust capabilities from a 400 message. Returns what changed, or None."""
        lowered = message.lower()
        unsupported = bool(_UNSUPPORTED.search(lowered))

        if "max_completion_tokens" in lowered and caps.token_param == "max_completion_tokens":
            caps.token_param = "max_tokens"
            return "deployment wants max_tokens"
        if "max_tokens" in lowered and caps.token_param == "max_tokens":
            caps.token_param = "max_completion_tokens"
            return "deployment wants max_completion_tokens"
        if "reasoning_effort" in lowered and caps.send_reasoning_effort:
            caps.send_reasoning_effort = False
            return "deployment does not accept reasoning_effort"
        if ("json_schema" in lowered or "response_format" in lowered) and unsupported:
            if caps.structured_output == "json_schema":
                caps.structured_output = "json_object"
                return "deployment or api-version does not support strict json_schema"
            if caps.structured_output == "json_object":
                caps.structured_output = "prompt_only"
                return "deployment does not support response_format at all"
        return None

    @staticmethod
    def _is_content_filter(exc: Exception) -> bool:
        body = getattr(exc, "body", None)
        if isinstance(body, dict):
            code = str(body.get("error", {}).get("code", "") or body.get("code", ""))
            if "content_filter" in code:
                return True
        return "content_filter" in str(exc).lower() or "responsibleai" in str(exc).lower()

    # -- response normalization -----------------------------------------------

    def _normalize(
        self,
        prompt: RenderedPrompt,
        deployment: str,
        caps: DeploymentCapabilities,
        response: Any,
        schema: dict | None,
    ) -> BackendResponse:
        if not response.choices:
            raise LLMError(
                f"{prompt.template_id}: Azure returned no choices. This usually means the "
                "prompt was filtered at the input stage; check the deployment's content "
                "filter configuration."
            )
        choice = response.choices[0]
        finish = choice.finish_reason

        if finish == "content_filter":
            raise LLMError(
                f"{prompt.template_id}: generation stopped by the Azure content filter. The "
                "partial text is not a usable result and is not returned."
            )
        if finish == "length":
            raise LLMError(
                f"{prompt.template_id}: response hit the token limit ({prompt.max_tokens}) and "
                "is truncated. Raise max_tokens in the template rather than using partial "
                "output. On a reasoning deployment the limit covers reasoning tokens too, so "
                "it needs more headroom than the visible answer suggests."
            )

        text = choice.message.content or ""
        if not text.strip():
            raise LLMError(
                f"{prompt.template_id}: empty response from {deployment!r} "
                f"(finish_reason={finish!r}). On a reasoning deployment this usually means "
                "the token budget was consumed by reasoning before any answer was emitted."
            )

        usage = response.usage
        prompt_details = getattr(usage, "prompt_tokens_details", None)
        cached = getattr(prompt_details, "cached_tokens", 0) or 0
        completion_details = getattr(usage, "completion_tokens_details", None)
        reasoning_tokens = getattr(completion_details, "reasoning_tokens", 0) or 0

        return BackendResponse(
            text=text,
            parsed=_parse_json(text, prompt.template_id) if schema is not None else None,
            model=deployment,
            stop_reason=finish,
            input_tokens=(getattr(usage, "prompt_tokens", 0) or 0) - cached,
            output_tokens=getattr(usage, "completion_tokens", 0) or 0,
            cached_input_tokens=cached,
            cache_write_tokens=0,  # Azure caching is implicit; there is no write charge.
            # Azure pricing varies by region, tier, and agreement. Reporting a first-party
            # list price here would be wrong, so tokens are recorded and cost is left unset.
            cost_usd=None,
            request_id=getattr(response, "id", None),
            provider=self.name,
            extra={
                "deployment": deployment,
                "capabilities": caps.describe(),
                "reasoning_tokens": reasoning_tokens,
                "api_version": os.environ.get("AZURE_OPENAI_API_VERSION", "unset"),
            },
        )


# ------------------------------------------------------------------------ shared


def _parse_json(text: str, template_id: str) -> Any:
    try:
        return json.loads(text)
    except json.JSONDecodeError as exc:
        raise LLMError(
            f"{template_id}: structured output is not valid JSON. First 200 chars: "
            f"{text[:200]!r}"
        ) from exc


def resolve_backend(name: str | None = None, **kwargs: Any) -> Backend:
    """Construct a backend by name, or auto-detect from the environment.

    Explicit ``name`` (or ``PANSPATIAL_BACKEND``) wins. Otherwise Azure is selected when its
    variables are present, so a configured ``.env`` is enough to switch providers without
    touching a prompt template or a command line.
    """
    requested = (name or os.environ.get("PANSPATIAL_BACKEND") or "").strip().lower()
    if requested in {"anthropic", "claude"}:
        return AnthropicBackend(**kwargs)
    if requested in {"azure", "azure-openai", "azure_openai"}:
        return AzureOpenAIBackend(**kwargs)
    if requested:
        raise LLMError(f"unknown backend {requested!r}; choose 'anthropic' or 'azure'")

    if os.environ.get("AZURE_OPENAI_API_KEY") and os.environ.get("AZURE_OPENAI_ENDPOINT"):
        log.info("auto-selected the Azure OpenAI backend from the environment")
        return AzureOpenAIBackend(**kwargs)
    return AnthropicBackend(**kwargs)
