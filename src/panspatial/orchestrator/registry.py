"""Versioned prompt registry with strict variable binding and provenance hashing.

Every LLM call in this pipeline is driven by a file under ``prompts/``, never by an
ad-hoc string. Three properties matter:

1. **Strict binding.** A template with an unbound ``{{variable}}`` raises rather than
   sending a literal ``{{variable}}`` to the model. An unnoticed placeholder in a
   manuscript prompt is how a draft ends up containing invented numbers.
2. **Provenance.** The SHA-256 of (system prompt + rendered prompt + model + effort)
   identifies the exact call. Run records key on it, so any sentence in the manuscript
   can be traced back to the prompt and inputs that produced it.
3. **Gating.** Templates marked ``requires_verified_results`` refuse to render unless the
   caller asserts the input numbers came from an executed analysis.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

_FRONTMATTER = re.compile(r"\A---\n(.*?)\n---\n(.*)\Z", re.DOTALL)
_PLACEHOLDER = re.compile(r"\{\{\s*([A-Za-z_][A-Za-z0-9_]*)\s*\}\}")

SYSTEM_PROMPT_ID = "system.meta"


class PromptError(RuntimeError):
    """Raised when a template cannot be loaded or safely rendered."""


@dataclass(frozen=True)
class VariableSpec:
    name: str
    required: bool = True
    default: str | None = None
    description: str = ""


@dataclass(frozen=True)
class PromptTemplate:
    id: str
    version: str
    path: Path
    body: str
    meta: dict[str, Any]
    variables: dict[str, VariableSpec] = field(default_factory=dict)

    @property
    def phase(self) -> str:
        return str(self.meta.get("phase", "-"))

    @property
    def model(self) -> str:
        return str(self.meta.get("model", "claude-opus-5"))

    @property
    def effort(self) -> str:
        return str(self.meta.get("effort", "high"))

    @property
    def max_tokens(self) -> int:
        return int(self.meta.get("max_tokens", 16000))

    @property
    def structured_output(self) -> str | None:
        value = self.meta.get("structured_output")
        return str(value) if value else None

    @property
    def requires_verified_results(self) -> bool:
        return bool(self.meta.get("requires_verified_results", False))

    def declared_placeholders(self) -> set[str]:
        return set(_PLACEHOLDER.findall(self.body))

    def render(self, values: dict[str, str], *, results_verified: bool = False) -> str:
        """Substitute variables, refusing anything that could silently ship a placeholder."""
        if self.requires_verified_results and not results_verified:
            raise PromptError(
                f"{self.id} writes manuscript text from result values and is gated: pass "
                "results_verified=True only when the inputs came from an executed analysis "
                "(a real results table), not from an example or a guess."
            )

        unknown = set(values) - set(self.variables)
        if unknown:
            raise PromptError(f"{self.id}: undeclared variables supplied: {sorted(unknown)}")

        resolved: dict[str, str] = {}
        missing: list[str] = []
        for name, spec in self.variables.items():
            if name in values and values[name] is not None:
                resolved[name] = str(values[name])
            elif spec.default is not None:
                resolved[name] = spec.default
            elif spec.required:
                missing.append(name)
        if missing:
            detail = "; ".join(f"{n} ({self.variables[n].description})" for n in sorted(missing))
            raise PromptError(f"{self.id}: missing required variables: {detail}")

        rendered = _PLACEHOLDER.sub(lambda m: resolved.get(m.group(1), m.group(0)), self.body)

        leftover = sorted(set(_PLACEHOLDER.findall(rendered)))
        if leftover:
            raise PromptError(
                f"{self.id}: unresolved placeholders after render: {leftover}. "
                "Declare them in the template frontmatter or supply values."
            )
        return rendered


@dataclass(frozen=True)
class RenderedPrompt:
    template_id: str
    template_version: str
    system: str
    user: str
    model: str
    effort: str
    max_tokens: int
    structured_output: str | None
    bound_values: dict[str, str]

    @property
    def fingerprint(self) -> str:
        """Stable identifier for this exact call; goes into every run record."""
        h = hashlib.sha256()
        for part in (
            self.template_id,
            self.template_version,
            self.model,
            self.effort,
            self.system,
            self.user,
        ):
            h.update(part.encode("utf-8"))
            h.update(b"\x00")
        return h.hexdigest()[:16]


def _parse_variables(raw: Any) -> dict[str, VariableSpec]:
    if not raw:
        return {}
    if not isinstance(raw, dict):
        raise PromptError(f"`variables` must be a mapping, got {type(raw).__name__}")
    out: dict[str, VariableSpec] = {}
    for name, spec in raw.items():
        if isinstance(spec, str):
            out[name] = VariableSpec(name=name, description=spec)
        elif isinstance(spec, dict):
            default = spec.get("default")
            out[name] = VariableSpec(
                name=name,
                required=bool(spec.get("required", default is None)),
                default=None if default is None else str(default),
                description=str(spec.get("description", "")),
            )
        else:
            raise PromptError(f"variable {name!r}: expected str or mapping")
    return out


def _load_template(path: Path) -> PromptTemplate:
    text = path.read_text(encoding="utf-8")
    match = _FRONTMATTER.match(text)
    if not match:
        raise PromptError(f"{path}: missing YAML frontmatter delimited by '---'")
    meta = yaml.safe_load(match.group(1)) or {}
    if not isinstance(meta, dict):
        raise PromptError(f"{path}: frontmatter must be a mapping")
    for key in ("id", "version"):
        if key not in meta:
            raise PromptError(f"{path}: frontmatter missing required key {key!r}")
    return PromptTemplate(
        id=str(meta["id"]),
        version=str(meta["version"]),
        path=path,
        body=match.group(2).strip(),
        meta=meta,
        variables=_parse_variables(meta.get("variables")),
    )


class PromptRegistry:
    """Loads every template under ``prompts/`` and renders calls against the system anchor."""

    def __init__(self, root: str | Path = "prompts") -> None:
        self.root = Path(root)
        if not self.root.is_dir():
            raise PromptError(f"prompt root not found: {self.root}")
        self._templates: dict[str, PromptTemplate] = {}
        for path in sorted(self.root.rglob("*.md")):
            template = _load_template(path)
            if template.id in self._templates:
                raise PromptError(
                    f"duplicate prompt id {template.id!r}: "
                    f"{self._templates[template.id].path} and {path}"
                )
            self._templates[template.id] = template
        if SYSTEM_PROMPT_ID not in self._templates:
            raise PromptError(f"system anchor {SYSTEM_PROMPT_ID!r} not found under {self.root}")

    def __len__(self) -> int:
        return len(self._templates)

    def ids(self) -> list[str]:
        return sorted(self._templates)

    def get(self, template_id: str) -> PromptTemplate:
        try:
            return self._templates[template_id]
        except KeyError:
            raise PromptError(
                f"unknown prompt id {template_id!r}; available: {', '.join(self.ids())}"
            ) from None

    @property
    def system_prompt(self) -> str:
        return self._templates[SYSTEM_PROMPT_ID].body

    def render(
        self,
        template_id: str,
        values: dict[str, str] | None = None,
        *,
        results_verified: bool = False,
        model: str | None = None,
        effort: str | None = None,
    ) -> RenderedPrompt:
        template = self.get(template_id)
        if template.meta.get("kind") == "system":
            raise PromptError(f"{template_id} is the system anchor, not a task prompt")
        values = values or {}
        user = template.render(values, results_verified=results_verified)
        return RenderedPrompt(
            template_id=template.id,
            template_version=template.version,
            system=self.system_prompt,
            user=user,
            model=model or template.model,
            effort=effort or template.effort,
            max_tokens=template.max_tokens,
            structured_output=template.structured_output,
            bound_values={k: str(v) for k, v in values.items()},
        )
