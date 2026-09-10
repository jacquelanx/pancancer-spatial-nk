"""The contract every analysis module implements.

The deck specifies sixteen analysis layers. Implemented as sixteen independent scripts
they would drift: different replication units, different evidence standards, different
provenance. This module makes the shared rules structural rather than advisory.

Every module declares four things and cannot produce a result without them:

* **Evidence level** — :class:`Evidence`. Whether a number was measured, inferred from a
  model, or imputed from another modality. An imputed result may generate hypotheses; it
  may never validate the model that produced it, and :func:`assert_can_validate` enforces
  that.
* **Replication unit** — :class:`ReplicationUnit`. Anything reported at cell or spot level
  is marked as such, so the aggregation layer knows it must be pooled before a
  cross-patient claim is made.
* **Platform support** — which platform classes the result is valid on. A rare-cell claim
  from a 55 µm spot is refused at the point of production, not caught in review.
* **Provenance** — parameters, seed, and the versions of every external tool that ran.

Modules return tidy per-sample frames, which is what lets
:mod:`panspatial.stats.spatial_stats` pool them to patients uniformly.
"""

from __future__ import annotations

import hashlib
import json
import logging
import platform
import sys
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from typing import Any, Callable, Iterable, Sequence

log = logging.getLogger("panspatial.modules")


class Evidence(str, Enum):
    """How strongly a result is tied to a measurement.

    The ordering is meaningful: ``MEASURED > INFERRED > IMPUTED``.
    """

    MEASURED = "measured"    # read out directly (segmented cell counts, expression)
    INFERRED = "inferred"    # model output over measured input (deconvolution, CNV, GRN)
    IMPUTED = "imputed"      # predicted from a different modality (expression from H&E)

    @property
    def rank(self) -> int:
        return {"measured": 2, "inferred": 1, "imputed": 0}[self.value]


class ReplicationUnit(str, Enum):
    """The unit an n refers to. Only PATIENT may be quoted in a cross-cohort claim."""

    CELL = "cell"
    SPOT = "spot"
    SECTION = "section"
    PATIENT = "patient"
    COHORT = "cohort"

    @property
    def is_independent(self) -> bool:
        return self in (ReplicationUnit.PATIENT, ReplicationUnit.COHORT)


class PlatformClass(str, Enum):
    SINGLE_CELL = "single-cell"   # Xenium, CosMx, MERFISH — segmented cells
    SPOT = "spot"                 # Visium, Slide-seq — mixtures
    REGION = "region"             # GeoMx — pooled ROIs
    DISSOCIATED = "dissociated"   # scRNA-seq, no coordinates


SINGLE_CELL_PLATFORMS = {"Xenium", "CosMx", "MERFISH", "seqFISH", "STARmap"}
SPOT_PLATFORMS = {"Visium", "Visium HD", "Slide-seq", "Stereo-seq", "ST (legacy)"}
REGION_PLATFORMS = {"GeoMx"}


def platform_class(platform_name: str | None) -> PlatformClass:
    """Classify a platform name into the class that governs what it can support."""
    if platform_name is None:
        return PlatformClass.DISSOCIATED
    if platform_name in SINGLE_CELL_PLATFORMS:
        return PlatformClass.SINGLE_CELL
    if platform_name in SPOT_PLATFORMS:
        return PlatformClass.SPOT
    if platform_name in REGION_PLATFORMS:
        return PlatformClass.REGION
    log.warning("unknown platform %r; treating as dissociated (no spatial claims)", platform_name)
    return PlatformClass.DISSOCIATED


class ModuleError(RuntimeError):
    """A module could not run, or was asked to produce a result it cannot support."""


class DependencyMissing(ModuleError):
    """An optional analysis dependency is not installed."""


def require(package: str, *, extra: str, purpose: str) -> Any:
    """Import an optional dependency or fail with an instruction, not a traceback."""
    import importlib

    try:
        return importlib.import_module(package)
    except ImportError as exc:
        raise DependencyMissing(
            f"{purpose} needs `{package}`, which is not installed. "
            f"Install it with: pip install 'panspatial[{extra}]'"
        ) from exc


# --------------------------------------------------------------------------- results


@dataclass
class Provenance:
    """Everything needed to reproduce one module run."""

    module: str
    module_version: str
    seed: int
    params: dict[str, Any]
    tool_versions: dict[str, str] = field(default_factory=dict)
    created_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    python: str = field(default_factory=lambda: sys.version.split()[0])
    host: str = field(default_factory=platform.node)

    @property
    def fingerprint(self) -> str:
        payload = json.dumps(
            {
                "module": self.module,
                "version": self.module_version,
                "seed": self.seed,
                "params": self.params,
                "tools": self.tool_versions,
            },
            sort_keys=True,
            default=str,
        )
        return hashlib.sha256(payload.encode()).hexdigest()[:16]

    def to_dict(self) -> dict[str, Any]:
        return {
            "module": self.module,
            "module_version": self.module_version,
            "fingerprint": self.fingerprint,
            "seed": self.seed,
            "params": self.params,
            "tool_versions": self.tool_versions,
            "created_at": self.created_at,
            "python": self.python,
            "host": self.host,
        }


@dataclass
class ModuleResult:
    """Tidy output of one module run, carrying the claims it is allowed to support."""

    module: str
    table: Any                                   # pandas DataFrame, one row per unit
    evidence: Evidence
    replication_unit: ReplicationUnit
    supported_platforms: frozenset[PlatformClass]
    provenance: Provenance
    reliable_entities: frozenset[str] | None = None   # e.g. cell types that passed a gate
    notes: list[str] = field(default_factory=list)
    metadata: dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.evidence is Evidence.IMPUTED and not self.notes:
            self.notes.append(
                "Imputed from another modality. Usable for hypothesis generation and cohort "
                "extension; not usable as validation of the model that produced it."
            )

    @property
    def n_rows(self) -> int:
        return 0 if self.table is None else len(self.table)

    def supports(self, platform_name: str | None) -> bool:
        return platform_class(platform_name) in self.supported_platforms

    def assert_supports(self, platform_name: str | None, claim: str) -> None:
        if not self.supports(platform_name):
            raise ModuleError(
                f"{self.module}: {claim} is not supported on {platform_name} "
                f"({platform_class(platform_name).value}). This module's results are valid on: "
                f"{', '.join(sorted(p.value for p in self.supported_platforms))}."
            )

    def assert_entity_reliable(self, entity: str) -> None:
        """Refuse a claim about an entity that failed the module's own reliability gate."""
        if self.reliable_entities is None:
            return
        if entity not in self.reliable_entities:
            raise ModuleError(
                f"{self.module}: {entity!r} did not pass this module's reliability gate and "
                f"may not carry a conclusion. Passed: "
                f"{', '.join(sorted(self.reliable_entities)) or '(none)'}."
            )

    def summary(self) -> dict[str, Any]:
        return {
            "module": self.module,
            "rows": self.n_rows,
            "evidence": self.evidence.value,
            "replication_unit": self.replication_unit.value,
            "platforms": sorted(p.value for p in self.supported_platforms),
            "reliable_entities": sorted(self.reliable_entities) if self.reliable_entities else None,
            "notes": self.notes,
            "provenance": self.provenance.to_dict(),
            **self.metadata,
        }


def assert_can_validate(claim_result: ModuleResult, validating_result: ModuleResult) -> None:
    """Refuse a validation whose evidence is weaker than the claim it validates.

    This is the guard against the circularity the deck warns about: expression predicted
    from histology by a model trained on spatial data cannot validate a spatial finding.
    """
    # The imputed rule is checked first because it is absolute: it holds regardless of
    # what the imputed result is being compared against, and it gives the clearer message.
    if validating_result.evidence is Evidence.IMPUTED:
        raise ModuleError(
            f"{validating_result.module} is imputed from another modality and cannot validate "
            "anything, regardless of what it is compared against."
        )
    if validating_result.evidence.rank < claim_result.evidence.rank:
        raise ModuleError(
            f"{validating_result.module} ({validating_result.evidence.value}) cannot validate "
            f"{claim_result.module} ({claim_result.evidence.value}): the validating evidence is "
            "weaker than the claim. Validation must come from an independent measurement."
        )


# --------------------------------------------------------------------------- modules


@dataclass
class ModuleContext:
    """Inputs shared by every module in a run."""

    sample_id: str
    platform: str | None = None
    patient_id: str | None = None
    cancer_type: str | None = None
    seed: int = 42
    params: dict[str, Any] = field(default_factory=dict)
    paths: dict[str, str] = field(default_factory=dict)

    @property
    def platform_class(self) -> PlatformClass:
        return platform_class(self.platform)

    def param(self, key: str, default: Any = None) -> Any:
        return self.params.get(key, default)


class AnalysisModule(ABC):
    """Base class for every analysis layer in the pipeline."""

    name: str = "unnamed"
    version: str = "0.1.0"
    phase: str = "-"
    description: str = ""
    evidence: Evidence = Evidence.INFERRED
    replication_unit: ReplicationUnit = ReplicationUnit.SECTION
    supported_platforms: frozenset[PlatformClass] = frozenset(
        {PlatformClass.SINGLE_CELL, PlatformClass.SPOT}
    )
    optional_dependencies: tuple[str, ...] = ()

    def check_platform(self, ctx: ModuleContext) -> None:
        if ctx.platform_class not in self.supported_platforms:
            raise ModuleError(
                f"{self.name} does not support {ctx.platform!r} "
                f"({ctx.platform_class.value}). Supported: "
                f"{', '.join(sorted(p.value for p in self.supported_platforms))}."
            )

    def provenance(self, ctx: ModuleContext, tool_versions: dict[str, str] | None = None) -> Provenance:
        return Provenance(
            module=self.name,
            module_version=self.version,
            seed=ctx.seed,
            params=dict(ctx.params),
            tool_versions=tool_versions or {},
        )

    def result(
        self,
        ctx: ModuleContext,
        table: Any,
        *,
        evidence: Evidence | None = None,
        replication_unit: ReplicationUnit | None = None,
        reliable_entities: Iterable[str] | None = None,
        notes: Sequence[str] = (),
        tool_versions: dict[str, str] | None = None,
        **metadata: Any,
    ) -> ModuleResult:
        return ModuleResult(
            module=self.name,
            table=table,
            evidence=evidence or self.evidence,
            replication_unit=replication_unit or self.replication_unit,
            supported_platforms=self.supported_platforms,
            provenance=self.provenance(ctx, tool_versions),
            reliable_entities=frozenset(reliable_entities) if reliable_entities is not None else None,
            notes=list(notes),
            metadata=metadata,
        )

    @abstractmethod
    def run(self, ctx: ModuleContext, **inputs: Any) -> ModuleResult:
        """Execute the module and return a tidy, provenance-carrying result."""

    def __repr__(self) -> str:  # pragma: no cover - display only
        return f"<{type(self).__name__} {self.name} v{self.version} phase {self.phase}>"


# -------------------------------------------------------------------------- registry

_REGISTRY: dict[str, type[AnalysisModule]] = {}


def register(cls: type[AnalysisModule]) -> type[AnalysisModule]:
    """Class decorator adding a module to the global registry."""
    if cls.name in _REGISTRY and _REGISTRY[cls.name] is not cls:
        raise ModuleError(f"duplicate module name {cls.name!r}")
    _REGISTRY[cls.name] = cls
    return cls


def get_module(name: str) -> AnalysisModule:
    try:
        return _REGISTRY[name]()
    except KeyError:
        raise ModuleError(
            f"unknown module {name!r}; registered: {', '.join(sorted(_REGISTRY))}"
        ) from None


def registered_modules() -> dict[str, type[AnalysisModule]]:
    return dict(sorted(_REGISTRY.items()))


def modules_for_phase(phase: str) -> list[type[AnalysisModule]]:
    return [c for c in registered_modules().values() if c.phase == str(phase)]


def describe_registry() -> Any:
    """A DataFrame describing every registered module — the pipeline's table of contents."""
    import pandas as pd

    return pd.DataFrame(
        [
            {
                "module": c.name,
                "version": c.version,
                "phase": c.phase,
                "evidence": c.evidence.value,
                "replication_unit": c.replication_unit.value,
                "platforms": ",".join(sorted(p.value for p in c.supported_platforms)),
                "optional_deps": ",".join(c.optional_dependencies),
                "description": c.description,
            }
            for c in registered_modules().values()
        ]
    )
