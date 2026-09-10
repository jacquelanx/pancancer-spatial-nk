"""Fate and trajectory inference, with the corroboration rule velocity requires.

Two assumptions underpinning RNA velocity are readily violated: a common splicing rate
across genes, and that the system's equilibria are observed during the experiment. Where
they fail, inference yields incorrect results. Velocity methods also fit genes
independently and so do not preserve gene–gene coherence — the same cell can sit at the
start of one gene's trajectory and the end of another's — and in systems dominated by
post-transcriptional regulation, splicing kinetics need not reflect cellular dynamics at
all. Tumour tissue is such a system.

CellRank 2 is the primary framework precisely because it unifies several views and does not
require velocity. This module enforces the consequence: a direction of change is reported
only when at least two independent views agree.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

import numpy as np

from panspatial.modules.base import (
    AnalysisModule,
    Evidence,
    ModuleContext,
    ModuleError,
    ModuleResult,
    PlatformClass,
    ReplicationUnit,
    register,
)

log = logging.getLogger("panspatial.dynamics")

# Views that can drive a CellRank kernel. Velocity is deliberately not independent of
# itself: two velocity flavours agreeing is one view, not two.
VIEW_FAMILIES = {
    "velocity": ("scvelo", "unitvelo", "velovi", "dynamo"),
    "pseudotime": ("dpt", "palantir", "monocle3", "slingshot", "paga"),
    "realtime": ("timepoint", "moscot"),
    "similarity": ("connectivity", "cytotrace"),
    "spatial": ("spatial_gradient", "veloagent"),
}


def view_family(view: str) -> str:
    for family, members in VIEW_FAMILIES.items():
        if view.lower() in members or view.lower() == family:
            return family
    log.warning("unrecognised view %r; treating it as its own family", view)
    return view.lower()


@dataclass
class Corroboration:
    transition: str
    views: dict[str, str]                  # view -> direction ("forward"/"reverse"/"none")
    families_agreeing: list[str] = field(default_factory=list)
    direction: str | None = None
    corroborated: bool = False
    reason: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "transition": self.transition,
            "direction": self.direction,
            "corroborated": self.corroborated,
            "n_independent_families": len(self.families_agreeing),
            "families": ",".join(sorted(self.families_agreeing)),
            "views": ";".join(f"{k}={v}" for k, v in sorted(self.views.items())),
            "reason": self.reason,
        }


def corroborate_direction(
    transition: str,
    views: Mapping[str, str],
    *,
    min_independent_families: int = 2,
) -> Corroboration:
    """Report a direction only when independent method families agree on it.

    Two velocity implementations agreeing is one line of evidence, because they share the
    assumptions that fail. ``min_independent_families`` counts *families*, not tools.
    """
    if not views:
        raise ValueError("at least one view is required")
    normalised = {k: str(v).lower() for k, v in views.items()}
    bad = {k: v for k, v in normalised.items() if v not in {"forward", "reverse", "none"}}
    if bad:
        raise ValueError(f"directions must be forward|reverse|none, got {bad}")

    by_direction: dict[str, set[str]] = {}
    for view, direction in normalised.items():
        if direction == "none":
            continue
        by_direction.setdefault(direction, set()).add(view_family(view))

    if not by_direction:
        return Corroboration(transition, normalised, [], None, False,
                             "no view reported a direction")

    best = max(by_direction, key=lambda d: len(by_direction[d]))
    families = sorted(by_direction[best])
    conflicting = {d: sorted(f) for d, f in by_direction.items() if d != best}

    if conflicting:
        return Corroboration(
            transition, normalised, families, None, False,
            f"views disagree: {best} supported by {families}, but "
            + "; ".join(f"{d} supported by {f}" for d, f in conflicting.items())
            + ". A contested direction is not reported.",
        )
    if len(families) < min_independent_families:
        return Corroboration(
            transition, normalised, families, None, False,
            f"only {len(families)} independent method family ({families[0]}) supports "
            f"{best}; {min_independent_families} are required. Two velocity tools agreeing "
            "is one view, because they share the assumptions that fail in tumour tissue.",
        )
    return Corroboration(
        transition, normalised, families, best, True,
        f"{best} supported by {len(families)} independent families: {', '.join(families)}",
    )


def corroboration_table(results: Sequence[Corroboration]) -> Any:
    import pandas as pd
    return pd.DataFrame([r.to_dict() for r in results])


@register
class Dynamics(AnalysisModule):
    name = "dynamics"
    version = "0.1.0"
    phase = "5"
    description = "Trajectory and fate mapping with a multi-view corroboration requirement."
    evidence = Evidence.INFERRED
    replication_unit = ReplicationUnit.CELL
    supported_platforms = frozenset(
        {PlatformClass.SINGLE_CELL, PlatformClass.SPOT, PlatformClass.DISSOCIATED}
    )
    optional_dependencies = ("cellrank", "scvelo", "scFates")

    def run(
        self, ctx: ModuleContext,
        *, transitions: Mapping[str, Mapping[str, str]] | None = None, **inputs: Any,
    ) -> ModuleResult:
        if not transitions:
            raise ModuleError(
                "dynamics.run requires transitions={transition: {view: direction}}. Run "
                "CellRank 2 kernels (pseudotime, real-time, similarity, and velocity where "
                "defensible) first, then pass their directions here for corroboration."
            )
        results = [
            corroborate_direction(
                name, views,
                min_independent_families=ctx.param("min_independent_families", 2),
            )
            for name, views in transitions.items()
        ]
        frame = corroboration_table(results)
        confirmed = frozenset(r.transition for r in results if r.corroborated)
        return self.result(
            ctx, frame, reliable_entities=confirmed,
            notes=[
                "A direction of change is reported only where at least two independent "
                "method families agree; velocity alone is never sufficient.",
                f"{len(confirmed)}/{len(results)} transitions corroborated.",
            ],
        )
