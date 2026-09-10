"""Niche calling adapters and the module wrapper.

NicheCompass (Nature Genetics 2025) is the default on segmented platforms: in benchmark it
uniquely recovered spatially contiguous niches and led on spatial consistency and niche
coherence against BANKSY, GraphST and CellCharter. BANKSY, UTAG, CellCharter and scNiche
are retained as scalable alternatives — BANKSY, UTAG, CellCharter and scNiche are the
methods shown to scale beyond three million cells, which matters on Visium HD and large
Xenium runs.
"""

from __future__ import annotations

import logging
from typing import Any, Mapping, Sequence

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
from panspatial.niches.recurrence import composition_matrix, match_niches, recurrence_test

log = logging.getLogger("panspatial.niches")

METHODS = {
    "nichecompass": "Graph deep learning over communication signals; leads on niche coherence",
    "banksy": "Neighbour-augmented features with a tunable microenvironment weight; scales >3M cells",
    "cellcharter": "Neighbourhood aggregation with explicit cluster-number selection",
    "utag": "Linear neighbour weighting; fast baseline that scales >3M cells",
    "scniche": "Multi-view learning over cell neighbourhoods; scales >3M cells",
}
SCALABLE = ("banksy", "utag", "cellcharter", "scniche")


def recommend_method(n_cells: int, *, platform: str | None = None) -> tuple[str, str]:
    """Pick a niche caller from cohort size and platform."""
    if n_cells > 3_000_000:
        return "banksy", (
            f"{n_cells:,} cells exceeds what NicheCompass has been demonstrated on; BANKSY, "
            "UTAG, CellCharter and scNiche are the methods shown to scale past 3M"
        )
    return "nichecompass", (
        "NicheCompass uniquely recovered spatially contiguous niches and led on spatial "
        "consistency and niche coherence against BANKSY, GraphST and CellCharter"
    )


@register
class Niches(AnalysisModule):
    name = "niches"
    version = "0.1.0"
    phase = "4"
    description = "Niche calling, cross-sample matching, and patient-level recurrence testing."
    evidence = Evidence.INFERRED
    replication_unit = ReplicationUnit.PATIENT
    supported_platforms = frozenset({PlatformClass.SINGLE_CELL, PlatformClass.SPOT})
    optional_dependencies = ("nichecompass", "squidpy", "banksy")

    def run(
        self,
        ctx: ModuleContext,
        *,
        sample_compositions: Mapping[str, tuple[Any, Sequence[str]]] | None = None,
        reference: Any = None,
        reference_names: Sequence[str] | None = None,
        patient_of: Mapping[str, str] | None = None,
        cancer_type_of: Mapping[str, str] | None = None,
        **inputs: Any,
    ) -> ModuleResult:
        self.check_platform(ctx)
        if not sample_compositions or reference is None or reference_names is None:
            raise ModuleError(
                "niches.run requires sample_compositions=, reference= and reference_names=. "
                "Call niches per section with the recommended method first, then pass the "
                "per-section compositions here for matching and recurrence testing."
            )
        if not patient_of or not cancer_type_of:
            raise ModuleError(
                "patient_of= and cancer_type_of= are required: recurrence is tested at "
                "patient level, and a section-level count would be pseudoreplication"
            )
        matches = match_niches(
            reference, reference_names, sample_compositions,
            min_similarity=ctx.param("min_similarity", 0.7),
        )
        frame = recurrence_test(
            matches, patient_of, cancer_type_of,
            null_prevalence=ctx.param("null_prevalence", 0.25),
            min_cancer_types=ctx.param("min_cancer_types", 3),
            min_patient_fraction=ctx.param("min_patient_fraction", 0.5),
        )
        conserved = frozenset(frame.loc[frame["conserved"], "niche"]) if len(frame) else frozenset()
        return self.result(
            ctx, frame, reliable_entities=conserved,
            notes=[
                "Recurrence is tested over patients, not sections.",
                "A niche in fewer than the required number of cancer types is reported as "
                "cancer-type-restricted, never as pan-cancer.",
            ],
            n_matched=int(sum(m.matched for m in matches)),
            n_local_niches=len(matches),
        )
