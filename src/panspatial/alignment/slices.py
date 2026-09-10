"""Multi-section alignment into a common frame.

A Genome Biology 2024 benchmark found SPACEL and PASTE gave the highest layer-wise
alignment accuracy, with STAligner best among the embedding-based approaches. PASTE2 solves
a partial fused Gromov–Wasserstein problem for slices that only partly overlap and
estimates the overlap fraction rather than assuming it; STalign uses diffeomorphic metric
mapping to absorb local non-linear distortion; SLAT is graph-based and aligns across
distinct technologies and modalities.

The choice among them is driven by two facts about the sections, not by preference:
whether they overlap fully, and whether they come from the same technology.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Sequence

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

log = logging.getLogger("panspatial.alignment")


@dataclass(frozen=True)
class AlignerChoice:
    method: str
    rationale: str
    estimate_overlap: bool

    def to_dict(self) -> dict[str, Any]:
        return {"method": self.method, "rationale": self.rationale,
                "estimate_overlap": self.estimate_overlap}


def select_aligner(
    *, same_technology: bool, full_overlap: bool, nonlinear_distortion: bool = False
) -> AlignerChoice:
    """Pick an alignment method from the properties of the sections being aligned."""
    if not same_technology:
        return AlignerChoice(
            "slat",
            "Sections come from different technologies or modalities; SLAT is the "
            "graph-based method designed for heterogeneous alignment.",
            estimate_overlap=False,
        )
    if nonlinear_distortion:
        return AlignerChoice(
            "stalign",
            "Local non-linear distortion is expected; STalign's diffeomorphic metric "
            "mapping models it rather than approximating it with a rigid transform.",
            estimate_overlap=False,
        )
    if full_overlap:
        return AlignerChoice(
            "paste",
            "Same technology, fully overlapping sections: PASTE was among the most "
            "accurate for layer-wise alignment in benchmark.",
            estimate_overlap=False,
        )
    return AlignerChoice(
        "paste2",
        "Sections overlap only partially; PASTE2 solves the partial fused "
        "Gromov-Wasserstein problem and estimates the overlap fraction rather than "
        "assuming complete correspondence.",
        estimate_overlap=True,
    )


def check_alignment(
    source: np.ndarray, target: np.ndarray, mapping: Sequence[int], *,
    max_residual_um: float = 100.0,
) -> dict[str, Any]:
    """Residual statistics for a computed alignment, so a bad one is visible.

    An alignment is accepted silently far too often. Median residual distance and the
    fraction of points beyond ``max_residual_um`` make failure legible.
    """
    source = np.asarray(source, dtype=float)
    target = np.asarray(target, dtype=float)
    mapping = np.asarray(mapping, dtype=int)
    if len(mapping) != len(source):
        raise ValueError("mapping must give one target index per source point")
    if mapping.min() < 0 or mapping.max() >= len(target):
        raise ValueError("mapping contains out-of-range target indices")

    residual = np.linalg.norm(source - target[mapping], axis=1)
    out = {
        "n_points": int(len(source)),
        "median_residual_um": float(np.median(residual)),
        "p90_residual_um": float(np.percentile(residual, 90)),
        "fraction_beyond_threshold": float(np.mean(residual > max_residual_um)),
        "acceptable": bool(np.median(residual) <= max_residual_um),
    }
    if not out["acceptable"]:
        log.warning(
            "alignment median residual %.1f um exceeds %.1f um; do not build a 3D "
            "reconstruction on it", out["median_residual_um"], max_residual_um,
        )
    return out


@register
class SliceAlignment(AnalysisModule):
    name = "alignment"
    version = "0.1.0"
    phase = "5"
    description = "Multi-section alignment with property-driven method choice and residual QC."
    evidence = Evidence.INFERRED
    replication_unit = ReplicationUnit.SECTION
    supported_platforms = frozenset({PlatformClass.SINGLE_CELL, PlatformClass.SPOT})
    optional_dependencies = ("paste-bio", "STalign", "scSLAT")

    def run(
        self, ctx: ModuleContext, *, source: Any = None, target: Any = None,
        mapping: Sequence[int] | None = None, **inputs: Any,
    ) -> ModuleResult:
        import pandas as pd

        choice = select_aligner(
            same_technology=ctx.param("same_technology", True),
            full_overlap=ctx.param("full_overlap", True),
            nonlinear_distortion=ctx.param("nonlinear_distortion", False),
        )
        if source is None or target is None or mapping is None:
            raise ModuleError(
                f"alignment.run requires source=, target= and mapping=. Run {choice.method} "
                f"first, then pass the correspondence here for residual QC. "
                f"Rationale: {choice.rationale}"
            )
        stats = check_alignment(
            source, target, mapping,
            max_residual_um=ctx.param("max_residual_um", 100.0),
        )
        notes = [choice.rationale]
        if not stats["acceptable"]:
            notes.append("Alignment failed residual QC; downstream 3D claims are blocked.")
        return self.result(
            ctx, pd.DataFrame([{**stats, **choice.to_dict()}]), notes=notes,
            aligner=choice.method, acceptable=stats["acceptable"],
        )
