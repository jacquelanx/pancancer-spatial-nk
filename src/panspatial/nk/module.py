"""The NK spatial niche analysis, exposed through the common module contract.

The NK layer predates the module framework and keeps its own tested interface in
:mod:`panspatial.nk.nk_spatial`. This wrapper registers it alongside the other analysis
layers so it appears in the registry, carries the same provenance, and is bound by the same
platform gate.
"""

from __future__ import annotations

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
from panspatial.nk.nk_spatial import (
    cohort_meta_analysis,
    per_section_niche_test,
    require_single_cell_resolution,
)


@register
class NKSpatialNiche(AnalysisModule):
    name = "nk_spatial"
    version = "0.2.0"
    phase = "4"
    description = "NK identity, states, and torus-shift spatial niche association."
    evidence = Evidence.MEASURED
    replication_unit = ReplicationUnit.PATIENT
    # Single-cell resolution only: an NK weight on a 55 um spot is a similarity score.
    supported_platforms = frozenset({PlatformClass.SINGLE_CELL})
    optional_dependencies = ("scanpy", "squidpy", "scikit-learn")

    def run(
        self,
        ctx: ModuleContext,
        *,
        coords: Any = None,
        labels: Sequence[str] | None = None,
        sections: Sequence[str] | None = None,
        patients: Sequence[str] | None = None,
        anchor: str | None = None,
        **inputs: Any,
    ) -> ModuleResult:
        self.check_platform(ctx)
        require_single_cell_resolution(ctx.platform or "unknown")
        if coords is None or labels is None or sections is None or patients is None:
            raise ModuleError(
                "nk_spatial.run requires coords=, labels=, sections= and patients=. "
                "Coordinates must already be in microns."
            )
        anchor = anchor or ctx.param("anchor")
        if not anchor:
            raise ModuleError("nk_spatial.run requires anchor= (the reference structure)")

        frame = per_section_niche_test(
            np.asarray(coords, dtype=float), labels, sections, patients,
            focal_label=ctx.param("focal", "NK"), anchor_label=anchor,
            n_perms=ctx.param("n_perms", 999),
            alternative=ctx.param("alternative", "less"),
            min_focal=ctx.param("min_focal", 30), seed=ctx.seed,
        )
        notes = [
            "Null model is a torus shift within the tissue mask, preserving spatial "
            "autocorrelation; a label shuffle would inflate significance.",
            f"{len(frame.attrs.get('skipped_sections', []))} section(s) excluded for too "
            "few NK cells and counted as excluded.",
        ]
        meta = None
        if len(frame) >= 2:
            meta = cohort_meta_analysis(frame)
            notes.append(str(meta))
        return self.result(
            ctx, frame, notes=notes, anchor=anchor,
            n_sections=int(len(frame)),
            pooled_estimate=None if meta is None else meta.estimate,
            pooled_p_value=None if meta is None else meta.p_value,
            n_patients=None if meta is None else meta.k,
        )
