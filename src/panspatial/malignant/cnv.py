"""Choosing a CNV caller from what the data actually supports.

Benchmark findings this module encodes (Brief Bioinform 2025; Nat Commun 2025):

* Numbat, which uses expression *and* B-allele frequency, performed best overall at
  separating tumour from normal — but its accuracy degrades sharply as sequencing depth
  falls, and it is the most CPU-expensive of the tools tested.
* CopyKAT was the best of the expression-only tools and is notably depth-robust: it still
  separated tumour from normal at a median of 1,000 UMIs per cell.
* CopyKAT and SCEVAN were fastest; inferCNV and CaSpER were intermediate.

The consequence for a pan-cancer cohort is that no single caller is right for every sample,
and depth must be reported alongside every malignant-fraction estimate.
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

log = logging.getLogger("panspatial.malignant.cnv")

LOW_DEPTH_UMI = 3000     # below this, Numbat's advantage erodes
VERY_LOW_DEPTH_UMI = 1000


@dataclass(frozen=True)
class CallerChoice:
    primary: str
    cross_check: str
    rationale: str
    depth_caveat: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "primary": self.primary, "cross_check": self.cross_check,
            "rationale": self.rationale, "depth_caveat": self.depth_caveat,
        }


def select_cnv_caller(
    *, has_allele_info: bool, median_umi_per_cell: float, prioritise_speed: bool = False
) -> CallerChoice:
    """Pick a primary caller and a cross-check from allele availability and depth."""
    if median_umi_per_cell <= 0:
        raise ValueError("median_umi_per_cell must be positive")

    if has_allele_info and median_umi_per_cell >= LOW_DEPTH_UMI:
        return CallerChoice(
            "numbat", "copykat",
            "Allele information is available at adequate depth; Numbat performed best "
            "overall at separating tumour from normal in benchmark.",
            depth_caveat=None if median_umi_per_cell >= 2 * LOW_DEPTH_UMI else
            f"median {median_umi_per_cell:.0f} UMI/cell is close to where Numbat begins to degrade",
        )
    if has_allele_info:
        return CallerChoice(
            "copykat", "numbat",
            f"Allele information exists but median depth is {median_umi_per_cell:.0f} "
            "UMI/cell. Numbat's classification degrades sharply with depth while CopyKAT "
            "remains robust, so CopyKAT leads and Numbat cross-checks.",
            depth_caveat=f"low depth ({median_umi_per_cell:.0f} UMI/cell)",
        )
    if prioritise_speed:
        return CallerChoice(
            "scevan", "copykat",
            "Expression only, runtime-constrained: SCEVAN and CopyKAT were fastest in "
            "benchmark.",
            depth_caveat=None,
        )
    return CallerChoice(
        "copykat", "infercnv",
        "Expression only: CopyKAT was the best-performing expression-only tool and stayed "
        "accurate down to ~1,000 UMI/cell. inferCNV cross-checks and gives comparability "
        "with published atlases.",
        depth_caveat=(f"very low depth ({median_umi_per_cell:.0f} UMI/cell); report the "
                      "malignant fraction with this caveat attached")
        if median_umi_per_cell < VERY_LOW_DEPTH_UMI else None,
    )


def malignant_fraction_table(
    calls: Sequence[str], samples: Sequence[str], median_umi: dict[str, float], caller: str
) -> Any:
    """Per-sample malignant fraction, always reported next to depth.

    Cross-cohort comparison of tumour purity is only valid between samples of comparable
    depth, so the depth column travels with the estimate rather than living in a methods
    paragraph.
    """
    import pandas as pd

    calls = np.asarray(calls, dtype=object)
    samples = np.asarray(samples, dtype=object)
    if calls.shape != samples.shape:
        raise ValueError("calls and samples must describe the same cells")

    rows = []
    for s in sorted({str(x) for x in samples}):
        mask = samples.astype(str) == s
        n = int(mask.sum())
        n_mal = int(np.sum(calls[mask] == "malignant"))
        depth = float(median_umi.get(s, float("nan")))
        rows.append({
            "sample_id": s, "caller": caller, "n_cells": n,
            "n_malignant": n_mal,
            "malignant_fraction": round(n_mal / n, 4) if n else np.nan,
            "median_umi_per_cell": depth,
            "depth_comparable": bool(np.isfinite(depth) and depth >= LOW_DEPTH_UMI),
        })
    return pd.DataFrame(rows)


@register
class MalignantCompartment(AnalysisModule):
    name = "malignant"
    version = "0.1.0"
    phase = "4"
    description = "CNV-based malignant cell identification with depth-aware caller selection."
    evidence = Evidence.INFERRED
    replication_unit = ReplicationUnit.CELL
    supported_platforms = frozenset(
        {PlatformClass.SINGLE_CELL, PlatformClass.SPOT, PlatformClass.DISSOCIATED}
    )
    optional_dependencies = ("infercnvpy",)

    def run(
        self, ctx: ModuleContext, *, calls: Sequence[str] | None = None,
        samples: Sequence[str] | None = None, median_umi: dict[str, float] | None = None,
        **inputs: Any,
    ) -> ModuleResult:
        choice = select_cnv_caller(
            has_allele_info=ctx.param("has_allele_info", False),
            median_umi_per_cell=ctx.param("median_umi_per_cell", 5000.0),
            prioritise_speed=ctx.param("prioritise_speed", False),
        )
        if calls is None or samples is None:
            raise ModuleError(
                f"malignant.run requires calls= and samples=. Run {choice.primary} "
                f"(cross-check {choice.cross_check}) first — these are external R tools — "
                f"then pass the per-cell calls here. Rationale: {choice.rationale}"
            )
        table = malignant_fraction_table(
            calls, samples, median_umi or {}, choice.primary
        )
        notes = [choice.rationale]
        if choice.depth_caveat:
            notes.append("Depth caveat: " + choice.depth_caveat)
        notes.append(
            "Cross-cohort purity comparisons are valid only between samples flagged "
            "depth_comparable."
        )
        return self.result(ctx, table, notes=notes, caller=choice.to_dict())
