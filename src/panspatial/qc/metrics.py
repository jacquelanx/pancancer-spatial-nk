"""Quality control that adapts to the sample instead of imposing a global threshold.

A flat mitochondrial cut is the most common QC error in pan-cancer work: 10% deletes real
cells in kidney, liver and heart, where baseline mitochondrial content is high, and retains
debris in tissues where it is low. Thresholds here are derived per sample from the median
absolute deviation, and every threshold and every cell lost is recorded.

The second failure this module prevents is silent unit mixing. Visium coordinates arrive in
full-resolution image pixels; Xenium arrives in microns. A distance computed across that
boundary is meaningless, and nothing downstream can detect it. :func:`to_microns` refuses to
guess a scale factor.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
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
    require,
)

log = logging.getLogger("panspatial.qc")

# Platforms whose coordinates are natively in microns and need no conversion.
NATIVE_MICRON_PLATFORMS = {"Xenium", "CosMx", "MERFISH", "seqFISH", "STARmap"}


@dataclass
class Thresholds:
    """Per-metric cut points with the statistics that produced them."""

    metric: str
    lower: float
    upper: float
    median: float
    mad: float
    n_below: int = 0
    n_above: int = 0

    def to_dict(self) -> dict[str, Any]:
        return {
            "metric": self.metric, "lower": self.lower, "upper": self.upper,
            "median": self.median, "mad": self.mad,
            "n_below": self.n_below, "n_above": self.n_above,
        }


def mad_thresholds(
    values: Sequence[float],
    *,
    metric: str = "metric",
    nmads: float = 3.0,
    direction: str = "both",
    log_transform: bool = False,
) -> Thresholds:
    """Median ± ``nmads`` × MAD cut points for one metric within one sample.

    Args:
        values: the metric across cells in a single sample. Pooling samples defeats the point.
        direction: ``"both"``, ``"upper"`` (only cut high values, e.g. mitochondrial percent),
            or ``"lower"``.
        log_transform: apply ``log10(x + 1)`` before computing statistics, appropriate for
            count and feature distributions, which are right-skewed.

    A zero MAD — a constant metric, or one where over half the cells share a value — would
    make every cell an outlier. It is treated as no filtering, and logged.
    """
    if direction not in {"both", "upper", "lower"}:
        raise ValueError(f"direction must be both|upper|lower, got {direction!r}")
    x = np.asarray(values, dtype=float)
    if x.size == 0:
        raise ValueError(f"{metric}: no values supplied")
    finite = x[np.isfinite(x)]
    if finite.size == 0:
        raise ValueError(f"{metric}: no finite values")
    work = np.log10(finite + 1.0) if log_transform else finite

    med = float(np.median(work))
    mad = float(np.median(np.abs(work - med)) * 1.4826)  # scaled to a normal-consistent SD
    if mad == 0.0:
        log.warning(
            "%s: MAD is zero (over half the cells share a value); disabling this filter "
            "rather than flagging every cell as an outlier", metric,
        )
        return Thresholds(metric, -np.inf, np.inf, med, 0.0)

    lower = med - nmads * mad if direction in {"both", "lower"} else -np.inf
    upper = med + nmads * mad if direction in {"both", "upper"} else np.inf
    if log_transform:
        lower = 10**lower - 1.0 if np.isfinite(lower) else lower
        upper = 10**upper - 1.0 if np.isfinite(upper) else upper
        med = 10**med - 1.0
    return Thresholds(
        metric, float(lower), float(upper), float(med), mad,
        n_below=int(np.sum(finite < lower)), n_above=int(np.sum(finite > upper)),
    )


@dataclass
class QCReport:
    """What QC did to one sample, in a form that becomes a supplementary table."""

    sample_id: str
    n_before: int
    n_after: int
    thresholds: list[Thresholds] = field(default_factory=list)
    flags: dict[str, int] = field(default_factory=dict)

    @property
    def fraction_removed(self) -> float:
        return 1.0 - self.n_after / self.n_before if self.n_before else 0.0

    def to_row(self) -> dict[str, Any]:
        row: dict[str, Any] = {
            "sample_id": self.sample_id,
            "n_before": self.n_before,
            "n_after": self.n_after,
            "fraction_removed": round(self.fraction_removed, 4),
        }
        for t in self.thresholds:
            row[f"{t.metric}_lower"] = t.lower
            row[f"{t.metric}_upper"] = t.upper
            row[f"{t.metric}_median"] = t.median
        row.update({f"flag_{k}": v for k, v in self.flags.items()})
        return row

    def summary(self) -> str:
        return (
            f"{self.sample_id}: {self.n_before} -> {self.n_after} cells "
            f"({100*self.fraction_removed:.1f}% removed); "
            + ", ".join(f"{k}={v}" for k, v in sorted(self.flags.items()))
        )


def qc_mask(
    *,
    sample_id: str,
    percent_mito: Sequence[float],
    n_features: Sequence[float],
    n_counts: Sequence[float],
    nmads: float = 3.0,
    min_features: int = 200,
    max_mito_ceiling: float = 50.0,
) -> tuple[np.ndarray, QCReport]:
    """Build a keep-mask from per-sample adaptive thresholds.

    ``max_mito_ceiling`` is a backstop, not the primary filter: in a sample that is mostly
    dying cells, the MAD threshold adapts to the damage and would keep it. The adaptive cut
    does the work; the ceiling catches the pathological case.
    """
    mito = np.asarray(percent_mito, dtype=float)
    feats = np.asarray(n_features, dtype=float)
    counts = np.asarray(n_counts, dtype=float)
    if not (mito.shape == feats.shape == counts.shape):
        raise ValueError("percent_mito, n_features and n_counts must describe the same cells")

    t_mito = mad_thresholds(mito, metric="percent_mito", nmads=nmads, direction="upper")
    t_feat = mad_thresholds(feats, metric="n_features", nmads=nmads, log_transform=True)
    t_cnt = mad_thresholds(counts, metric="n_counts", nmads=nmads, log_transform=True)

    fail_mito = (mito > t_mito.upper) | (mito > max_mito_ceiling)
    fail_feat = (feats < t_feat.lower) | (feats > t_feat.upper) | (feats < min_features)
    fail_cnt = (counts < t_cnt.lower) | (counts > t_cnt.upper)
    keep = ~(fail_mito | fail_feat | fail_cnt)

    report = QCReport(
        sample_id=sample_id,
        n_before=int(mito.size),
        n_after=int(keep.sum()),
        thresholds=[t_mito, t_feat, t_cnt],
        flags={
            "high_mito": int(fail_mito.sum()),
            "feature_outlier": int(fail_feat.sum()),
            "count_outlier": int(fail_cnt.sum()),
        },
    )
    log.info(report.summary())
    if report.fraction_removed > 0.5:
        log.warning(
            "%s: QC removed %.0f%% of cells. Inspect the sample before using it; a filter "
            "this aggressive usually indicates a failed dissociation rather than a strict "
            "threshold.", sample_id, 100 * report.fraction_removed,
        )
    return keep, report


def to_microns(
    coords: np.ndarray,
    *,
    platform: str,
    microns_per_unit: float | None = None,
) -> np.ndarray:
    """Convert spatial coordinates to microns, refusing to guess a scale factor.

    Imaging platforms report microns natively. Visium reports full-resolution image pixels,
    and the conversion depends on the scalefactors JSON for that specific capture area —
    there is no universal constant. Passing the wrong one silently rescales every distance
    in the study, so an unknown scale is an error rather than a default.
    """
    arr = np.asarray(coords, dtype=float)
    if arr.ndim != 2 or arr.shape[1] < 2:
        raise ValueError(f"coords must be (n, >=2), got {arr.shape}")
    arr = arr[:, :2]

    if platform in NATIVE_MICRON_PLATFORMS:
        if microns_per_unit not in (None, 1.0):
            log.warning(
                "%s reports microns natively; applying the supplied factor %.4f anyway",
                platform, microns_per_unit,
            )
            return arr * microns_per_unit
        return arr
    if microns_per_unit is None:
        raise ModuleError(
            f"{platform} coordinates are not in microns and no scale factor was given. "
            "Read microns_per_unit from the platform's scalefactors (for Visium: spot "
            "diameter in full-resolution pixels divided by 55 µm). Guessing here would "
            "silently rescale every distance in the study."
        )
    if microns_per_unit <= 0:
        raise ValueError(f"microns_per_unit must be positive, got {microns_per_unit}")
    return arr * microns_per_unit


def spatial_qc(
    coords: np.ndarray,
    *,
    in_tissue: Sequence[bool] | None = None,
    min_neighbour_distance_um: float = 0.0,
    max_isolation_um: float | None = None,
) -> tuple[np.ndarray, dict[str, int]]:
    """Flag off-tissue spots, coincident points, and isolated debris.

    Coincident points (distance zero) usually indicate a segmentation failure on imaging
    platforms; isolated points far from any neighbour are typically tissue fragments or
    debris that distort neighbourhood graphs.
    """
    from scipy.spatial import cKDTree

    arr = np.asarray(coords, dtype=float)
    keep = np.ones(len(arr), dtype=bool)
    flags = {"off_tissue": 0, "coincident": 0, "isolated": 0}

    if in_tissue is not None:
        on = np.asarray(in_tissue, dtype=bool)
        flags["off_tissue"] = int((~on).sum())
        keep &= on

    if len(arr) > 1:
        tree = cKDTree(arr)
        nn = tree.query(arr, k=2)[0][:, 1]
        if min_neighbour_distance_um > 0:
            bad = nn < min_neighbour_distance_um
            flags["coincident"] = int(bad.sum())
            keep &= ~bad
        if max_isolation_um is not None:
            lonely = nn > max_isolation_um
            flags["isolated"] = int(lonely.sum())
            keep &= ~lonely

    log.info("spatial QC: %d/%d retained; %s", int(keep.sum()), len(arr),
             ", ".join(f"{k}={v}" for k, v in flags.items()))
    return keep, flags


@register
class QualityControl(AnalysisModule):
    name = "qc"
    version = "0.1.0"
    phase = "2"
    description = "Per-sample adaptive QC, doublet and ambient-RNA removal, coordinate normalisation."
    evidence = Evidence.MEASURED
    replication_unit = ReplicationUnit.SECTION
    supported_platforms = frozenset(
        {PlatformClass.SINGLE_CELL, PlatformClass.SPOT, PlatformClass.REGION,
         PlatformClass.DISSOCIATED}
    )
    optional_dependencies = ("scanpy", "scikit-learn")

    def run(self, ctx: ModuleContext, *, adata: Any = None, **inputs: Any) -> ModuleResult:
        import pandas as pd

        sc = require("scanpy", extra="spatial", purpose="quality control on an AnnData")
        if adata is None:
            raise ModuleError("qc.run requires adata=")

        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
        adata.var["ribo"] = adata.var_names.str.upper().str.match(r"^RP[SL]")
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=["mt", "ribo"], inplace=True, percent_top=None, log1p=False
        )

        reports: list[QCReport] = []
        sample_key = ctx.param("sample_key", "sample_id")
        samples = (
            adata.obs[sample_key].astype(str).to_numpy()
            if sample_key in adata.obs
            else np.full(adata.n_obs, ctx.sample_id)
        )
        keep_all = np.zeros(adata.n_obs, dtype=bool)
        for s in np.unique(samples):
            idx = np.flatnonzero(samples == s)
            keep, rep = qc_mask(
                sample_id=str(s),
                percent_mito=adata.obs["pct_counts_mt"].to_numpy()[idx],
                n_features=adata.obs["n_genes_by_counts"].to_numpy()[idx],
                n_counts=adata.obs["total_counts"].to_numpy()[idx],
                nmads=ctx.param("nmads", 3.0),
                min_features=ctx.param("min_features", 200),
            )
            keep_all[idx] = keep
            reports.append(rep)

        adata.obs["qc_pass"] = keep_all
        table = pd.DataFrame([r.to_row() for r in reports])
        return self.result(
            ctx, table,
            tool_versions={"scanpy": sc.__version__},
            n_cells_in=int(adata.n_obs),
            n_cells_pass=int(keep_all.sum()),
            notes=["Thresholds are per sample; the table is the audit trail for the supplement."],
        )
