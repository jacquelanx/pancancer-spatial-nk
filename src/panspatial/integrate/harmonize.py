"""Integration that refuses to remove the variable the study is about.

In a pan-cancer cohort assembled from public depositions, batch and cancer type are very
nearly the same variable: each tumour type arrives from different labs, chemistries and
years. Correcting "batch" globally therefore removes the biology.

:func:`confounding` measures the association with Cramér's V and
:func:`plan_integration` refuses a global plan when it is high, returning a within-cancer-type
plan instead. That refusal is the module's main contribution — the integration calls
themselves are thin.
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
    require,
)

log = logging.getLogger("panspatial.integrate")

METHODS = ("harmony", "scvi", "scanvi", "none")


def cramers_v(a: Sequence[Any], b: Sequence[Any], *, bias_correct: bool = True) -> float:
    """Cramér's V between two categorical vectors, in [0, 1].

    1.0 means one variable is perfectly predictable from the other — for batch and cancer
    type, that no correction can separate them. Bias correction (Bergsma) is applied by
    default because the uncorrected statistic is inflated for small samples with many levels,
    which is exactly the regime of a public pan-cancer cohort.
    """
    a = np.asarray(a)
    b = np.asarray(b)
    if a.shape != b.shape:
        raise ValueError(f"vectors must be the same length, got {a.shape} and {b.shape}")
    if a.size == 0:
        raise ValueError("cannot compute association on empty vectors")

    a_levels, a_idx = np.unique(a, return_inverse=True)
    b_levels, b_idx = np.unique(b, return_inverse=True)
    r, k = len(a_levels), len(b_levels)
    if r < 2 or k < 2:
        return 0.0

    table = np.zeros((r, k), dtype=float)
    np.add.at(table, (a_idx, b_idx), 1.0)
    n = table.sum()
    expected = np.outer(table.sum(1), table.sum(0)) / n
    with np.errstate(divide="ignore", invalid="ignore"):
        chi2 = np.nansum(np.where(expected > 0, (table - expected) ** 2 / expected, 0.0))
    phi2 = chi2 / n

    if bias_correct:
        phi2 = max(0.0, phi2 - (k - 1) * (r - 1) / max(n - 1, 1))
        r = r - (r - 1) ** 2 / max(n - 1, 1)
        k = k - (k - 1) ** 2 / max(n - 1, 1)
    denom = min(k - 1, r - 1)
    return float(np.sqrt(phi2 / denom)) if denom > 0 else 0.0


@dataclass
class IntegrationPlan:
    """What may be corrected, over what, and what the resulting claims may say."""

    scope: str                    # "within_cancer_type" or "global"
    batch_key: str
    groups: dict[str, int]        # group label -> n patients
    confounding: float
    method: str
    rationale: str
    cross_group_comparison: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "scope": self.scope, "batch_key": self.batch_key,
            "n_groups": len(self.groups), "groups": self.groups,
            "confounding_cramers_v": round(self.confounding, 4),
            "method": self.method, "rationale": self.rationale,
            "cross_group_comparison": self.cross_group_comparison,
        }


def plan_integration(
    batch: Sequence[Any],
    cancer_type: Sequence[Any],
    patient: Sequence[Any],
    *,
    method: str = "harmony",
    confounding_threshold: float = 0.5,
    min_patients_per_group: int = 2,
) -> IntegrationPlan:
    """Decide the integration scope from the actual confounding in this cohort.

    Above ``confounding_threshold`` a global correction over ``batch`` would also remove
    cancer-type signal, so the plan is forced to within-cancer-type. Below it, the cohort
    genuinely has batches spanning cancer types and global correction is defensible.
    """
    if method not in METHODS:
        raise ValueError(f"method must be one of {METHODS}, got {method!r}")
    batch = np.asarray(batch)
    cancer_type = np.asarray(cancer_type)
    patient = np.asarray(patient)
    if not (len(batch) == len(cancer_type) == len(patient)):
        raise ValueError("batch, cancer_type and patient must describe the same cells")

    v = cramers_v(batch, cancer_type)
    groups = {
        str(ct): int(len(np.unique(patient[cancer_type == ct])))
        for ct in np.unique(cancer_type)
    }
    thin = {g: n for g, n in groups.items() if n < min_patients_per_group}
    if thin:
        log.warning(
            "cancer type(s) with fewer than %d patients will not be integrated, only "
            "normalised: %s", min_patients_per_group, ", ".join(sorted(thin)),
        )

    if v >= confounding_threshold:
        return IntegrationPlan(
            scope="within_cancer_type", batch_key="patient", groups=groups, confounding=v,
            method=method,
            rationale=(
                f"Batch and cancer type are associated at Cramer's V = {v:.2f} "
                f"(threshold {confounding_threshold}). A global correction over batch would "
                "also remove cancer-type biology, so correction is restricted to within "
                "cancer type, over patient."
            ),
            cross_group_comparison=(
                "Compare across cancer types on shared cell-type labels and per-sample "
                "summary statistics — replication across cohorts, not a shared latent space."
            ),
        )
    return IntegrationPlan(
        scope="global", batch_key="batch", groups=groups, confounding=v, method=method,
        rationale=(
            f"Batch and cancer type are only weakly associated (Cramer's V = {v:.2f}); "
            "batches span cancer types, so a global correction does not remove the "
            "biological contrast."
        ),
        cross_group_comparison="Cross-cancer-type comparison in the shared embedding is defensible.",
    )


def scib_summary(adata: Any, *, label_key: str, batch_key: str, embed: str) -> dict[str, float]:
    """Batch-mixing and biological-conservation metrics for one integration.

    Both halves are reported. An integration that maximises mixing while collapsing known
    cell-type structure has failed, however good the embedding looks.
    """
    scib = require("scib_metrics", extra="spatial", purpose="integration benchmarking")
    from scib_metrics.benchmark import Benchmarker

    bm = Benchmarker(adata, batch_key=batch_key, label_key=label_key, embedding_obsm_keys=[embed])
    bm.benchmark()
    frame = bm.get_results(min_max_scale=False)
    row = frame.loc[embed].to_dict()
    out = {str(k): float(v) for k, v in row.items() if isinstance(v, (int, float))}
    if out.get("Bio conservation", 1.0) < 0.5:
        log.warning(
            "biological conservation is %.2f: this integration is mixing batches by "
            "destroying cell-type structure", out.get("Bio conservation", float("nan")),
        )
    return out


@register
class Integration(AnalysisModule):
    name = "integrate"
    version = "0.1.0"
    phase = "2"
    description = "Confounding-aware batch integration, within cancer type over patient."
    evidence = Evidence.INFERRED
    replication_unit = ReplicationUnit.COHORT
    supported_platforms = frozenset(
        {PlatformClass.SINGLE_CELL, PlatformClass.SPOT, PlatformClass.DISSOCIATED}
    )
    optional_dependencies = ("scanpy", "harmonypy", "scvi-tools", "scib-metrics")

    def run(self, ctx: ModuleContext, *, adata: Any = None, **inputs: Any) -> ModuleResult:
        import pandas as pd

        sc = require("scanpy", extra="spatial", purpose="integration")
        if adata is None:
            raise ModuleError("integrate.run requires adata=")

        batch_key = ctx.param("batch_key", "batch")
        cancer_key = ctx.param("cancer_key", "cancer_type")
        patient_key = ctx.param("patient_key", "patient")
        for key in (batch_key, cancer_key, patient_key):
            if key not in adata.obs:
                raise ModuleError(f"adata.obs is missing {key!r}, required to plan integration")

        plan = plan_integration(
            adata.obs[batch_key], adata.obs[cancer_key], adata.obs[patient_key],
            method=ctx.param("method", "harmony"),
            confounding_threshold=ctx.param("confounding_threshold", 0.5),
        )
        log.info("integration plan: %s — %s", plan.scope, plan.rationale)

        sc.pp.pca(adata, n_comps=ctx.param("n_pcs", 30), random_state=ctx.seed)
        versions = {"scanpy": sc.__version__}
        if plan.method == "harmony":
            hp = require("harmonypy", extra="spatial", purpose="Harmony integration")
            sc.external.pp.harmony_integrate(
                adata, key=plan.batch_key, basis="X_pca", adjusted_basis="X_integrated",
                random_state=ctx.seed,
            )
            versions["harmonypy"] = getattr(hp, "__version__", "unknown")
        elif plan.method in {"scvi", "scanvi"}:
            scvi = require("scvi", extra="spatial", purpose="scVI integration")
            scvi.settings.seed = ctx.seed
            scvi.model.SCVI.setup_anndata(adata, batch_key=plan.batch_key)
            model = scvi.model.SCVI(adata)
            model.train()
            adata.obsm["X_integrated"] = model.get_latent_representation()
            versions["scvi-tools"] = scvi.__version__
        else:
            adata.obsm["X_integrated"] = adata.obsm["X_pca"]

        table = pd.DataFrame([{"group": g, "n_patients": n, **plan.to_dict()}
                              for g, n in plan.groups.items()])
        return self.result(
            ctx, table, tool_versions=versions,
            notes=[plan.rationale, plan.cross_group_comparison],
            scope=plan.scope, confounding_cramers_v=plan.confounding,
        )
