"""Co-expression modules and whether they survive in another cohort.

hdWGCNA builds context-specific networks per cell population, using metacell aggregation to
counter single-cell sparsity. The pan-cancer question is not what modules exist in one
cancer type but which ones are *preserved* in others — and a module can have a perfectly
good eigengene in a second cohort while its internal structure has fallen apart.

:func:`module_preservation` implements a two-statistic reduction of Langfelder's Zsummary:
density preservation (do module genes stay correlated?) and connectivity preservation (do
the same genes stay hub genes?), each z-scored against permuted module memberships. The
conventional reading of Zsummary applies: below 2 is no preservation, 2–10 weak to
moderate, above 10 strong.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
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

log = logging.getLogger("panspatial.comodules")


def _adjacency(expr: np.ndarray, power: int = 6) -> np.ndarray:
    """Signed-hybrid WGCNA adjacency from a gene x sample expression matrix."""
    if expr.shape[0] < 2:
        raise ValueError("need at least two genes")
    with np.errstate(invalid="ignore"):
        corr = np.corrcoef(expr)
    corr = np.nan_to_num(corr, nan=0.0)
    return np.clip(corr, 0.0, None) ** power


def _preservation_statistics(expr: np.ndarray, power: int = 6) -> tuple[float, float]:
    """(mean adjacency, mean intramodular connectivity) for a set of genes."""
    adj = _adjacency(expr, power)
    n = adj.shape[0]
    off = ~np.eye(n, dtype=bool)
    mean_adj = float(adj[off].mean()) if n > 1 else 0.0
    k = adj.sum(axis=1) - np.diag(adj)
    return mean_adj, k


@dataclass
class PreservationResult:
    module: str
    n_genes: int
    z_density: float
    z_connectivity: float
    z_summary: float
    interpretation: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "module": self.module, "n_genes": self.n_genes,
            "z_density": round(self.z_density, 3),
            "z_connectivity": round(self.z_connectivity, 3),
            "z_summary": round(self.z_summary, 3),
            "preserved": self.z_summary >= 2.0,
            "strongly_preserved": self.z_summary >= 10.0,
            "interpretation": self.interpretation,
        }


def module_preservation(
    reference_expr: Any,
    test_expr: Any,
    modules: Mapping[str, Sequence[str]],
    *,
    gene_names: Sequence[str] | None = None,
    n_permutations: int = 100,
    power: int = 6,
    rng: np.random.Generator | None = None,
) -> Any:
    """Zsummary-style preservation of reference modules in a test cohort.

    Args:
        reference_expr / test_expr: gene × sample matrices with the same gene ordering,
            or DataFrames indexed by gene.
        modules: ``{module_name: [genes]}`` from the reference cohort.
        n_permutations: random gene sets of matching size form the null. 100 is enough to
            separate "preserved" from "not"; raise it for a published Z value.
    """
    import pandas as pd

    if hasattr(reference_expr, "index"):
        gene_names = list(reference_expr.index)
        reference_expr = reference_expr.to_numpy(dtype=float)
    if hasattr(test_expr, "index"):
        test_expr = test_expr.to_numpy(dtype=float)
    reference_expr = np.asarray(reference_expr, dtype=float)
    test_expr = np.asarray(test_expr, dtype=float)
    if gene_names is None:
        raise ValueError("gene_names is required when matrices are passed as arrays")
    if reference_expr.shape[0] != test_expr.shape[0] != len(gene_names):
        raise ValueError("reference and test must have the same genes in the same order")
    rng = rng or np.random.default_rng(0)
    index = {g: i for i, g in enumerate(gene_names)}
    n_genes = len(gene_names)

    rows = []
    for name, genes in modules.items():
        idx = np.array([index[g] for g in genes if g in index], dtype=int)
        if len(idx) < 5:
            log.warning("module %s has %d genes present; too small to assess", name, len(idx))
            continue
        obs_density, obs_k = _preservation_statistics(test_expr[idx], power)
        ref_density, ref_k = _preservation_statistics(reference_expr[idx], power)
        obs_conn = float(np.corrcoef(ref_k, obs_k)[0, 1]) if len(idx) > 2 else 0.0
        obs_conn = 0.0 if not np.isfinite(obs_conn) else obs_conn

        null_density, null_conn = [], []
        for _ in range(n_permutations):
            perm = rng.choice(n_genes, size=len(idx), replace=False)
            d, k_test = _preservation_statistics(test_expr[perm], power)
            _, k_ref = _preservation_statistics(reference_expr[perm], power)
            null_density.append(d)
            c = np.corrcoef(k_ref, k_test)[0, 1] if len(perm) > 2 else 0.0
            null_conn.append(0.0 if not np.isfinite(c) else float(c))

        def z(obs: float, null: list[float]) -> float:
            sd = float(np.std(null, ddof=1))
            return float((obs - float(np.mean(null))) / sd) if sd > 0 else 0.0

        zd, zc = z(obs_density, null_density), z(obs_conn, null_conn)
        zs = float(np.mean([zd, zc]))
        interp = (
            "strong preservation" if zs >= 10 else
            "weak to moderate preservation" if zs >= 2 else
            "no evidence of preservation; do not report this module as pan-cancer"
        )
        rows.append(PreservationResult(str(name), len(idx), zd, zc, zs, interp))

    frame = pd.DataFrame([r.to_dict() for r in rows])
    if not frame.empty:
        frame = frame.sort_values("z_summary", ascending=False).reset_index(drop=True)
        log.info("%d/%d modules preserved (Zsummary >= 2)",
                 int(frame["preserved"].sum()), len(frame))
    return frame


@register
class CoexpressionModules(AnalysisModule):
    name = "comodules"
    version = "0.1.0"
    phase = "5"
    description = "hdWGCNA/cNMF co-expression modules with cross-cohort preservation testing."
    evidence = Evidence.INFERRED
    replication_unit = ReplicationUnit.COHORT
    supported_platforms = frozenset(
        {PlatformClass.SINGLE_CELL, PlatformClass.SPOT, PlatformClass.DISSOCIATED}
    )
    optional_dependencies = ("hdWGCNA (R)", "cnmf", "mofapy2")

    def run(
        self, ctx: ModuleContext, *, reference_expr: Any = None, test_expr: Any = None,
        modules: Mapping[str, Sequence[str]] | None = None, **inputs: Any,
    ) -> ModuleResult:
        if reference_expr is None or test_expr is None or not modules:
            raise ModuleError(
                "comodules.run requires reference_expr=, test_expr= and modules=. Build "
                "modules with hdWGCNA (R) or cNMF first; this module tests whether they "
                "survive in another cohort, which is what makes them pan-cancer."
            )
        frame = module_preservation(
            reference_expr, test_expr, modules,
            n_permutations=ctx.param("n_permutations", 100),
            power=ctx.param("soft_power", 6),
            rng=np.random.default_rng(ctx.seed),
        )
        preserved = frozenset(frame.loc[frame["preserved"], "module"]) if len(frame) else frozenset()
        return self.result(
            ctx, frame, reliable_entities=preserved,
            notes=[
                "Zsummary below 2 means no evidence of preservation; such modules are "
                "cohort-specific and must not be described as pan-cancer.",
                "Preservation is assessed on module structure, not on eigengene "
                "correlation: a module can keep its eigengene while losing its topology.",
            ],
        )
