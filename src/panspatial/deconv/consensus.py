"""Consensus deconvolution across methods, because the benchmarks disagree.

A 2023 Nature Communications benchmark of 18 methods over 50 datasets put CARD,
cell2location and Tangram first. A 2023 Bioinformatics benchmark put RCTD and cell2location
first on AUPR and Jensen–Shannon divergence. Only cell2location appears at the top of both.

Rankings shift with tissue, platform, reference quality and metric, so picking one winner
from one benchmark imports that benchmark's tissue bias into every result. This module runs
several and reports the *disagreement* as a first-class output: where methods diverge, the
estimate is uncertain, and that is information rather than noise to be averaged away.
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

log = logging.getLogger("panspatial.deconv")

# The union of the top performers across published benchmarks.
DEFAULT_METHODS = ("cell2location", "rctd", "card", "tangram")


@dataclass
class ConsensusResult:
    weights: Any                  # DataFrame: spots x cell types, median across methods
    dispersion: Any               # DataFrame: spots x cell types, inter-method IQR
    concordance: Any              # DataFrame: per cell type, mean pairwise agreement
    methods: tuple[str, ...]

    @property
    def concordant_types(self) -> frozenset[str]:
        frame = self.concordance
        return frozenset(frame.loc[frame["concordant"], "cell_type"])


def _align(method_weights: Mapping[str, Any]) -> tuple[list[str], list[str], np.ndarray]:
    """Stack per-method weight frames onto a common spot × cell-type grid."""
    import pandas as pd

    if len(method_weights) < 2:
        raise ValueError(
            "consensus needs at least two methods; running one method and calling it a "
            "consensus is the failure mode this function exists to prevent"
        )
    frames = {m: pd.DataFrame(w) for m, w in method_weights.items()}
    spots = sorted(set.intersection(*(set(f.index) for f in frames.values())))
    types = sorted(set.intersection(*(set(f.columns) for f in frames.values())))
    if not spots or not types:
        raise ValueError("methods share no common spots or cell types")

    dropped_types = set.union(*(set(f.columns) for f in frames.values())) - set(types)
    if dropped_types:
        log.warning(
            "cell types missing from at least one method and excluded from the consensus: %s",
            ", ".join(sorted(map(str, dropped_types))),
        )
    stack = np.stack([frames[m].loc[spots, types].to_numpy(dtype=float)
                      for m in frames], axis=0)
    return list(spots), list(types), stack


def consensus_weights(
    method_weights: Mapping[str, Any],
    *,
    min_concordance: float = 0.5,
    normalise: bool = True,
) -> ConsensusResult:
    """Median weight per spot and cell type, with inter-method agreement.

    Concordance is the mean pairwise Spearman correlation across methods for a cell type,
    computed over spots. A type on which methods disagree is flagged: its consensus weight
    may still be reported, but not as a confident estimate.
    """
    import pandas as pd
    from scipy.stats import spearmanr

    spots, types, stack = _align(method_weights)
    methods = tuple(method_weights)

    median = np.median(stack, axis=0)
    q75, q25 = np.percentile(stack, [75, 25], axis=0)
    if normalise:
        totals = median.sum(axis=1, keepdims=True)
        median = np.divide(median, totals, out=np.zeros_like(median), where=totals > 0)

    rows = []
    for j, ct in enumerate(types):
        pair_rhos = []
        for a in range(len(methods)):
            for b in range(a + 1, len(methods)):
                x, y = stack[a, :, j], stack[b, :, j]
                if np.std(x) == 0 or np.std(y) == 0:
                    continue
                rho = spearmanr(x, y)[0]
                if np.isfinite(rho):
                    pair_rhos.append(float(rho))
        mean_rho = float(np.mean(pair_rhos)) if pair_rhos else float("nan")
        rows.append({
            "cell_type": ct,
            "mean_pairwise_spearman": round(mean_rho, 4) if np.isfinite(mean_rho) else np.nan,
            "n_method_pairs": len(pair_rhos),
            "median_weight": round(float(np.mean(median[:, j])), 5),
            "mean_iqr": round(float(np.mean(q75[:, j] - q25[:, j])), 5),
            "concordant": bool(np.isfinite(mean_rho) and mean_rho >= min_concordance),
        })
    concordance = pd.DataFrame(rows).sort_values("mean_pairwise_spearman", ascending=False)
    discordant = concordance.loc[~concordance["concordant"], "cell_type"].tolist()
    if discordant:
        log.warning(
            "methods disagree (mean pairwise Spearman < %.2f) for: %s. Their consensus "
            "weights are reported with the disagreement, not as confident estimates.",
            min_concordance, ", ".join(map(str, discordant)),
        )
    return ConsensusResult(
        weights=pd.DataFrame(median, index=spots, columns=types),
        dispersion=pd.DataFrame(q75 - q25, index=spots, columns=types),
        concordance=concordance.reset_index(drop=True),
        methods=methods,
    )


def mode_agreement(full_weights: Any, doublet_calls: Mapping[str, Sequence[str]]) -> Any:
    """Compare RCTD full-mode weights against doublet-mode calls, per cell type.

    Full mode returns a weight for every type; doublet mode commits to at most two. A type
    with substantial full-mode weight but a near-zero doublet call rate is being carried by
    reference similarity rather than by evidence of presence.
    """
    import pandas as pd

    frame = pd.DataFrame(full_weights)
    n = len(frame)
    called = {ct: 0 for ct in frame.columns}
    for calls in doublet_calls.values():
        for ct in set(calls):
            if ct in called:
                called[ct] += 1

    rows = []
    for ct in frame.columns:
        rate = called[ct] / n if n else 0.0
        mean_w = float(frame[ct].mean())
        rows.append({
            "cell_type": ct,
            "full_mean_weight": round(mean_w, 5),
            "doublet_call_rate": round(rate, 5),
            "similarity_carried": bool(mean_w > 0.05 and rate < 0.02),
        })
    out = pd.DataFrame(rows)
    flagged = out.loc[out["similarity_carried"], "cell_type"].tolist()
    if flagged:
        log.warning(
            "full-mode weight without doublet-mode support (reference similarity, not "
            "presence): %s", ", ".join(map(str, flagged)),
        )
    return out


@register
class Deconvolution(AnalysisModule):
    name = "deconv"
    version = "0.1.0"
    phase = "4"
    description = "Consensus spot deconvolution across cell2location, RCTD, CARD and Tangram."
    evidence = Evidence.INFERRED
    replication_unit = ReplicationUnit.SPOT
    supported_platforms = frozenset({PlatformClass.SPOT})
    optional_dependencies = ("cell2location", "scanpy", "tangram-sc")

    def run(
        self,
        ctx: ModuleContext,
        *,
        method_weights: Mapping[str, Any] | None = None,
        reliability: Any = None,
        **inputs: Any,
    ) -> ModuleResult:
        if not method_weights:
            raise ModuleError(
                "deconv.run requires method_weights={method: DataFrame}. Run the individual "
                "deconvolution tools first (they are heavy, external, and best driven by the "
                "workflow), then pass their outputs here for consensus."
            )
        self.check_platform(ctx)
        result = consensus_weights(
            method_weights, min_concordance=ctx.param("min_concordance", 0.5)
        )
        reliable = result.concordant_types
        notes = [
            f"Consensus of {len(result.methods)} methods: {', '.join(result.methods)}.",
            "Weights are proportions of a spot's signal, never cell counts.",
        ]
        if reliability is not None:
            gate = set(getattr(reliability, "reliable", reliability))
            reliable = frozenset(reliable & gate)
            notes.append(
                "Intersected with the simulated-mixture reliability gate; a type must both "
                "be recovered accurately and be agreed on by methods."
            )
        return self.result(
            ctx, result.concordance, reliable_entities=reliable, notes=notes,
            methods=list(result.methods),
            n_spots=int(len(result.weights)), n_cell_types=int(result.weights.shape[1]),
        )
