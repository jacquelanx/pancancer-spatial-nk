"""Spatial domain detection: platform-aware method choice and a disagreement-aware consensus.

A 2025 Nucleic Acids Research benchmark ranked GraphST first overall, closely followed by
BayesSpace, SpaGCN and STAGATE — but the ordering flips by platform. On 10x Visium, GraphST
led (mean ARI 0.552) ahead of STAGATE (0.515) and CCST (0.481); on legacy ST data BayesSpace
led decisively (0.642). No method wins everywhere, so the method is chosen per platform and
the choice travels with the result.

Where methods disagree about a spot, this module labels it boundary-uncertain rather than
picking a side. Domain boundaries are exactly where the biology is interesting and where
methods are least reliable, so silently resolving the disagreement discards the signal that
a boundary is uncertain.
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

log = logging.getLogger("panspatial.domains")

BOUNDARY_UNCERTAIN = "boundary_uncertain"

# Benchmark-informed defaults. (primary, second opinion, evidence)
PLATFORM_METHODS: dict[str, tuple[str, str, str]] = {
    "Visium":      ("graphst", "stagate", "GraphST led the Visium benchmark (mean ARI 0.552), STAGATE second (0.515)"),
    "Visium HD":   ("graphst", "banksy",  "GraphST default for 10x grids; BANKSY scales to HD bin counts"),
    "ST (legacy)": ("bayesspace", "graphst", "BayesSpace led decisively on legacy ST (ARI 0.642) ahead of Leiden (0.562)"),
    "Slide-seq":   ("bayesspace", "stagate", "Bayesian spatial prior suits low counts per bead"),
    "Stereo-seq":  ("graphst", "banksy",  "Grid-like geometry; BANKSY for very large bin counts"),
    "Xenium":      ("banksy", "cellcharter", "Segmented cells: neighbour-augmented features rather than spot smoothing"),
    "CosMx":       ("banksy", "cellcharter", "Segmented cells: neighbour-augmented features rather than spot smoothing"),
    "MERFISH":     ("banksy", "cellcharter", "Segmented cells: neighbour-augmented features rather than spot smoothing"),
}
FALLBACK = ("graphst", "stagate", "No platform-specific benchmark evidence; GraphST default with STAGATE as second opinion")


@dataclass(frozen=True)
class MethodChoice:
    platform: str
    primary: str
    second_opinion: str
    rationale: str

    def to_dict(self) -> dict[str, str]:
        return {
            "platform": self.platform, "primary": self.primary,
            "second_opinion": self.second_opinion, "rationale": self.rationale,
        }


def select_method(platform: str | None) -> MethodChoice:
    """Choose primary and second-opinion domain callers from benchmark evidence."""
    if platform is None:
        raise ModuleError("domain detection needs a platform; coordinates alone are not enough")
    primary, second, why = PLATFORM_METHODS.get(platform, FALLBACK)
    if platform not in PLATFORM_METHODS:
        log.warning("no benchmark entry for %r; using the default pair", platform)
    return MethodChoice(platform, primary, second, why)


def _align_labels(reference: np.ndarray, other: np.ndarray) -> np.ndarray:
    """Relabel ``other`` to maximise agreement with ``reference`` (Hungarian on contingency).

    Cluster labels are arbitrary integers; comparing them directly measures nothing.
    """
    from scipy.optimize import linear_sum_assignment

    ref_levels, ref_idx = np.unique(reference, return_inverse=True)
    oth_levels, oth_idx = np.unique(other, return_inverse=True)
    table = np.zeros((len(oth_levels), len(ref_levels)), dtype=float)
    np.add.at(table, (oth_idx, ref_idx), 1.0)
    rows, cols = linear_sum_assignment(-table)
    mapping = {int(r): ref_levels[c] for r, c in zip(rows, cols)}
    # Unmatched source clusters keep a distinct label rather than collapsing together.
    out = np.empty(len(other), dtype=object)
    for i, oi in enumerate(oth_idx):
        out[i] = mapping.get(int(oi), f"unmatched_{oth_levels[oi]}")
    return out


@dataclass
class DomainConsensus:
    labels: np.ndarray
    agreement: np.ndarray
    uncertain: np.ndarray
    methods: tuple[str, ...]
    mean_agreement: float

    @property
    def fraction_uncertain(self) -> float:
        return float(self.uncertain.mean()) if len(self.uncertain) else 0.0

    def table(self, index: Sequence[Any] | None = None):
        import pandas as pd
        return pd.DataFrame(
            {"domain": self.labels, "agreement": self.agreement,
             "boundary_uncertain": self.uncertain},
            index=index if index is not None else np.arange(len(self.labels)),
        )


def consensus_domains(
    labelings: Mapping[str, Sequence[Any]],
    *,
    min_agreement: float = 0.6,
) -> DomainConsensus:
    """Majority-vote domain labels across methods, flagging spots where methods disagree.

    Labelings are aligned to the first method before voting, since cluster indices are
    arbitrary. A spot whose winning label is supported by fewer than ``min_agreement`` of
    the methods is marked ``boundary_uncertain`` and carries that label instead.
    """
    if len(labelings) < 2:
        raise ValueError(
            "a consensus needs at least two labelings; with one method there is no "
            "disagreement to measure and no basis for calling a boundary uncertain"
        )
    names = list(labelings)
    arrays = [np.asarray(labelings[n], dtype=object) for n in names]
    n = len(arrays[0])
    if any(len(a) != n for a in arrays):
        raise ValueError("all labelings must cover the same spots")

    reference = arrays[0]
    aligned = [reference] + [_align_labels(reference, a) for a in arrays[1:]]

    labels = np.empty(n, dtype=object)
    agreement = np.zeros(n, dtype=float)
    for i in range(n):
        votes: dict[Any, int] = {}
        for a in aligned:
            votes[a[i]] = votes.get(a[i], 0) + 1
        winner, count = max(votes.items(), key=lambda kv: (kv[1], str(kv[0])))
        labels[i] = winner
        agreement[i] = count / len(aligned)

    uncertain = agreement < min_agreement
    labels_out = labels.copy()
    labels_out[uncertain] = BOUNDARY_UNCERTAIN
    mean_agree = float(agreement.mean())
    log.info(
        "domain consensus over %d methods: mean agreement %.2f, %d/%d spots "
        "boundary-uncertain", len(names), mean_agree, int(uncertain.sum()), n,
    )
    return DomainConsensus(labels_out, agreement, uncertain, tuple(names), mean_agree)


@register
class SpatialDomains(AnalysisModule):
    name = "domains"
    version = "0.1.0"
    phase = "4"
    description = "Platform-aware spatial domain detection with a disagreement-aware consensus."
    evidence = Evidence.INFERRED
    replication_unit = ReplicationUnit.SECTION
    supported_platforms = frozenset({PlatformClass.SINGLE_CELL, PlatformClass.SPOT})
    optional_dependencies = ("scanpy", "squidpy", "GraphST", "STAGATE", "banksy")

    def run(
        self, ctx: ModuleContext, *, labelings: Mapping[str, Sequence[Any]] | None = None,
        **inputs: Any,
    ) -> ModuleResult:
        self.check_platform(ctx)
        choice = select_method(ctx.platform)
        if not labelings:
            raise ModuleError(
                f"domains.run requires labelings={{method: labels}}. For {ctx.platform} run "
                f"{choice.primary} and {choice.second_opinion} first (external, GPU-bound), "
                "then pass both here. Rationale: " + choice.rationale
            )
        cons = consensus_domains(labelings, min_agreement=ctx.param("min_agreement", 0.6))
        return self.result(
            ctx, cons.table(index=inputs.get("index")),
            notes=[
                choice.rationale,
                f"{cons.fraction_uncertain:.1%} of spots are boundary-uncertain and are "
                "excluded from domain-specific statistics rather than assigned.",
            ],
            method_choice=choice.to_dict(),
            methods=list(cons.methods),
            mean_agreement=cons.mean_agreement,
            fraction_uncertain=cons.fraction_uncertain,
        )
