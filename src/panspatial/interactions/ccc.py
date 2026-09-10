"""Cell–cell communication: spatial constraint first, then consensus across scorers.

The LIANA comparison of seven methods and sixteen resources found low overlap among
top-ranked interactions — the tools disagree because their scoring strategies differ
(specificity-focused for CellPhoneDB and CellChat, intensity-focused for others), not
because one is correct. Two things reduce that noise:

1. **A physical constraint.** Spatial data lets us require that sender and receiver are
   actually within signalling range. An interaction whose participants never neighbour each
   other is not a candidate however high it scores.
2. **Rank aggregation.** Robust rank aggregation across scorers surfaces the interactions
   that every method agrees on, which is a different and more defensible claim than any
   single tool's top hit.

Nothing here converts co-expression into communication. That distinction is preserved in
the output vocabulary.
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

log = logging.getLogger("panspatial.interactions")

METHODS = {
    "cellphonedb": "Permutation-based specificity; robust to noisy input and resource error",
    "cellchat": "Pathway-level aggregation with tissue-architecture visualisation",
    "nichenet": "Ligand activity from downstream target expression in the receiver",
    "commot": "Collective optimal transport over coordinates, modelling ligand competition",
    "spacci": "Spatially informed scoring; highest mean F1 in a recent spatial CCC benchmark",
    "liana": "Consensus framework across methods and resources",
}

# Typical upper bound for juxtacrine/paracrine range. Diffusible factors travel further;
# state the assumption rather than hiding it in a default.
DEFAULT_MAX_DISTANCE_UM = 30.0


def contact_frequency(
    coords: np.ndarray,
    cell_types: Sequence[str],
    *,
    max_distance_um: float = DEFAULT_MAX_DISTANCE_UM,
) -> Any:
    """Observed neighbour frequency for every ordered cell-type pair in one section.

    This is the physical prior an LR score is filtered against: a sender-receiver pair that
    never comes within range cannot be signalling by contact, whatever its expression.
    """
    import pandas as pd
    from scipy.spatial import cKDTree

    coords = np.asarray(coords, dtype=float)
    cell_types = np.asarray(cell_types, dtype=object)
    if len(coords) != len(cell_types):
        raise ValueError("coords and cell_types must describe the same cells")
    if max_distance_um <= 0:
        raise ValueError("max_distance_um must be positive")

    tree = cKDTree(coords)
    pairs = tree.query_pairs(max_distance_um, output_type="ndarray")
    types = sorted({str(t) for t in cell_types})
    index = {t: i for i, t in enumerate(types)}
    counts = np.zeros((len(types), len(types)), dtype=float)
    for i, j in pairs:
        a, b = index[str(cell_types[i])], index[str(cell_types[j])]
        counts[a, b] += 1
        counts[b, a] += 1

    n_by_type = np.array([np.sum(cell_types.astype(str) == t) for t in types], dtype=float)
    rows = []
    for a, ta in enumerate(types):
        for b, tb in enumerate(types):
            rows.append({
                "sender": ta, "receiver": tb,
                "n_contacts": int(counts[a, b]),
                "contacts_per_sender_cell": round(counts[a, b] / n_by_type[a], 4)
                if n_by_type[a] else 0.0,
                "in_contact": bool(counts[a, b] > 0),
            })
    return pd.DataFrame(rows)


def spatial_constraint(
    interactions: Any,
    contacts: Any,
    *,
    min_contacts: int = 1,
    sender_col: str = "sender",
    receiver_col: str = "receiver",
) -> Any:
    """Drop interactions whose sender and receiver are never within signalling range.

    Adds a ``spatially_supported`` column rather than silently filtering, so the number of
    expression-only calls that failed the physical check is reportable.
    """
    import pandas as pd

    interactions = pd.DataFrame(interactions).copy()
    contacts = pd.DataFrame(contacts)
    for col in (sender_col, receiver_col):
        if col not in interactions.columns:
            raise ValueError(f"interactions is missing the {col!r} column")

    lookup = {
        (str(r[sender_col]), str(r[receiver_col])): int(r["n_contacts"])
        for _, r in contacts.iterrows()
    }
    n = interactions[[sender_col, receiver_col]].apply(
        lambda r: lookup.get((str(r[sender_col]), str(r[receiver_col])), 0), axis=1
    )
    interactions["n_contacts"] = n
    interactions["spatially_supported"] = n >= min_contacts
    dropped = int((~interactions["spatially_supported"]).sum())
    log.info(
        "spatial constraint: %d/%d interactions have sender and receiver in contact "
        "(%d expression-only calls flagged)",
        len(interactions) - dropped, len(interactions), dropped,
    )
    return interactions


def robust_rank_aggregate(method_ranks: Mapping[str, Mapping[Any, int]]) -> Any:
    """Robust rank aggregation (Kolde) over per-method interaction rankings.

    For each interaction, normalised ranks across methods are compared against the beta
    order statistics expected under a uniform null; the smallest tail probability, corrected
    for the number of methods, is the aggregate score. An interaction ranked highly by every
    method scores far better than one ranked first by a single tool.
    """
    import pandas as pd
    from scipy.stats import beta

    if len(method_ranks) < 2:
        raise ValueError(
            "rank aggregation needs at least two methods; the point is to find what "
            "several scorers agree on"
        )
    methods = list(method_ranks)
    sizes = {m: max(len(method_ranks[m]), 1) for m in methods}
    all_items = sorted({item for m in methods for item in method_ranks[m]}, key=str)

    rows = []
    for item in all_items:
        norm = []
        for m in methods:
            r = method_ranks[m].get(item)
            # Absent from a method's list => worst possible normalised rank for that method.
            norm.append((r / sizes[m]) if r is not None else 1.0)
        norm_sorted = np.sort(np.asarray(norm, dtype=float))
        k = len(norm_sorted)
        probs = [beta.cdf(norm_sorted[i], i + 1, k - i) for i in range(k)]
        rho = float(min(probs))
        rows.append({
            "interaction": item,
            "rra_score": min(1.0, rho * k),          # Bonferroni over the order statistics
            "mean_normalised_rank": round(float(np.mean(norm)), 4),
            "n_methods_reporting": int(sum(1 for m in methods if item in method_ranks[m])),
            "n_methods": k,
        })
    frame = pd.DataFrame(rows).sort_values("rra_score").reset_index(drop=True)
    frame["consensus_rank"] = np.arange(1, len(frame) + 1)
    return frame


@register
class CellCellInteraction(AnalysisModule):
    name = "interactions"
    version = "0.1.0"
    phase = "5"
    description = "Spatially constrained ligand-receptor inference with cross-method rank consensus."
    evidence = Evidence.INFERRED
    replication_unit = ReplicationUnit.SECTION
    supported_platforms = frozenset({PlatformClass.SINGLE_CELL, PlatformClass.SPOT})
    optional_dependencies = ("squidpy", "liana", "cellphonedb", "commot")

    def run(
        self, ctx: ModuleContext, *, method_ranks: Mapping[str, Mapping[Any, int]] | None = None,
        interactions: Any = None, contacts: Any = None, **inputs: Any,
    ) -> ModuleResult:
        self.check_platform(ctx)
        if method_ranks is None:
            raise ModuleError(
                "interactions.run requires method_ranks={method: {interaction: rank}}. Run "
                "the scorers first (CellPhoneDB, CellChat, NicheNet, COMMOT), then pass "
                "their rankings for aggregation."
            )
        frame = robust_rank_aggregate(method_ranks)
        if interactions is not None and contacts is not None:
            constrained = spatial_constraint(interactions, contacts)
            supported = set(map(str, constrained.loc[
                constrained["spatially_supported"], "interaction"
            ])) if "interaction" in constrained.columns else None
            if supported is not None:
                frame["spatially_supported"] = frame["interaction"].astype(str).isin(supported)
        return self.result(
            ctx, frame,
            notes=[
                "A ranked ligand-receptor pair is co-expression in adjacent neighbourhoods. "
                "It is a hypothesis about signalling, not evidence of it.",
                f"Consensus over {len(method_ranks)} scorers by robust rank aggregation; "
                "single-method top hits are down-weighted by construction.",
            ],
            methods=list(method_ranks),
            n_interactions=int(len(frame)),
        )
