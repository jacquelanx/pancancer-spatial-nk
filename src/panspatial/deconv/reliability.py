"""The gate that decides which cell types a spot-based study may draw conclusions about.

Deconvolution of a rare cell type is the weak link in every spot-based pan-cancer study.
NK cells are 0.5–3% of a tumour and a 55 µm Visium spot holds 1–10 cells, so an "NK weight"
is a similarity score against a reference profile whose nearest competitor is CD8 T.

Rather than assume reliability from a published benchmark run on other tissue, this module
measures it on *this* cohort's reference: spots of known composition are simulated by
summing counts from a known number of reference cells, so the ground truth is exact.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Sequence

import numpy as np

log = logging.getLogger("panspatial.deconv.reliability")


def simulate_mixtures(
    counts: np.ndarray,
    labels: Sequence[str],
    *,
    n_spots: int = 500,
    cells_per_spot: int = 8,
    rng: np.random.Generator | None = None,
    rare_boost: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Build synthetic spots of known composition from a single-cell reference.

    Args:
        counts: (n_genes, n_cells) reference count matrix.
        labels: cell-type label per reference cell.
        cells_per_spot: cells pooled into each synthetic spot; 8 approximates a Visium spot.
        rare_boost: in ``[0, 1]``, mixes uniform-over-types sampling into the natural
            abundance sampling. Zero reproduces the cohort's real composition; raising it
            guarantees rare types appear often enough to be evaluated at all, which is the
            regime the gate is about.

    Returns ``(sim_counts (n_genes, n_spots), truth (n_spots, n_types), type_names)``.
    """
    counts = np.asarray(counts)
    labels = np.asarray(labels)
    if counts.shape[1] != len(labels):
        raise ValueError(
            f"counts has {counts.shape[1]} cells but {len(labels)} labels were given"
        )
    if not 0.0 <= rare_boost <= 1.0:
        raise ValueError(f"rare_boost must be in [0, 1], got {rare_boost}")
    if cells_per_spot < 1:
        raise ValueError("cells_per_spot must be at least 1")
    rng = rng or np.random.default_rng(0)

    types = sorted(set(map(str, labels)))
    n_types = len(types)
    type_index = {t: i for i, t in enumerate(types)}
    by_type = {t: np.flatnonzero(labels.astype(str) == t) for t in types}

    natural = np.array([len(by_type[t]) for t in types], dtype=float)
    natural /= natural.sum()
    uniform = np.full(n_types, 1.0 / n_types)
    p_type = (1.0 - rare_boost) * natural + rare_boost * uniform

    sim = np.zeros((counts.shape[0], n_spots), dtype=float)
    truth = np.zeros((n_spots, n_types), dtype=float)
    for s in range(n_spots):
        chosen_types = rng.choice(n_types, size=cells_per_spot, p=p_type)
        for ti in chosen_types:
            pool = by_type[types[ti]]
            cell = pool[rng.integers(len(pool))]
            col = counts[:, cell]
            sim[:, s] += np.asarray(col.todense()).ravel() if hasattr(col, "todense") else np.asarray(col).ravel()
            truth[s, ti] += 1.0
    truth /= cells_per_spot
    log.info(
        "simulated %d spots x %d cells from %d reference cells across %d types "
        "(rare_boost=%.2f)", n_spots, cells_per_spot, counts.shape[1], n_types, rare_boost,
    )
    return sim, truth, types


@dataclass
class ReliabilityRow:
    cell_type: str
    true_mean_proportion: float
    pearson_r: float
    spearman_r: float
    rmse: float
    n_spots_present: int
    reliable: bool
    reason: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "cell_type": self.cell_type,
            "true_mean_proportion": round(self.true_mean_proportion, 5),
            "pearson_r": round(self.pearson_r, 4),
            "spearman_r": round(self.spearman_r, 4),
            "rmse": round(self.rmse, 5),
            "n_spots_present": self.n_spots_present,
            "reliable": self.reliable,
            "reason": self.reason,
        }


@dataclass
class ReliabilityReport:
    rows: list[ReliabilityRow] = field(default_factory=list)
    min_r: float = 0.7
    min_spots: int = 20

    @property
    def reliable(self) -> frozenset[str]:
        return frozenset(r.cell_type for r in self.rows if r.reliable)

    @property
    def unreliable(self) -> list[str]:
        return sorted(r.cell_type for r in self.rows if not r.reliable)

    def table(self):
        import pandas as pd
        return pd.DataFrame([r.to_dict() for r in self.rows]).sort_values(
            "pearson_r", ascending=False
        ).reset_index(drop=True)

    def summary(self) -> str:
        return (
            f"{len(self.reliable)}/{len(self.rows)} cell types passed "
            f"(Pearson r >= {self.min_r}); deconvolution-limited: "
            f"{', '.join(self.unreliable) or 'none'}"
        )


def evaluate_deconvolution(
    truth: np.ndarray,
    estimate: np.ndarray,
    cell_types: Sequence[str],
    *,
    min_r: float = 0.7,
    min_spots_present: int = 20,
) -> ReliabilityReport:
    """Score estimated against known proportions, per cell type.

    A type absent from nearly every simulated spot cannot be evaluated, and is failed with
    that reason rather than given a meaningless correlation on three points.
    """
    from scipy.stats import pearsonr, spearmanr

    truth = np.asarray(truth, dtype=float)
    estimate = np.asarray(estimate, dtype=float)
    if truth.shape != estimate.shape:
        raise ValueError(f"truth {truth.shape} and estimate {estimate.shape} must match")
    if truth.shape[1] != len(cell_types):
        raise ValueError(
            f"{truth.shape[1]} columns but {len(cell_types)} cell types were given"
        )

    report = ReliabilityReport(min_r=min_r, min_spots=min_spots_present)
    for i, ct in enumerate(cell_types):
        t, e = truth[:, i], estimate[:, i]
        present = int(np.sum(t > 0))
        rmse = float(np.sqrt(np.mean((t - e) ** 2)))
        if present < min_spots_present:
            report.rows.append(ReliabilityRow(
                str(ct), float(t.mean()), 0.0, 0.0, rmse, present, False,
                f"present in only {present} simulated spots (< {min_spots_present}); "
                "too rare in this reference to evaluate, so not evaluable and not usable",
            ))
            continue
        if np.std(t) == 0 or np.std(e) == 0:
            report.rows.append(ReliabilityRow(
                str(ct), float(t.mean()), 0.0, 0.0, rmse, present, False,
                "no variance in truth or estimate; correlation undefined",
            ))
            continue
        r = float(pearsonr(t, e)[0])
        rho = float(spearmanr(t, e)[0])
        ok = r >= min_r
        report.rows.append(ReliabilityRow(
            str(ct), float(t.mean()), r, rho, rmse, present, ok,
            f"Pearson r = {r:.2f} against known mixture proportions"
            + ("" if ok else f", below the {min_r} threshold; weights are similarity scores, "
                             "not abundances, and may not carry a conclusion"),
        ))
    log.info(report.summary())
    if report.unreliable:
        log.warning(
            "deconvolution-limited cell types on this reference: %s",
            ", ".join(report.unreliable),
        )
    return report
