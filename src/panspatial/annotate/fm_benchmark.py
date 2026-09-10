"""The gate a foundation model must pass before it is allowed into the pipeline.

A 2025 Genome Biology zero-shot evaluation found scGPT and Geneformer performed worse than
scVI, Harmony and a plain highly-variable-gene baseline for cell-type clustering — in some
cases worse than the same architectures with random weights. Nicheformer, pretrained on
spatial data, does beat those baselines on spatial tasks.

The lesson is not "foundation models are bad", it is "measure on your own held-out patients
before adopting one". This module makes that measurement a precondition.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Mapping, Sequence

import numpy as np

log = logging.getLogger("panspatial.annotate.fm")

# Baselines a candidate must beat. HVG is included deliberately: it is the embarrassing one.
REQUIRED_BASELINES = ("hvg", "pca", "harmony", "scvi")


@dataclass
class BenchmarkVerdict:
    candidate: str
    metric: str
    candidate_score: float
    baseline_scores: dict[str, float]
    best_baseline: str
    margin: float
    admitted: bool
    reason: str
    missing_baselines: list[str] = field(default_factory=list)
    n_held_out_patients: int | None = None

    def to_dict(self) -> dict[str, object]:
        return {
            "candidate": self.candidate, "metric": self.metric,
            "candidate_score": self.candidate_score,
            "best_baseline": self.best_baseline,
            "best_baseline_score": self.baseline_scores.get(self.best_baseline),
            "margin": round(self.margin, 5), "admitted": self.admitted,
            "reason": self.reason,
            "n_held_out_patients": self.n_held_out_patients,
            **{f"baseline_{k}": v for k, v in self.baseline_scores.items()},
        }

    def __str__(self) -> str:  # pragma: no cover - display only
        verdict = "ADMITTED" if self.admitted else "REJECTED"
        return (f"{verdict}: {self.candidate} {self.metric}={self.candidate_score:.4f} vs "
                f"best baseline {self.best_baseline}={self.baseline_scores.get(self.best_baseline, float('nan')):.4f} "
                f"(margin {self.margin:+.4f})")


def must_beat_baselines(
    candidate: str,
    candidate_score: float,
    baseline_scores: Mapping[str, float],
    *,
    metric: str = "ARI",
    margin: float = 0.02,
    higher_is_better: bool = True,
    required_baselines: Sequence[str] = REQUIRED_BASELINES,
    n_held_out_patients: int | None = None,
    min_held_out_patients: int = 5,
) -> BenchmarkVerdict:
    """Admit a candidate embedding only if it beats every required baseline by ``margin``.

    The comparison must be on held-out patients: an embedding evaluated on the patients it
    was fitted to will win regardless of merit, so too few held-out patients is itself a
    rejection.
    """
    if not baseline_scores:
        raise ValueError("at least one baseline score is required")
    if not np.isfinite(candidate_score):
        raise ValueError(f"candidate score for {candidate!r} is not finite")

    scores = {str(k): float(v) for k, v in baseline_scores.items() if np.isfinite(v)}
    missing = [b for b in required_baselines if b not in scores]
    best = max(scores, key=scores.get) if higher_is_better else min(scores, key=scores.get)
    diff = (candidate_score - scores[best]) if higher_is_better else (scores[best] - candidate_score)

    if n_held_out_patients is not None and n_held_out_patients < min_held_out_patients:
        verdict = BenchmarkVerdict(
            candidate, metric, float(candidate_score), scores, best, diff, False,
            f"evaluated on only {n_held_out_patients} held-out patient(s); at least "
            f"{min_held_out_patients} are required for the comparison to mean anything",
            missing, n_held_out_patients,
        )
    elif missing:
        verdict = BenchmarkVerdict(
            candidate, metric, float(candidate_score), scores, best, diff, False,
            "required baseline(s) were not run: " + ", ".join(missing)
            + ". A candidate cannot be admitted against a partial comparison.",
            missing, n_held_out_patients,
        )
    elif diff >= margin:
        verdict = BenchmarkVerdict(
            candidate, metric, float(candidate_score), scores, best, diff, True,
            f"beats the best baseline ({best}) by {diff:.4f} {metric}, at or above the "
            f"required margin of {margin}",
            missing, n_held_out_patients,
        )
    else:
        verdict = BenchmarkVerdict(
            candidate, metric, float(candidate_score), scores, best, diff, False,
            f"does not beat {best} by the required margin ({diff:+.4f} < {margin}). "
            "Use the baseline: it is cheaper, faster and better understood.",
            missing, n_held_out_patients,
        )
    log.info("%s — %s", verdict, verdict.reason)
    return verdict


def benchmark_table(verdicts: Sequence[BenchmarkVerdict]):
    """Tidy frame of verdicts, for the methods section."""
    import pandas as pd

    return pd.DataFrame([v.to_dict() for v in verdicts])
