"""Outcome models with leakage prevention built into the evaluation loop.

Two spots 100 µm apart share microenvironment, patient, batch and — on spot platforms —
sometimes the same cells. Put one in training and one in test and the model memorises. An
AUC of 0.95 obtained that way means nothing, and it is the default outcome of a random
split.

The protocol here is not advisory. Splits are made over patients; a spatial buffer is
applied where sections are split internally; preprocessing is fitted inside the fold; and
:func:`panspatial.stats.leakage.assert_no_leakage` runs on every fold, not once.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Callable, Sequence

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
from panspatial.stats.leakage import assert_no_leakage, patient_blocked_splits, spatial_buffer_filter

log = logging.getLogger("panspatial.predict")


@dataclass
class FoldResult:
    fold: int
    n_train: int
    n_test: int
    n_train_after_buffer: int
    score: float
    n_test_patients: int


@dataclass
class CVReport:
    metric: str
    folds: list[FoldResult] = field(default_factory=list)
    leakage_checks: int = 0
    buffer_um: float = 0.0

    @property
    def mean(self) -> float:
        return float(np.mean([f.score for f in self.folds])) if self.folds else float("nan")

    @property
    def sd(self) -> float:
        return float(np.std([f.score for f in self.folds], ddof=1)) if len(self.folds) > 1 else 0.0

    def table(self):
        import pandas as pd
        return pd.DataFrame([{
            "fold": f.fold, "n_train": f.n_train,
            "n_train_after_buffer": f.n_train_after_buffer,
            "n_test": f.n_test, "n_test_patients": f.n_test_patients,
            self.metric: round(f.score, 4),
        } for f in self.folds])

    def summary(self) -> str:
        return (f"{self.metric} = {self.mean:.4f} +/- {self.sd:.4f} over {len(self.folds)} "
                f"patient-level folds; {self.leakage_checks} leakage assertions passed")


def patient_level_cv(
    X: np.ndarray,
    y: np.ndarray,
    patients: Sequence[str],
    *,
    fit_predict: Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray],
    score: Callable[[np.ndarray, np.ndarray], float],
    metric: str = "auc",
    sections: Sequence[str] | None = None,
    coords: np.ndarray | None = None,
    buffer_um: float = 0.0,
    n_splits: int = 5,
    stratify_by: Sequence[str] | None = None,
    seed: int = 0,
) -> CVReport:
    """Cross-validate with patient-level folds, a spatial buffer, and per-fold assertions.

    ``fit_predict`` receives ``(X_train, y_train, X_test)`` and must do *all* fitting —
    scaling, feature selection, imputation — inside itself. Anything fitted before this
    call has already seen the test fold.
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y)
    patients = np.asarray(patients, dtype=object).astype(str)
    if not (len(X) == len(y) == len(patients)):
        raise ValueError("X, y and patients must describe the same observations")
    if buffer_um > 0 and (coords is None or sections is None):
        raise ValueError("a spatial buffer needs coords= and sections=")

    report = CVReport(metric=metric, buffer_um=buffer_um)
    splits = patient_blocked_splits(
        patients, n_splits=n_splits, rng=np.random.default_rng(seed), stratify_by=stratify_by
    )
    for i, (train, test) in enumerate(splits):
        assert_no_leakage(train, test, patients)
        report.leakage_checks += 1
        n_train_before = len(train)
        if buffer_um > 0:
            train = spatial_buffer_filter(train, test, coords, sections, buffer_um=buffer_um)
            assert_no_leakage(
                train, test, patients, coords=coords, section_ids=sections,
                buffer_um=buffer_um,
            )
            report.leakage_checks += 1
        if len(train) == 0:
            raise ModuleError(
                f"fold {i}: the spatial buffer removed every training observation. The "
                "buffer is too large relative to the section, or the split is degenerate."
            )
        preds = fit_predict(X[train], y[train], X[test])
        report.folds.append(FoldResult(
            fold=i, n_train=n_train_before, n_test=len(test),
            n_train_after_buffer=len(train), score=float(score(y[test], preds)),
            n_test_patients=len(set(patients[test])),
        ))
    log.info(report.summary())
    return report


@register
class OutcomePrediction(AnalysisModule):
    name = "predict"
    version = "0.1.0"
    phase = "6"
    description = "Niche-to-outcome models evaluated under patient-level, buffered CV."
    evidence = Evidence.INFERRED
    replication_unit = ReplicationUnit.PATIENT
    supported_platforms = frozenset({PlatformClass.SINGLE_CELL, PlatformClass.SPOT})
    optional_dependencies = ("scikit-learn", "torch", "torch-geometric", "lifelines")

    def run(
        self, ctx: ModuleContext, *, X: Any = None, y: Any = None,
        patients: Sequence[str] | None = None, fit_predict: Callable | None = None,
        score: Callable | None = None, **inputs: Any,
    ) -> ModuleResult:
        if X is None or y is None or patients is None or fit_predict is None or score is None:
            raise ModuleError(
                "predict.run requires X=, y=, patients=, fit_predict= and score=. "
                "fit_predict must perform all preprocessing internally: anything fitted "
                "before the split has already seen the test fold."
            )
        report = patient_level_cv(
            X, y, patients, fit_predict=fit_predict, score=score,
            metric=ctx.param("metric", "auc"),
            sections=inputs.get("sections"), coords=inputs.get("coords"),
            buffer_um=ctx.param("buffer_um", 0.0),
            n_splits=ctx.param("n_splits", 5),
            stratify_by=inputs.get("stratify_by"), seed=ctx.seed,
        )
        return self.result(
            ctx, report.table(),
            notes=[
                report.summary(),
                "Folds are patient-level with a spatial buffer; leakage is asserted inside "
                "the loop, on every fold.",
                "Final performance must additionally be confirmed on a cohort held out "
                "entirely and never touched during development.",
            ],
            mean_score=report.mean, sd_score=report.sd,
            leakage_checks=report.leakage_checks,
        )
