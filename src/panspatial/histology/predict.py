"""Extending reach with histology, without letting imputation become evidence.

H&E slides vastly outnumber spatial assays. HEST-1k pairs 1,108 ST samples with whole-slide
images — 1.5 million spots and over 60 million cells — and its benchmark evaluates 11
pathology foundation models on gene-expression prediction across nine organs and eight
cancer types. Performance varies widely by task and organ, which is the fact that governs
how these predictions may be used.

Everything this module produces is marked :attr:`Evidence.IMPUTED`. That is not a caveat in
prose; it is a property the framework checks, so a predicted profile cannot be used to
validate the spatial model it was derived from.
"""

from __future__ import annotations

import logging
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

log = logging.getLogger("panspatial.histology")

FOUNDATION_MODELS = {
    "uni": "Pathology foundation model (Nature Medicine)",
    "uni2-h": "Successor UNI checkpoint",
    "virchow": "Pathology foundation model",
    "virchow2": "Successor Virchow checkpoint",
    "conch": "Vision-language pathology model",
    "prov-gigapath": "Whole-slide foundation model",
    "stpath": "Generative model predicting spatial expression from WSIs",
    "omiclip": "Visual-omics model linking H&E and transcriptomics",
}


def per_organ_performance(
    predicted: Any,
    observed: Any,
    organs: Sequence[str],
    patients: Sequence[str],
    *,
    min_patients: int = 3,
    good_r: float = 0.4,
) -> Any:
    """Prediction accuracy per organ, on held-out patients, naming where it fails.

    Pooled metrics hide the organs where a model does not work. Reporting per organ, and
    listing the failures explicitly, is the difference between a usable extension and an
    unfalsifiable one.
    """
    import pandas as pd
    from scipy.stats import pearsonr

    predicted = np.asarray(predicted, dtype=float)
    observed = np.asarray(observed, dtype=float)
    organs = np.asarray(organs, dtype=object).astype(str)
    patients = np.asarray(patients, dtype=object).astype(str)
    if not (predicted.shape == observed.shape):
        raise ValueError(f"predicted {predicted.shape} and observed {observed.shape} must match")
    if len(organs) != predicted.shape[0] or len(patients) != predicted.shape[0]:
        raise ValueError("organs and patients must have one entry per spot")

    rows = []
    for organ in sorted(set(organs)):
        mask = organs == organ
        n_patients = len(set(patients[mask]))
        if n_patients < min_patients:
            rows.append({
                "organ": organ, "n_spots": int(mask.sum()), "n_patients": n_patients,
                "mean_gene_r": np.nan, "usable": False,
                "verdict": f"only {n_patients} held-out patient(s); not evaluable",
            })
            continue
        p, o = predicted[mask], observed[mask]
        rs = []
        for j in range(p.shape[1]):
            if np.std(p[:, j]) > 0 and np.std(o[:, j]) > 0:
                r = pearsonr(p[:, j], o[:, j])[0]
                if np.isfinite(r):
                    rs.append(r)
        mean_r = float(np.mean(rs)) if rs else float("nan")
        usable = bool(np.isfinite(mean_r) and mean_r >= good_r)
        rows.append({
            "organ": organ, "n_spots": int(mask.sum()), "n_patients": n_patients,
            "mean_gene_r": round(mean_r, 4) if np.isfinite(mean_r) else np.nan,
            "n_genes_scored": len(rs), "usable": usable,
            "verdict": (f"usable for cohort extension (mean gene r = {mean_r:.2f})" if usable
                        else f"prediction is poor here (mean gene r = {mean_r:.2f} < {good_r}); "
                             "exclude this organ from histology-extended analyses and say so"),
        })
    frame = pd.DataFrame(rows)
    failing = frame.loc[~frame["usable"], "organ"].tolist()
    if failing:
        log.warning("histology prediction unusable in: %s", ", ".join(failing))
    return frame


@register
class HistologyExtension(AnalysisModule):
    name = "histology"
    version = "0.1.0"
    phase = "5"
    description = "Foundation-model expression prediction from H&E, flagged as imputation."
    evidence = Evidence.IMPUTED
    replication_unit = ReplicationUnit.SPOT
    supported_platforms = frozenset({PlatformClass.SPOT, PlatformClass.SINGLE_CELL})
    optional_dependencies = ("torch", "timm", "hest")

    def run(
        self, ctx: ModuleContext, *, predicted: Any = None, observed: Any = None,
        organs: Sequence[str] | None = None, patients: Sequence[str] | None = None,
        **inputs: Any,
    ) -> ModuleResult:
        if predicted is None or observed is None or organs is None or patients is None:
            raise ModuleError(
                "histology.run requires predicted=, observed=, organs= and patients= for "
                "held-out evaluation. A prediction model with no held-out evaluation "
                "cannot be admitted to the pipeline."
            )
        frame = per_organ_performance(
            predicted, observed, organs, patients,
            min_patients=ctx.param("min_patients", 3),
            good_r=ctx.param("good_r", 0.4),
        )
        usable = frozenset(frame.loc[frame["usable"], "organ"]) if len(frame) else frozenset()
        return self.result(
            ctx, frame, reliable_entities=usable,
            notes=[
                "Expression predicted from morphology is an imputation. It may generate "
                "hypotheses and extend cohorts; it may never serve as the measurement that "
                "validates the model that produced it.",
                "Performance is reported per organ on held-out patients; organs where "
                "prediction fails are named rather than omitted.",
            ],
            model=ctx.param("model", "unspecified"),
        )
