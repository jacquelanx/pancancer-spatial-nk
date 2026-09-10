"""From a spatial niche to a ranked, falsifiable target list.

Everything upstream of the validation experiments is observational. A hazard ratio derived
from bulk deconvolution of a different cohort is a separate and weaker line of evidence
than the spatial observation it is invoked to support — not a validation of it. This module
keeps that distinction in the data structure: each candidate carries its evidence tier and
the experiment that would refute it, and no candidate can be ranked without one.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
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

log = logging.getLogger("panspatial.translate")

# Weights for the composite priority score. Spatial effect and independent genetic evidence
# dominate; tractability breaks ties rather than driving the ranking.
DEFAULT_WEIGHTS = {
    "spatial_effect": 0.30,
    "conservation": 0.25,
    "dependency": 0.20,
    "survival": 0.15,
    "tractability": 0.10,
}


def score_signature(expression: Any, signature: Sequence[str], *, method: str = "zscore") -> np.ndarray:
    """Score a gene signature per sample from a gene × sample expression matrix.

    Genes are z-scored across samples before averaging, so a single highly expressed gene
    cannot dominate the score.
    """
    import pandas as pd

    frame = pd.DataFrame(expression)
    present = [g for g in signature if g in frame.index]
    missing = sorted(set(signature) - set(present))
    if not present:
        raise ValueError("none of the signature genes are present in the expression matrix")
    if missing:
        log.warning("%d/%d signature genes absent: %s",
                    len(missing), len(signature), ", ".join(missing[:8]))
    sub = frame.loc[present].to_numpy(dtype=float)
    if method == "zscore":
        sd = sub.std(axis=1, keepdims=True)
        z = np.divide(sub - sub.mean(axis=1, keepdims=True), sd,
                      out=np.zeros_like(sub), where=sd > 0)
        return z.mean(axis=0)
    if method == "mean":
        return sub.mean(axis=0)
    raise ValueError(f"method must be zscore|mean, got {method!r}")


def survival_association(
    scores: Sequence[float],
    time: Sequence[float],
    event: Sequence[int],
    covariates: Any = None,
    *,
    required_covariates: Sequence[str] = ("stage",),
) -> dict[str, Any]:
    """Cox association between a signature score and outcome, adjustment enforced.

    An unadjusted hazard ratio for a tumour-microenvironment signature is largely a
    restatement of stage. Adjustment is required rather than recommended: omitting the
    required covariates raises.
    """
    import pandas as pd

    lifelines = None
    try:
        import lifelines  # noqa: F811
    except ImportError:
        pass

    frame = pd.DataFrame({"score": np.asarray(scores, dtype=float),
                          "time": np.asarray(time, dtype=float),
                          "event": np.asarray(event, dtype=int)})
    if covariates is not None:
        cov = pd.DataFrame(covariates).reset_index(drop=True)
        missing = [c for c in required_covariates if c not in cov.columns]
        if missing:
            raise ValueError(
                f"required covariate(s) absent: {missing}. An unadjusted hazard ratio for a "
                "microenvironment signature mostly restates stage; adjustment is not optional."
            )
        frame = pd.concat([frame, cov], axis=1)
    elif required_covariates:
        raise ValueError(
            f"covariates= is required and must include {list(required_covariates)}. "
            "Reporting an unadjusted association here would overstate the evidence."
        )
    frame = frame.dropna()
    if frame["event"].sum() < 10:
        raise ValueError(
            f"only {int(frame['event'].sum())} events; a Cox model on this few is not "
            "interpretable. Report the sample size and decline the analysis."
        )
    if lifelines is None:
        raise ModuleError(
            "survival analysis needs `lifelines`: pip install 'panspatial[translate]'"
        )

    model = lifelines.CoxPHFitter()
    model.fit(frame, duration_col="time", event_col="event")
    row = model.summary.loc["score"]
    return {
        "hazard_ratio": float(np.exp(row["coef"])),
        "ci_low": float(np.exp(row["coef lower 95%"])),
        "ci_high": float(np.exp(row["coef upper 95%"])),
        "p_value": float(row["p"]),
        "n": int(len(frame)),
        "n_events": int(frame["event"].sum()),
        "adjusted_for": [c for c in frame.columns if c not in {"score", "time", "event"}],
        "concordance": float(model.concordance_index_),
    }


@dataclass
class TargetCandidate:
    """One nominated target, with the experiment that would kill it."""

    gene: str
    niche: str
    spatial_effect: float          # standardised effect in the niche, 0-1 after scaling
    conservation: float            # fraction of cancer types where it recurs, 0-1
    dependency: float | None = None       # e.g. DepMap essentiality, 0-1
    survival: float | None = None         # scaled strength of adjusted association, 0-1
    tractability: float | None = None     # existing chemical matter / surface accessibility
    evidence_tier: str = "observational"
    falsifier: str = ""
    notes: list[str] = field(default_factory=list)

    def validate(self) -> None:
        if not self.falsifier.strip():
            raise ValueError(
                f"{self.gene}: a target cannot be nominated without a falsifier — the "
                "observation or experiment that would refute it."
            )
        for name in ("spatial_effect", "conservation", "dependency", "survival", "tractability"):
            v = getattr(self, name)
            if v is not None and not 0.0 <= v <= 1.0:
                raise ValueError(f"{self.gene}: {name} must be scaled to [0, 1], got {v}")


def prioritise_targets(
    candidates: Sequence[TargetCandidate],
    *,
    weights: Mapping[str, float] = DEFAULT_WEIGHTS,
    require_conservation: float = 0.0,
) -> Any:
    """Rank candidates by a weighted composite, with missing evidence penalised not ignored.

    A missing component scores zero rather than being dropped from the average: a candidate
    with no dependency data is genuinely weaker evidence than one with a strong dependency
    score, and averaging over present components only would hide that.
    """
    import pandas as pd

    if not candidates:
        raise ValueError("no candidates supplied")
    unknown = set(weights) - set(DEFAULT_WEIGHTS)
    if unknown:
        raise ValueError(f"unknown weight key(s): {sorted(unknown)}")
    total_w = sum(weights.values())
    if not np.isclose(total_w, 1.0, atol=1e-6):
        raise ValueError(f"weights must sum to 1.0, got {total_w:.4f}")

    rows = []
    for c in candidates:
        c.validate()
        if c.conservation < require_conservation:
            continue
        components = {
            "spatial_effect": c.spatial_effect,
            "conservation": c.conservation,
            "dependency": c.dependency,
            "survival": c.survival,
            "tractability": c.tractability,
        }
        missing = [k for k, v in components.items() if v is None]
        score = sum(weights.get(k, 0.0) * (v if v is not None else 0.0)
                    for k, v in components.items())
        rows.append({
            "gene": c.gene, "niche": c.niche,
            "priority_score": round(float(score), 4),
            **{k: v for k, v in components.items()},
            "evidence_tier": c.evidence_tier,
            "missing_evidence": ",".join(missing),
            "falsifier": c.falsifier,
            "notes": "; ".join(c.notes),
        })
    frame = pd.DataFrame(rows).sort_values("priority_score", ascending=False).reset_index(drop=True)
    frame["rank"] = np.arange(1, len(frame) + 1)
    log.info(
        "%d candidates ranked; %d carry complete evidence",
        len(frame), int((frame["missing_evidence"] == "").sum()),
    )
    return frame


@register
class Translation(AnalysisModule):
    name = "translate"
    version = "0.1.0"
    phase = "6"
    description = "Signature scoring, adjusted survival association, and falsifiable target ranking."
    evidence = Evidence.INFERRED
    replication_unit = ReplicationUnit.PATIENT
    supported_platforms = frozenset(
        {PlatformClass.SINGLE_CELL, PlatformClass.SPOT, PlatformClass.REGION,
         PlatformClass.DISSOCIATED}
    )
    optional_dependencies = ("lifelines", "scikit-survival")

    def run(
        self, ctx: ModuleContext, *, candidates: Sequence[TargetCandidate] | None = None,
        **inputs: Any,
    ) -> ModuleResult:
        if not candidates:
            raise ModuleError(
                "translate.run requires candidates=[TargetCandidate...]. Each must carry a "
                "falsifier; a target without one cannot be nominated."
            )
        frame = prioritise_targets(
            candidates,
            weights=ctx.param("weights", DEFAULT_WEIGHTS),
            require_conservation=ctx.param("require_conservation", 0.0),
        )
        return self.result(
            ctx, frame,
            notes=[
                "All evidence upstream of the validation experiments is observational; "
                "verbs in the output are matched to that.",
                "A bulk-deconvolution hazard ratio from a different cohort is a separate, "
                "weaker line of evidence — never a validation of the spatial observation.",
                f"{int((frame['missing_evidence'] != '').sum())} candidate(s) have "
                "incomplete evidence and are penalised, not excluded.",
            ],
            n_candidates=int(len(frame)),
        )
