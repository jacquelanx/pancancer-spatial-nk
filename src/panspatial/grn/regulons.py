"""Transcription-factor regulons, projected onto tissue and tested for conservation.

BEELINE (Nature Methods 2020) evaluated GRN inference on synthetic and experimental
single-cell data and found heterogeneous performance with moderate AUPRC and moderate early
precision; methods not requiring pseudotime-ordered cells were generally more accurate. A
regulon is therefore a prioritised hypothesis, not a validated network, and this module is
built to say so.

Two things a pan-cancer cohort can add that a single study cannot: whether a regulon
recurs across cancer types at patient level, and whether its activity actually localises to
the niche it is claimed to drive.
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
from panspatial.stats.spatial_stats import benjamini_hochberg, random_effects_meta

log = logging.getLogger("panspatial.grn")


def regulon_niche_enrichment(
    activity: Any,
    niche: Sequence[str],
    patient: Sequence[str],
    *,
    target_niche: str,
    min_cells_per_patient: int = 20,
) -> Any:
    """Per-patient enrichment of each regulon's activity inside one niche, then pooled.

    Effect size is the standardised mean difference (Hedges-corrected Cohen's d) between
    cells inside and outside the niche, computed per patient and combined by random-effects
    meta-analysis. Pooling cells across patients first would make every regulon significant.
    """
    import pandas as pd

    frame = pd.DataFrame(activity)
    niche = np.asarray(niche, dtype=object).astype(str)
    patient = np.asarray(patient, dtype=object).astype(str)
    if not (len(frame) == len(niche) == len(patient)):
        raise ValueError("activity, niche and patient must describe the same cells")
    if target_niche not in set(niche):
        raise ValueError(f"target niche {target_niche!r} does not occur in the niche labels")

    rows = []
    for regulon in frame.columns:
        values = frame[regulon].to_numpy(dtype=float)
        effects, variances, labels = [], [], []
        for p in sorted(set(patient)):
            in_p = patient == p
            inside = values[in_p & (niche == target_niche)]
            outside = values[in_p & (niche != target_niche)]
            if len(inside) < min_cells_per_patient or len(outside) < min_cells_per_patient:
                continue
            n1, n2 = len(inside), len(outside)
            pooled_sd = np.sqrt(
                ((n1 - 1) * np.var(inside, ddof=1) + (n2 - 1) * np.var(outside, ddof=1))
                / (n1 + n2 - 2)
            )
            if pooled_sd == 0:
                continue
            d = (inside.mean() - outside.mean()) / pooled_sd
            j = 1.0 - 3.0 / (4 * (n1 + n2) - 9)      # Hedges' small-sample correction
            g = d * j
            var = (n1 + n2) / (n1 * n2) + g**2 / (2 * (n1 + n2))
            effects.append(float(g))
            variances.append(float(var))
            labels.append(p)
        if len(effects) < 2:
            continue
        meta = random_effects_meta(effects, variances, labels=labels)
        rows.append({
            "regulon": str(regulon), "niche": target_niche,
            "hedges_g": round(meta.estimate, 4),
            "ci_low": round(meta.ci_low, 4), "ci_high": round(meta.ci_high, 4),
            "p_value": meta.p_value, "n_patients": meta.k, "i2": round(meta.i2, 3),
            "loo_min": round(meta.loo_min, 4), "loo_max": round(meta.loo_max, 4),
            "single_patient_driven": meta.is_driven_by_one_cohort,
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr"] = benjamini_hochberg(out["p_value"].to_numpy())
        out = out.sort_values("p_value").reset_index(drop=True)
    log.info("%d regulons tested for enrichment in niche %r", len(out), target_niche)
    return out


def conserved_regulons(
    per_cancer: Mapping[str, Any],
    *,
    min_cancer_types: int = 3,
    fdr_threshold: float = 0.05,
    require_same_direction: bool = True,
) -> Any:
    """Regulons enriched in the same niche, in the same direction, across cancer types."""
    import pandas as pd

    if len(per_cancer) < 2:
        raise ValueError("conservation needs results from at least two cancer types")

    records: dict[str, list[tuple[str, float]]] = {}
    for cancer, frame in per_cancer.items():
        frame = pd.DataFrame(frame)
        if frame.empty:
            continue
        hits = frame[frame.get("fdr", frame["p_value"]) < fdr_threshold]
        for _, r in hits.iterrows():
            records.setdefault(str(r["regulon"]), []).append((str(cancer), float(r["hedges_g"])))

    rows = []
    for regulon, hits in records.items():
        cancers = [c for c, _ in hits]
        effects = np.array([g for _, g in hits])
        same_dir = bool(np.all(effects > 0) or np.all(effects < 0))
        conserved = len(set(cancers)) >= min_cancer_types and (same_dir or not require_same_direction)
        rows.append({
            "regulon": regulon,
            "n_cancer_types": len(set(cancers)),
            "cancer_types": ",".join(sorted(set(cancers))),
            "mean_hedges_g": round(float(effects.mean()), 4),
            "consistent_direction": same_dir,
            "conserved": conserved,
            "verdict": (
                f"conserved across {len(set(cancers))} cancer types"
                if conserved else
                (f"direction flips across cancer types ({', '.join(sorted(set(cancers)))})"
                 if not same_dir else
                 f"only {len(set(cancers))} cancer type(s); report as cancer-type-specific")
            ),
        })
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["conserved", "n_cancer_types"], ascending=False).reset_index(drop=True)
    return out


@register
class RegulatoryNetworks(AnalysisModule):
    name = "grn"
    version = "0.1.0"
    phase = "5"
    description = "SCENIC regulon activity projected onto niches, tested for conservation."
    evidence = Evidence.INFERRED
    replication_unit = ReplicationUnit.PATIENT
    supported_platforms = frozenset(
        {PlatformClass.SINGLE_CELL, PlatformClass.SPOT, PlatformClass.DISSOCIATED}
    )
    optional_dependencies = ("pyscenic", "arboreto", "ctxcore")

    def run(
        self, ctx: ModuleContext, *, activity: Any = None, niche: Sequence[str] | None = None,
        patient: Sequence[str] | None = None, **inputs: Any,
    ) -> ModuleResult:
        if activity is None or niche is None or patient is None:
            raise ModuleError(
                "grn.run requires activity= (cells x regulons AUCell matrix), niche= and "
                "patient=. Run pySCENIC/SCENIC+ first; this module handles the spatial "
                "projection and the patient-level statistics."
            )
        frame = regulon_niche_enrichment(
            activity, niche, patient,
            target_niche=ctx.param("target_niche"),
            min_cells_per_patient=ctx.param("min_cells_per_patient", 20),
        )
        significant = frozenset(
            frame.loc[(frame["fdr"] < 0.05) & (~frame["single_patient_driven"]), "regulon"]
        ) if not frame.empty else frozenset()
        return self.result(
            ctx, frame, reliable_entities=significant,
            notes=[
                "Regulons are prioritised hypotheses for perturbation, not a validated "
                "network: GRN inference benchmarks report only moderate AUPRC.",
                "Enrichment is computed per patient and pooled by random-effects "
                "meta-analysis; regulons driven by a single patient are excluded.",
            ],
            target_niche=ctx.param("target_niche"),
        )
