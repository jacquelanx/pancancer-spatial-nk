"""NK-cell identification, state resolution, and spatial niche testing.

Three things about NK cells make this harder than the equivalent analysis for T cells or
macrophages, and each is handled explicitly rather than assumed away:

1. **Identity is a discrimination problem, not a marker problem.** ``FCGR3A`` is also on
   CD16+ monocytes, ``NCAM1`` is on neurons and some tumors, and the transcriptomes of NK
   and CD8 effector T cells are close enough that the separating evidence is the *absence*
   of ``CD3D/E/G`` and ``TRAC`` -- an absence that is unreliable at low depth. ``KLRF1`` and
   ``NCR1`` carry most of the positive discrimination. NKT, gamma-delta T, and ILC1 sit in
   the gap and are reported, not silently absorbed into "NK".
2. **NK cells are rare.** At 0.5-3% of a tumor, a section with 20,000 cells yields a few
   hundred NK cells and a per-section spatial test on a handful of them is noise. Sections
   below a minimum count are dropped and *counted as dropped*.
3. **Spot platforms cannot resolve them.** A 55 um Visium spot holds 1-10 cells; an NK
   deconvolution weight there is a similarity score against a reference profile that is
   itself close to CD8 T. :func:`require_single_cell_resolution` gates NK-resolved claims
   to Xenium/CosMx/MERFISH; Visium contributes niche context only.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.spatial import cKDTree

from panspatial.ingest.geo_sra import st_resolution
from panspatial.stats.spatial_stats import (
    MetaAnalysisResult,
    mean_nn_distance,
    random_effects_meta,
    torus_shift_permutation_test,
)

log = logging.getLogger("panspatial.nk")

# Positive identity. KLRF1/NCR1 are the discriminating pair; the rest are corroborating.
NK_POSITIVE_MARKERS: dict[str, list[str]] = {
    "discriminating": ["KLRF1", "NCR1"],
    "corroborating": ["NCAM1", "KLRD1", "NKG7", "GNLY", "PRF1", "CD160", "FCGR3A", "KLRB1"],
}

# Absence required. Each entry names the population it excludes.
NK_NEGATIVE_MARKERS: dict[str, list[str]] = {
    "T cell": ["CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2"],
    "gamma-delta T": ["TRDC", "TRGC1", "TRGC2"],
    "myeloid": ["CD14", "LYZ", "ITGAM", "CSF1R"],
    "B cell": ["MS4A1", "CD79A"],
}

# Functional states. Scored as programs -- a single marker does not define a state.
NK_STATE_PROGRAMS: dict[str, list[str]] = {
    "cytotoxic": ["PRF1", "GZMB", "GZMH", "FGFBP2", "FCGR3A", "SPON2"],
    "dysfunctional": ["HAVCR2", "LAG3", "TIGIT", "KLRC1", "CD96", "ENTPD1"],
    "resident": ["ITGA1", "CD69", "ZNF683", "CXCR6", "ITGAE"],
    "cytokine": ["IFNG", "TNF", "XCL1", "XCL2", "CSF2", "IL2RA"],
}

# ILC1 vs NK: ILC1 lack the cytotoxic program while sharing surface identity.
ILC1_DISCRIMINANTS: list[str] = ["EOMES", "TBX21", "PRF1", "GZMB"]

MIN_NK_PER_SECTION = 30


class DeconvolutionWarning(UserWarning):
    """Raised as an error when a spot-resolution platform is used for an NK-resolved claim."""


def require_single_cell_resolution(platform: str, *, claim: str = "NK-resolved spatial") -> None:
    """Refuse an NK-resolved claim on a spot-based platform.

    Visium and Slide-seq data remain useful here for niche-level context -- stromal
    architecture, TLS location, CAF density -- but an "NK cell is 40 um from a CAF" claim
    requires segmented cells.
    """
    resolution, single_cell = st_resolution(platform)
    if not single_cell:
        detail = f"{resolution:.0f} um" if resolution else "region-level"
        raise DeconvolutionWarning(
            f"{claim} claims are not supportable on {platform} ({detail} resolution). "
            f"An NK signal there is a deconvolution weight, not a cell, and at 0.5-3% "
            f"abundance it is not separable from CD8 T. Use Xenium/CosMx/MERFISH for the "
            f"claim and {platform} for niche context only."
        )


# ------------------------------------------------------------------- identity (AnnData)


@dataclass
class IdentityReport:
    """Audit trail for NK gating. Every count here belongs in a supplementary table."""

    n_input: int
    n_positive: int
    n_dropped_by: dict[str, int]
    n_ambiguous: int
    n_nk: int
    missing_markers: list[str]
    nk_fraction_by_sample: dict[str, float]

    def summary(self) -> str:
        drops = ", ".join(f"{k}={v}" for k, v in sorted(self.n_dropped_by.items()))
        return (
            f"{self.n_nk}/{self.n_input} cells called NK "
            f"({100 * self.n_nk / max(self.n_input, 1):.2f}%); "
            f"positive={self.n_positive}, dropped[{drops}], ambiguous={self.n_ambiguous}"
        )


def _score(adata, genes: Sequence[str], name: str, *, seed: int = 0):
    """Score a gene program, tolerating genes absent from the panel."""
    import scanpy as sc

    present = [g for g in genes if g in adata.var_names]
    if not present:
        log.warning("no genes from program %r present; scoring as zero", name)
        adata.obs[name] = 0.0
        return []
    sc.tl.score_genes(adata, present, score_name=name, random_state=seed)
    return present


def nk_identity_report(
    adata,
    *,
    sample_key: str = "sample_id",
    positive_threshold: float = 0.0,
    negative_threshold: float = 0.0,
    seed: int = 0,
) -> IdentityReport:
    """Score NK identity and report exactly what was gated out and why.

    Cells scoring positive for both NK identity and an exclusion program are counted as
    ambiguous rather than assigned: at typical scRNA-seq depth a ``CD3E``-positive,
    ``KLRF1``-positive barcode is more often an NKT cell or a doublet than an NK cell.
    """
    n_input = adata.n_obs
    all_positive = NK_POSITIVE_MARKERS["discriminating"] + NK_POSITIVE_MARKERS["corroborating"]
    found = _score(adata, all_positive, "nk_score", seed=seed)
    missing = sorted(set(all_positive) - set(found))
    if not set(NK_POSITIVE_MARKERS["discriminating"]) & set(found):
        log.warning(
            "neither KLRF1 nor NCR1 is in the panel; NK calls will rest on markers shared "
            "with CD8 T and CD16+ monocytes and should be treated as provisional"
        )

    dropped: dict[str, int] = {}
    exclusion = np.zeros(n_input, dtype=bool)
    for population, genes in NK_NEGATIVE_MARKERS.items():
        key = f"exclude_{population.replace(' ', '_').replace('-', '_')}"
        if _score(adata, genes, key, seed=seed):
            hit = np.asarray(adata.obs[key] > negative_threshold)
            dropped[population] = int(hit.sum())
            exclusion |= hit

    positive = np.asarray(adata.obs["nk_score"] > positive_threshold)
    ambiguous = positive & exclusion
    is_nk = positive & ~exclusion
    adata.obs["is_nk"] = is_nk
    adata.obs["nk_ambiguous"] = ambiguous

    by_sample: dict[str, float] = {}
    if sample_key in adata.obs:
        samples = np.asarray(adata.obs[sample_key])
        for s in np.unique(samples):
            mask = samples == s
            by_sample[str(s)] = float(is_nk[mask].mean())

    report = IdentityReport(
        n_input=n_input,
        n_positive=int(positive.sum()),
        n_dropped_by=dropped,
        n_ambiguous=int(ambiguous.sum()),
        n_nk=int(is_nk.sum()),
        missing_markers=missing,
        nk_fraction_by_sample=by_sample,
    )
    log.info(report.summary())
    return report


def extract_nk_cells(adata, *, sample_key: str = "sample_id", seed: int = 0):
    """Return the NK subset plus its identity report, as ``(adata_nk, report)``."""
    report = nk_identity_report(adata, sample_key=sample_key, seed=seed)
    if report.n_nk == 0:
        raise ValueError("no cells passed NK gating; check that the panel includes NK markers")
    return adata[adata.obs["is_nk"]].copy(), report


def score_nk_states(
    adata,
    *,
    resolution: float = 0.6,
    min_margin: float = 0.1,
    seed: int = 0,
    n_bootstrap: int = 20,
):
    """Assign NK states by program score, leaving genuinely ambiguous cells unassigned.

    Cells whose top two program scores differ by less than ``min_margin`` are labelled
    ``unassigned``. Forcing every cell into a state manufactures the clean four-state
    structure that reviewers correctly distrust. Cluster stability is reported as the mean
    bootstrap adjusted Rand index against the full-data Leiden partition; a state count
    that is not stable to resampling should not be reported as a discovery.
    """
    import scanpy as sc
    from sklearn.metrics import adjusted_rand_score

    for state, genes in NK_STATE_PROGRAMS.items():
        _score(adata, genes, f"program_{state}", seed=seed)

    states = list(NK_STATE_PROGRAMS)
    scores = np.column_stack([np.asarray(adata.obs[f"program_{s}"]) for s in states])
    order = np.argsort(-scores, axis=1)
    top, second = order[:, 0], order[:, 1]
    margin = scores[np.arange(len(scores)), top] - scores[np.arange(len(scores)), second]
    assigned = np.array([states[i] for i in top], dtype=object)
    assigned[margin < min_margin] = "unassigned"
    adata.obs["nk_state"] = assigned
    adata.obs["nk_state_margin"] = margin

    sc.pp.neighbors(adata, random_state=seed)
    sc.tl.leiden(adata, resolution=resolution, random_state=seed, key_added="nk_leiden")
    reference = np.asarray(adata.obs["nk_leiden"])

    rng = np.random.default_rng(seed)
    aris = []
    for b in range(n_bootstrap):
        idx = rng.choice(adata.n_obs, adata.n_obs, replace=True)
        sub = adata[idx].copy()
        sc.pp.neighbors(sub, random_state=seed + b + 1)
        sc.tl.leiden(sub, resolution=resolution, random_state=seed + b + 1, key_added="boot")
        aris.append(adjusted_rand_score(reference[idx], np.asarray(sub.obs["boot"])))
    stability = float(np.mean(aris))

    counts = {s: int((assigned == s).sum()) for s in states + ["unassigned"]}
    log.info("NK states %s | bootstrap ARI %.2f over %d resamples", counts, stability, n_bootstrap)
    if stability < 0.6:
        log.warning(
            "cluster stability ARI=%.2f is low; report the program scores as a continuum "
            "rather than as discrete states",
            stability,
        )
    return adata, {"counts": counts, "bootstrap_ari": stability}


# --------------------------------------------------------------------- spatial (arrays)


def distance_to_anchor(
    coords: np.ndarray,
    labels: Sequence[str],
    section_ids: Sequence[str],
    *,
    focal_label: str,
    anchor_label: str,
) -> dict[str, np.ndarray]:
    """Per-section nearest-anchor distances for focal cells, in the coordinate unit given.

    Distances are computed strictly within a section: coordinates from different sections
    are in different frames, so a cross-section distance is meaningless. Sections lacking
    either population are omitted from the result and logged.
    """
    coords = np.asarray(coords, dtype=float)
    labels = np.asarray(labels)
    section_ids = np.asarray(section_ids)
    if not (len(coords) == len(labels) == len(section_ids)):
        raise ValueError(
            f"coords ({len(coords)}), labels ({len(labels)}) and section_ids "
            f"({len(section_ids)}) must describe the same observations"
        )

    out: dict[str, np.ndarray] = {}
    for section in np.unique(section_ids):
        in_section = section_ids == section
        focal = coords[in_section & (labels == focal_label)]
        anchor = coords[in_section & (labels == anchor_label)]
        if len(focal) == 0 or len(anchor) == 0:
            log.info(
                "section %s: skipped (%d %s, %d %s)",
                section, len(focal), focal_label, len(anchor), anchor_label,
            )
            continue
        distances, _ = cKDTree(anchor).query(focal, k=1)
        out[str(section)] = distances
    return out


def per_section_niche_test(
    coords: np.ndarray,
    labels: Sequence[str],
    section_ids: Sequence[str],
    patient_ids: Sequence[str],
    *,
    focal_label: str,
    anchor_label: str,
    n_perms: int = 999,
    alternative: str = "less",
    min_focal: int = MIN_NK_PER_SECTION,
    seed: int = 0,
):
    """Run the torus-shift association test per section and return a tidy DataFrame.

    One row per section, carrying the observed statistic, the section's own null, and the
    patient it belongs to -- which is what :func:`cohort_meta_analysis` needs to pool at the
    correct replication level. Sections with too few focal cells are excluded and the
    exclusion is returned in the frame's ``attrs`` so it can be reported.
    """
    import pandas as pd

    coords = np.asarray(coords, dtype=float)
    labels = np.asarray(labels)
    section_ids = np.asarray(section_ids)
    patient_ids = np.asarray(patient_ids)

    rows, skipped = [], []
    for i, section in enumerate(np.unique(section_ids)):
        in_section = section_ids == section
        focal = coords[in_section & (labels == focal_label)]
        anchor = coords[in_section & (labels == anchor_label)]
        patients = np.unique(patient_ids[in_section])
        if len(focal) < min_focal or len(anchor) < 3:
            skipped.append(
                {"section": str(section), "n_focal": len(focal), "n_anchor": len(anchor)}
            )
            continue
        result = torus_shift_permutation_test(
            focal,
            anchor,
            statistic=mean_nn_distance,
            alternative=alternative,
            n_perms=n_perms,
            rng=np.random.default_rng(seed + i),
        )
        # log ratio of observed to null distance: negative means closer than chance.
        effect = float(np.log(result.observed / result.null_mean))
        # Variance of that log ratio under the section's own null, by the delta method.
        # The permutation SD already describes the section-level statistic (a mean over
        # this section's focal cells), so it must NOT be divided by the cell count again --
        # doing so treats cells as independent replicates and shrinks the within-section
        # variance towards zero, which pushes all the spread into tau^2 and inflates I^2.
        variance = float((result.null_sd / result.null_mean) ** 2)
        rows.append(
            {
                "section": str(section),
                "patient": str(patients[0]) if len(patients) == 1 else "|".join(map(str, patients)),
                "n_focal": result.n_focal,
                "n_anchor": result.n_anchor,
                "observed_um": result.observed,
                "null_mean_um": result.null_mean,
                "null_sd_um": result.null_sd,
                "log_ratio": effect,
                "variance": variance,
                "z": result.z,
                "p_value": result.p_value,
            }
        )

    frame = pd.DataFrame(rows)
    frame.attrs["skipped_sections"] = skipped
    frame.attrs["focal_label"] = focal_label
    frame.attrs["anchor_label"] = anchor_label
    frame.attrs["min_focal"] = min_focal
    log.info(
        "%s vs %s: %d sections tested, %d skipped for insufficient cells",
        focal_label, anchor_label, len(frame), len(skipped),
    )
    return frame


def cohort_meta_analysis(section_frame, *, group_key: str = "patient") -> MetaAnalysisResult:
    """Pool per-section effects to the patient level, then across patients.

    Sections from one patient are averaged first (inverse-variance within patient), so a
    patient contributing three sections does not count three times. The pooled estimate is
    then a random-effects meta-analysis over patients, which is the n that may be quoted.
    """
    if len(section_frame) < 2:
        raise ValueError(
            f"{len(section_frame)} section(s) with usable data; a pan-cancer claim needs "
            "several patients, and with fewer than two there is no between-patient variance "
            "to estimate"
        )
    grouped = section_frame.groupby(group_key)
    effects, variances, labels = [], [], []
    for name, block in grouped:
        w = 1.0 / block["variance"].to_numpy()
        effects.append(float(np.sum(w * block["log_ratio"].to_numpy()) / np.sum(w)))
        variances.append(float(1.0 / np.sum(w)))
        labels.append(str(name))
    result = random_effects_meta(effects, variances, labels=labels)
    log.info(
        "%s vs %s pooled over %d %ss: %s",
        section_frame.attrs.get("focal_label", "focal"),
        section_frame.attrs.get("anchor_label", "anchor"),
        result.k,
        group_key,
        result,
    )
    if result.is_driven_by_one_cohort:
        log.warning(
            "pooled effect collapses when %s is removed; describe this as a single-%s "
            "finding, not a conserved one",
            result.loo_driver,
            group_key,
        )
    return result
