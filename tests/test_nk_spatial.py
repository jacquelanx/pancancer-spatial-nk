"""Tests for the NK spatial layer, using simulated tissue (no scanpy required)."""

import numpy as np
import pytest

from panspatial.nk.nk_spatial import (
    NK_NEGATIVE_MARKERS,
    NK_POSITIVE_MARKERS,
    NK_STATE_PROGRAMS,
    DeconvolutionWarning,
    cohort_meta_analysis,
    distance_to_anchor,
    per_section_niche_test,
    require_single_cell_resolution,
)


def simulate_section(rng, *, n_nk=120, n_caf=300, attraction=0.0, extent=2000.0):
    """One section: CAFs in stromal clusters, NK cells drawn towards them by `attraction`.

    attraction=0 -> NK cells independently clustered elsewhere (true null).
    attraction=1 -> NK cells sit on CAF positions.
    """
    caf_centers = rng.uniform(0, extent, size=(5, 2))
    caf = np.repeat(caf_centers, n_caf // 5, axis=0) + rng.normal(0, 40, size=(n_caf // 5 * 5, 2))

    nk_centers = rng.uniform(0, extent, size=(4, 2))
    nk = np.repeat(nk_centers, n_nk // 4, axis=0) + rng.normal(0, 40, size=(n_nk // 4 * 4, 2))
    if attraction > 0:
        targets = caf[rng.choice(len(caf), len(nk))]
        nk = (1 - attraction) * nk + attraction * targets + rng.normal(0, 15, size=nk.shape)

    coords = np.vstack([np.clip(nk, 0, extent), np.clip(caf, 0, extent)])
    labels = np.array(["NK"] * len(nk) + ["mCAF"] * len(caf))
    return coords, labels


def simulate_cohort(n_patients=8, sections_per_patient=1, attraction=0.0, seed=0):
    coords, labels, sections, patients = [], [], [], []
    for p in range(n_patients):
        for s in range(sections_per_patient):
            rng = np.random.default_rng(seed * 1000 + p * 10 + s)
            c, lab = simulate_section(rng, attraction=attraction)
            coords.append(c)
            labels.append(lab)
            sections.extend([f"P{p:02d}_S{s}"] * len(c))
            patients.extend([f"P{p:02d}"] * len(c))
    return np.vstack(coords), np.concatenate(labels), np.array(sections), np.array(patients)


# ----------------------------------------------------------------- marker panels


def test_discriminating_markers_are_not_shared_with_excluded_populations():
    excluded = {g for genes in NK_NEGATIVE_MARKERS.values() for g in genes}
    assert not set(NK_POSITIVE_MARKERS["discriminating"]) & excluded


def test_state_programs_are_multi_gene_and_distinct():
    for state, genes in NK_STATE_PROGRAMS.items():
        assert len(genes) >= 4, f"{state} is defined by too few genes to be a program"
    cytotoxic = set(NK_STATE_PROGRAMS["cytotoxic"])
    dysfunctional = set(NK_STATE_PROGRAMS["dysfunctional"])
    assert not cytotoxic & dysfunctional


# --------------------------------------------------------------- platform gating


@pytest.mark.parametrize("platform", ["Visium", "Visium HD", "Slide-seq", "GeoMx"])
def test_spot_platforms_are_blocked_for_nk_claims(platform):
    with pytest.raises(DeconvolutionWarning, match="deconvolution weight"):
        require_single_cell_resolution(platform)


@pytest.mark.parametrize("platform", ["Xenium", "CosMx", "MERFISH"])
def test_single_cell_platforms_are_allowed(platform):
    require_single_cell_resolution(platform)


# ------------------------------------------------------------------- distances


def test_distances_are_computed_within_section_only():
    coords = np.array([[0.0, 0.0], [1000.0, 0.0], [10.0, 0.0], [1010.0, 0.0]])
    labels = np.array(["NK", "NK", "mCAF", "mCAF"])
    sections = np.array(["A", "B", "A", "B"])
    out = distance_to_anchor(coords, labels, sections, focal_label="NK", anchor_label="mCAF")
    assert set(out) == {"A", "B"}
    assert out["A"][0] == pytest.approx(10.0)
    assert out["B"][0] == pytest.approx(10.0)


def test_sections_missing_a_population_are_skipped_not_crashed():
    coords = np.array([[0.0, 0.0], [5.0, 0.0], [100.0, 0.0]])
    labels = np.array(["NK", "NK", "mCAF"])
    sections = np.array(["A", "A", "B"])
    out = distance_to_anchor(coords, labels, sections, focal_label="NK", anchor_label="mCAF")
    assert out == {}


def test_distance_rejects_mismatched_lengths():
    with pytest.raises(ValueError, match="same observations"):
        distance_to_anchor(np.zeros((5, 2)), ["NK"] * 4, ["A"] * 5, focal_label="NK", anchor_label="X")


# ---------------------------------------------------------------- niche testing


def test_per_section_test_recovers_real_attraction():
    coords, labels, sections, patients = simulate_cohort(n_patients=4, attraction=0.8, seed=1)
    frame = per_section_niche_test(
        coords, labels, sections, patients,
        focal_label="NK", anchor_label="mCAF", n_perms=199, seed=5,
    )
    assert len(frame) == 4
    assert (frame["log_ratio"] < 0).all(), "attracted NK cells must be closer than the null"
    # A single section is a small, noisy experiment: the effect direction is consistent but
    # per-section significance is not guaranteed. Power comes from pooling, which is the
    # point of running the test per section and meta-analyzing rather than pooling cells.
    assert (frame["p_value"] < 0.05).sum() >= 3
    assert cohort_meta_analysis(frame).p_value < 0.01


def test_per_section_test_is_quiet_under_the_null():
    coords, labels, sections, patients = simulate_cohort(n_patients=8, attraction=0.0, seed=2)
    frame = per_section_niche_test(
        coords, labels, sections, patients,
        focal_label="NK", anchor_label="mCAF", n_perms=199, seed=7,
    )
    assert (frame["p_value"] < 0.05).sum() <= 2, "too many rejections under a true null"


def test_sparse_sections_are_excluded_and_recorded():
    coords, labels, sections, patients = simulate_cohort(n_patients=3, attraction=0.5, seed=3)
    frame = per_section_niche_test(
        coords, labels, sections, patients,
        focal_label="NK", anchor_label="mCAF", n_perms=99, min_focal=1000, seed=9,
    )
    assert len(frame) == 0
    assert len(frame.attrs["skipped_sections"]) == 3
    assert frame.attrs["skipped_sections"][0]["n_focal"] == 120


# --------------------------------------------------------------- meta-analysis


def test_cohort_meta_analysis_detects_conserved_attraction():
    coords, labels, sections, patients = simulate_cohort(n_patients=8, attraction=0.7, seed=4)
    frame = per_section_niche_test(
        coords, labels, sections, patients,
        focal_label="NK", anchor_label="mCAF", n_perms=199, seed=11,
    )
    result = cohort_meta_analysis(frame)
    assert result.k == 8
    assert result.estimate < 0
    assert result.p_value < 0.01
    assert not result.is_driven_by_one_cohort


def test_multiple_sections_from_one_patient_count_once():
    """Three sections from one patient must not become three units of evidence."""
    coords, labels, sections, patients = simulate_cohort(
        n_patients=4, sections_per_patient=3, attraction=0.7, seed=6
    )
    frame = per_section_niche_test(
        coords, labels, sections, patients,
        focal_label="NK", anchor_label="mCAF", n_perms=99, seed=13,
    )
    assert len(frame) == 12
    result = cohort_meta_analysis(frame)
    assert result.k == 4, "pooling must be over patients, not sections"


def test_meta_analysis_refuses_a_single_section():
    coords, labels, sections, patients = simulate_cohort(n_patients=1, attraction=0.7, seed=8)
    frame = per_section_niche_test(
        coords, labels, sections, patients,
        focal_label="NK", anchor_label="mCAF", n_perms=99, seed=15,
    )
    with pytest.raises(ValueError, match="between-patient variance"):
        cohort_meta_analysis(frame)
