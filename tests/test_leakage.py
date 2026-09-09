"""Tests for patient-level CV splitting and the leakage assertions."""

import numpy as np
import pytest

from panspatial.stats.leakage import (
    LeakageError,
    assert_no_leakage,
    patient_blocked_splits,
    spatial_buffer_filter,
)


def cohort(n_patients=10, per_patient=50, cancer_types=("LUAD", "BRCA")):
    patients = np.repeat([f"P{i:02d}" for i in range(n_patients)], per_patient)
    types = np.repeat(
        [cancer_types[i % len(cancer_types)] for i in range(n_patients)], per_patient
    )
    return patients, types


def test_no_patient_appears_in_both_folds():
    patients, _ = cohort()
    for train, test in patient_blocked_splits(patients, n_splits=5, rng=np.random.default_rng(0)):
        assert not set(patients[train]) & set(patients[test])


def test_every_observation_is_tested_exactly_once():
    patients, _ = cohort()
    splits = patient_blocked_splits(patients, n_splits=5, rng=np.random.default_rng(0))
    tested = np.concatenate([test for _, test in splits])
    np.testing.assert_array_equal(np.sort(tested), np.arange(len(patients)))


def test_stratification_spreads_cancer_types_across_folds():
    patients, types = cohort(n_patients=12, cancer_types=("LUAD", "BRCA", "CRC"))
    splits = patient_blocked_splits(
        patients, n_splits=3, rng=np.random.default_rng(1), stratify_by=types
    )
    for _, test in splits:
        assert len(set(types[test])) == 3


def test_stratify_rejects_label_varying_within_patient():
    patients = np.array(["P1", "P1", "P2", "P2"])
    bad = np.array(["LUAD", "BRCA", "CRC", "CRC"])
    with pytest.raises(ValueError, match="multiple strata"):
        patient_blocked_splits(patients, n_splits=2, stratify_by=bad)


def test_single_patient_cohort_is_refused():
    with pytest.raises(LeakageError, match="one patient"):
        patient_blocked_splits(np.array(["P1"] * 100))


def test_splits_are_reproducible_under_a_seed():
    patients, _ = cohort()
    a = patient_blocked_splits(patients, n_splits=4, rng=np.random.default_rng(42))
    b = patient_blocked_splits(patients, n_splits=4, rng=np.random.default_rng(42))
    for (tr_a, te_a), (tr_b, te_b) in zip(a, b):
        np.testing.assert_array_equal(tr_a, tr_b)
        np.testing.assert_array_equal(te_a, te_b)


# ------------------------------------------------------------------ assertions


def test_spot_level_split_is_caught():
    """The classic mistake: shuffling spots into folds ignoring the patient."""
    patients, _ = cohort(n_patients=4, per_patient=25)
    rng = np.random.default_rng(0)
    idx = rng.permutation(len(patients))
    train, test = idx[:80], idx[80:]
    with pytest.raises(LeakageError, match="both train and test"):
        assert_no_leakage(train, test, patients)


def test_overlapping_indices_are_caught_first():
    patients, _ = cohort(n_patients=4, per_patient=25)
    with pytest.raises(LeakageError, match="appear in both"):
        assert_no_leakage(np.arange(50), np.arange(40, 60), patients)


def test_clean_patient_split_passes():
    patients, _ = cohort()
    for train, test in patient_blocked_splits(patients, n_splits=5, rng=np.random.default_rng(0)):
        assert_no_leakage(train, test, patients)


def test_adjacent_spots_within_a_section_are_caught():
    """Region-held-out design inside one section: neighbours across the boundary leak."""
    coords = np.stack([np.arange(200, dtype=float) * 50.0, np.zeros(200)], axis=1)
    sections = np.array(["S1"] * 200)
    patients = np.array(["P1"] * 200)
    train, test = np.arange(0, 100), np.arange(100, 200)
    with pytest.raises(LeakageError, match="within 200.0 um"):
        assert_no_leakage(
            train,
            test,
            patients,
            coords=coords,
            section_ids=sections,
            buffer_um=200.0,
            allow_shared_patients=True,
        )


def test_region_holdout_without_a_buffer_is_refused():
    """Relaxing the patient check with no buffer leaves nothing guarding the split."""
    patients = np.array(["P1"] * 100)
    with pytest.raises(ValueError, match="disables every leakage guard"):
        assert_no_leakage(
            np.arange(50), np.arange(50, 100), patients, allow_shared_patients=True
        )


def test_buffer_filter_removes_only_the_boundary_strip():
    coords = np.stack([np.arange(200, dtype=float) * 50.0, np.zeros(200)], axis=1)
    sections = np.array(["S1"] * 200)
    train, test = np.arange(0, 100), np.arange(100, 200)
    kept = spatial_buffer_filter(train, test, coords, sections, buffer_um=200.0)
    # test starts at x=5000 um; training spots at 4850/4900/4950 are strictly within 200 um
    assert len(kept) == 97
    assert kept.max() == 96


def test_buffer_ignores_distances_across_different_sections():
    """Coordinates from two sections share no frame of reference."""
    coords = np.zeros((100, 2))
    sections = np.array(["S1"] * 50 + ["S2"] * 50)
    train, test = np.arange(0, 50), np.arange(50, 100)
    kept = spatial_buffer_filter(train, test, coords, sections, buffer_um=1000.0)
    np.testing.assert_array_equal(kept, train)


def test_buffer_of_zero_is_a_no_op():
    coords = np.zeros((10, 2))
    sections = np.array(["S1"] * 10)
    train = np.arange(5)
    np.testing.assert_array_equal(
        spatial_buffer_filter(train, np.arange(5, 10), coords, sections, buffer_um=0.0), train
    )
