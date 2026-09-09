"""Tests for the statistical layer, including the property that motivates it:
a torus-shift null is calibrated on autocorrelated patterns where a label shuffle is not.
"""

import numpy as np
import pytest

from panspatial.stats.spatial_stats import (
    benjamini_hochberg,
    fraction_within,
    mean_nn_distance,
    random_effects_meta,
    torus_shift_permutation_test,
)


def clustered_pattern(rng, n_clusters=6, per_cluster=40, spread=30.0, extent=1000.0):
    """A clustered point pattern, like immune cells in tissue."""
    centers = rng.uniform(0, extent, size=(n_clusters, 2))
    pts = np.repeat(centers, per_cluster, axis=0) + rng.normal(0, spread, size=(n_clusters * per_cluster, 2))
    return np.clip(pts, 0, extent)


# ----------------------------------------------------------------- statistics


def test_mean_nn_distance_simple_geometry():
    focal = np.array([[0.0, 0.0], [10.0, 0.0]])
    anchor = np.array([[0.0, 3.0], [10.0, 4.0]])
    assert mean_nn_distance(focal, anchor) == pytest.approx(3.5)


def test_mean_nn_distance_rejects_bad_shape():
    with pytest.raises(ValueError, match=r"\(n, 2\)"):
        mean_nn_distance(np.zeros((5, 3)), np.zeros((5, 2)))


def test_fraction_within_counts_correctly():
    stat = fraction_within(5.0)
    focal = np.array([[0.0, 0.0], [100.0, 0.0], [0.0, 4.0]])
    anchor = np.array([[0.0, 1.0]])
    assert stat(focal, anchor) == pytest.approx(2 / 3)


# ----------------------------------------------------------- permutation test


def test_torus_shift_detects_real_co_registration():
    """Focal points placed on top of anchors must come out significantly attracted."""
    rng = np.random.default_rng(0)
    anchor = clustered_pattern(rng)
    focal = anchor[rng.choice(len(anchor), 150, replace=False)] + rng.normal(0, 5, size=(150, 2))
    res = torus_shift_permutation_test(focal, anchor, n_perms=499, rng=np.random.default_rng(1))
    assert res.observed < res.null_mean
    assert res.p_value < 0.01
    assert res.n_focal == 150


def test_torus_shift_p_value_never_zero():
    rng = np.random.default_rng(3)
    anchor = clustered_pattern(rng)
    focal = anchor[:100] + rng.normal(0, 1, size=(100, 2))
    res = torus_shift_permutation_test(focal, anchor, n_perms=99, rng=np.random.default_rng(4))
    assert res.p_value >= 1 / 100


def test_torus_shift_is_calibrated_where_label_shuffle_is_not():
    """The reason this module exists.

    Under a true null -- two independently clustered patterns with no relationship -- a
    calibrated test rejects at roughly its nominal rate. The torus shift preserves each
    pattern's own clustering and randomizes only their registration, so it stays calibrated.
    A label shuffle over the pooled positions destroys that clustering: the null focal set
    is drawn from the union, so null distances collapse and the test reports spurious
    "exclusion" almost every time.
    """
    n_reps, alpha = 24, 0.05
    torus_rejects = 0
    shuffle_rejects = 0

    for rep in range(n_reps):
        rng = np.random.default_rng(100 + rep)
        anchor = clustered_pattern(rng, n_clusters=5, per_cluster=40)
        focal = clustered_pattern(rng, n_clusters=5, per_cluster=30)  # independent of anchor

        torus = torus_shift_permutation_test(
            focal,
            anchor,
            alternative="two-sided",
            n_perms=199,
            rng=np.random.default_rng(900 + rep),
        )
        torus_rejects += torus.p_value < alpha

        # The naive alternative: shuffle focal/anchor labels over the pooled positions.
        prng = np.random.default_rng(900 + rep)
        pooled = np.vstack([focal, anchor])
        n_focal = len(focal)
        observed = mean_nn_distance(focal, anchor)
        null = np.empty(199)
        for i in range(199):
            perm = prng.permutation(len(pooled))
            null[i] = mean_nn_distance(pooled[perm[:n_focal]], pooled[perm[n_focal:]])
        centered = np.abs(null - null.mean())
        naive_p = (1 + np.sum(centered >= abs(observed - null.mean()))) / 200
        shuffle_rejects += naive_p < alpha

    assert torus_rejects <= 0.25 * n_reps, f"torus-shift null over-rejects: {torus_rejects}/{n_reps}"
    assert shuffle_rejects >= 0.75 * n_reps, (
        f"expected the naive label shuffle to over-reject under a true null, got "
        f"{shuffle_rejects}/{n_reps} (torus: {torus_rejects}/{n_reps})"
    )


def test_torus_shift_respects_tissue_mask():
    """With a mask, shifted points outside the tissue are discarded rather than counted."""
    rng = np.random.default_rng(7)
    anchor = clustered_pattern(rng, n_clusters=4, per_cluster=50, extent=500.0)
    focal = anchor[:80] + rng.normal(0, 4, size=(80, 2))
    inside = lambda xy: (xy[:, 0] >= 0) & (xy[:, 0] <= 500) & (xy[:, 1] >= 0) & (xy[:, 1] <= 500)
    res = torus_shift_permutation_test(
        focal, anchor, n_perms=99, inside_tissue=inside, rng=np.random.default_rng(8)
    )
    assert np.isfinite(res.null_mean)


def test_torus_shift_refuses_tiny_focal_population():
    with pytest.raises(ValueError, match="noise"):
        torus_shift_permutation_test(
            np.zeros((2, 2)), np.ones((10, 2)), rng=np.random.default_rng(0)
        )


def test_torus_shift_rejects_bad_alternative():
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError, match="alternative"):
        torus_shift_permutation_test(
            clustered_pattern(rng), clustered_pattern(rng), alternative="smaller", rng=rng
        )


# ---------------------------------------------------------------- meta-analysis


def test_random_effects_recovers_homogeneous_effect():
    effects = [0.50, 0.48, 0.52, 0.49, 0.51]
    variances = [0.01] * 5
    res = random_effects_meta(effects, variances)
    assert res.estimate == pytest.approx(0.50, abs=0.02)
    assert res.tau2 == pytest.approx(0.0, abs=1e-6)
    assert res.i2 == pytest.approx(0.0, abs=1e-6)
    assert res.k == 5
    assert res.p_value < 0.001


def test_random_effects_widens_ci_under_heterogeneity():
    variances = [0.01] * 5
    tight = random_effects_meta([0.50, 0.48, 0.52, 0.49, 0.51], variances)
    spread = random_effects_meta([0.10, 0.90, 0.20, 0.80, 0.50], variances)
    assert spread.tau2 > tight.tau2
    assert spread.i2 > 0.8
    assert (spread.ci_high - spread.ci_low) > (tight.ci_high - tight.ci_low)


def test_random_effects_flags_single_cohort_driver():
    """One cohort carrying the whole effect must be visible, not averaged away."""
    res = random_effects_meta(
        effects=[2.0, 0.02, -0.01, 0.03, 0.0],
        variances=[0.02, 0.02, 0.02, 0.02, 0.02],
        labels=["LUAD", "BRCA", "CRC", "GBM", "PDAC"],
    )
    assert res.loo_driver == "LUAD"
    assert res.is_driven_by_one_cohort


def test_random_effects_requires_two_units():
    with pytest.raises(ValueError, match="at least 2 units"):
        random_effects_meta([0.5], [0.01])


def test_random_effects_drops_invalid_variances():
    res = random_effects_meta([0.5, 0.4, 0.6], [0.01, 0.0, 0.01])
    assert res.k == 2


# ------------------------------------------------------------- multiple testing


def test_bh_matches_hand_computed_values():
    p = np.array([0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205])
    expected = np.array([0.008, 0.032, 0.0672, 0.0672, 0.0672, 0.08, 0.074 * 8 / 7, 0.205])
    np.testing.assert_allclose(benjamini_hochberg(p), expected, rtol=1e-6)


def test_bh_is_monotone_and_bounded():
    rng = np.random.default_rng(0)
    p = rng.uniform(size=200)
    adj = benjamini_hochberg(p)
    assert np.all(adj >= p - 1e-12)
    assert np.all(adj <= 1.0)
    order = np.argsort(p)
    assert np.all(np.diff(adj[order]) >= -1e-12)


def test_bh_family_size_makes_correction_stricter():
    p = np.array([0.001, 0.01, 0.02])
    partial = benjamini_hochberg(p)
    full_family = benjamini_hochberg(p, family_size=1000)
    assert np.all(full_family >= partial)
    assert full_family[0] == pytest.approx(0.001 * 1000 / 1)


def test_bh_rejects_undersized_family():
    with pytest.raises(ValueError, match="family_size"):
        benjamini_hochberg([0.1, 0.2, 0.3], family_size=2)


def test_bh_rejects_out_of_range():
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        benjamini_hochberg([0.1, 1.5])
