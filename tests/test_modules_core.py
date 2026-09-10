"""Tests for the module contract and the phase-2 layers (QC, integration, annotation)."""

import numpy as np
import pandas as pd
import pytest

from panspatial.annotate import (
    AMBIGUOUS, UNASSIGNED, assign_with_ambiguity, must_beat_baselines,
)
from panspatial.deconv import (
    consensus_weights, evaluate_deconvolution, mode_agreement, simulate_mixtures,
)
from panspatial.integrate import cramers_v, plan_integration
from panspatial.modules import (
    Evidence, ModuleContext, ModuleError, PlatformClass, ReplicationUnit,
    assert_can_validate, describe_registry, get_module, platform_class, registered_modules,
)
from panspatial.modules.base import ModuleResult, Provenance
from panspatial.qc import mad_thresholds, qc_mask, spatial_qc, to_microns


def make_result(module="m", evidence=Evidence.INFERRED, reliable=None, platforms=None):
    return ModuleResult(
        module=module, table=pd.DataFrame({"a": [1]}), evidence=evidence,
        replication_unit=ReplicationUnit.PATIENT,
        supported_platforms=platforms or frozenset({PlatformClass.SINGLE_CELL}),
        provenance=Provenance(module=module, module_version="1", seed=0, params={}),
        reliable_entities=frozenset(reliable) if reliable is not None else None,
    )


# ------------------------------------------------------------------ contract


def test_every_module_is_registered_and_described():
    frame = describe_registry()
    assert len(frame) >= 16
    assert frame["module"].is_unique
    assert set(frame["evidence"]) <= {"measured", "inferred", "imputed"}


@pytest.mark.parametrize("name", sorted(registered_modules()))
def test_each_module_declares_its_contract(name):
    cls = registered_modules()[name]
    assert cls.description, f"{name} has no description"
    assert cls.supported_platforms, f"{name} declares no supported platforms"
    assert cls.phase != "-", f"{name} declares no phase"


def test_platform_classification():
    assert platform_class("Xenium") is PlatformClass.SINGLE_CELL
    assert platform_class("Visium") is PlatformClass.SPOT
    assert platform_class("GeoMx") is PlatformClass.REGION
    assert platform_class(None) is PlatformClass.DISSOCIATED
    assert platform_class("SomethingNew") is PlatformClass.DISSOCIATED


def test_imputed_results_carry_their_caveat_automatically():
    r = make_result(evidence=Evidence.IMPUTED)
    assert any("not usable as validation" in n for n in r.notes)


def test_imputed_evidence_cannot_validate_anything():
    claim = make_result("spatial", Evidence.MEASURED)
    imputed = make_result("histology", Evidence.IMPUTED)
    with pytest.raises(ModuleError, match="cannot validate anything"):
        assert_can_validate(claim, imputed)


def test_weaker_evidence_cannot_validate_stronger():
    claim = make_result("nk_spatial", Evidence.MEASURED)
    inferred = make_result("deconv", Evidence.INFERRED)
    with pytest.raises(ModuleError, match="weaker than the claim"):
        assert_can_validate(claim, inferred)
    assert_can_validate(inferred, claim)  # the reverse is fine


def test_unreliable_entity_cannot_carry_a_claim():
    r = make_result(reliable=["Tumor", "CAF"])
    r.assert_entity_reliable("Tumor")
    with pytest.raises(ModuleError, match="reliability gate"):
        r.assert_entity_reliable("NK")


def test_result_refuses_an_unsupported_platform_claim():
    r = make_result(platforms=frozenset({PlatformClass.SINGLE_CELL}))
    with pytest.raises(ModuleError, match="not supported on Visium"):
        r.assert_supports("Visium", "NK distance claim")


def test_provenance_fingerprint_is_parameter_sensitive():
    a = Provenance("m", "1", 0, {"x": 1})
    b = Provenance("m", "1", 0, {"x": 1})
    c = Provenance("m", "1", 0, {"x": 2})
    assert a.fingerprint == b.fingerprint
    assert a.fingerprint != c.fingerprint


def test_unknown_module_lists_the_available_ones():
    with pytest.raises(ModuleError, match="niches"):
        get_module("nichez")


# ------------------------------------------------------------------------ QC


def test_mad_threshold_adapts_to_the_sample():
    rng = np.random.default_rng(0)
    low = mad_thresholds(rng.normal(3, 1, 2000), metric="mt", direction="upper")
    high = mad_thresholds(rng.normal(25, 4, 2000), metric="mt", direction="upper")
    assert low.upper < 10 < high.upper, "a flat 10% cut would delete the high-mito tissue"


def test_zero_mad_disables_the_filter_instead_of_flagging_everything():
    t = mad_thresholds([5.0] * 100, metric="constant")
    assert t.upper == np.inf and t.lower == -np.inf and t.mad == 0.0


def test_upper_direction_never_cuts_the_low_tail():
    t = mad_thresholds([1, 2, 3, 4, 100], metric="mt", direction="upper")
    assert t.lower == -np.inf


def test_qc_mask_removes_the_damaged_cells():
    rng = np.random.default_rng(1)
    mito = np.concatenate([rng.normal(5, 1, 200), np.full(20, 60.0)])
    feats = np.concatenate([rng.lognormal(7.5, 0.3, 200), np.full(20, 50.0)])
    counts = np.concatenate([rng.lognormal(9, 0.4, 200), np.full(20, 120.0)])
    keep, report = qc_mask(sample_id="S", percent_mito=mito, n_features=feats, n_counts=counts)
    assert not keep[200:].any(), "the injected dying cells must be removed"
    assert keep[:200].mean() > 0.9
    assert report.flags["high_mito"] >= 20
    assert "n_before" in report.to_row()


def test_qc_mask_rejects_mismatched_inputs():
    with pytest.raises(ValueError, match="same cells"):
        qc_mask(sample_id="S", percent_mito=[1, 2], n_features=[1], n_counts=[1, 2])


def test_pixel_coordinates_without_a_scale_factor_are_refused():
    with pytest.raises(ModuleError, match="Guessing here would silently rescale"):
        to_microns(np.zeros((5, 2)), platform="Visium")


def test_micron_platforms_pass_through_unchanged():
    coords = np.array([[1.0, 2.0], [3.0, 4.0]])
    np.testing.assert_array_equal(to_microns(coords, platform="Xenium"), coords)


def test_scale_factor_is_applied_when_given():
    coords = np.array([[10.0, 20.0]])
    np.testing.assert_allclose(
        to_microns(coords, platform="Visium", microns_per_unit=0.5), [[5.0, 10.0]]
    )


def test_spatial_qc_flags_coincident_and_off_tissue_points():
    coords = np.array([[0, 0], [0, 0], [100, 100], [200, 200]], dtype=float)
    keep, flags = spatial_qc(
        coords, in_tissue=[True, True, True, False], min_neighbour_distance_um=1.0
    )
    assert flags["off_tissue"] == 1
    assert flags["coincident"] == 2
    assert keep.tolist() == [False, False, True, False]


# ---------------------------------------------------------------- integration


def test_perfect_confounding_forces_within_cancer_type_integration():
    ct = np.repeat(["LUAD", "BRCA", "CRC"], 100)
    batch = np.array(["lab_" + c for c in ct])
    pat = np.array([f"{c}_P{i%6}" for i, c in enumerate(ct)])
    assert cramers_v(batch, ct) == pytest.approx(1.0, abs=1e-6)
    plan = plan_integration(batch, ct, pat)
    assert plan.scope == "within_cancer_type"
    assert plan.batch_key == "patient"
    assert "Cramer's V" in plan.rationale


def test_crossed_design_allows_global_integration():
    ct = np.repeat(["LUAD", "BRCA", "CRC"], 100)
    batch = np.array([f"lab_{i%4}" for i in range(300)])
    pat = np.array([f"P{i%20}" for i in range(300)])
    assert cramers_v(batch, ct) < 0.2
    assert plan_integration(batch, ct, pat).scope == "global"


def test_cramers_v_is_symmetric_and_bounded():
    rng = np.random.default_rng(0)
    a = rng.integers(0, 4, 500)
    b = rng.integers(0, 3, 500)
    v = cramers_v(a, b)
    assert 0.0 <= v <= 1.0
    assert v == pytest.approx(cramers_v(b, a), abs=1e-9)


def test_single_level_variable_has_no_association():
    assert cramers_v(["x"] * 50, np.arange(50)) == 0.0


def test_plan_rejects_an_unknown_method():
    with pytest.raises(ValueError, match="method must be"):
        plan_integration(["a"] * 4, ["b"] * 4, ["p"] * 4, method="magic")


# ----------------------------------------------------------------- annotation


def test_ambiguous_cells_are_labelled_not_forced():
    scores = np.array([[2.0, 0.1], [1.00, 0.98], [-1.0, -1.2]])
    assignment, margin = assign_with_ambiguity(scores, ["NK", "T"], min_margin=0.1, min_score=0.0)
    assert assignment[0] == "NK"
    assert assignment[1] == AMBIGUOUS
    assert assignment[2] == UNASSIGNED
    assert margin[0] > margin[1]


def test_exclusion_program_forces_ambiguity():
    scores = np.array([[3.0, 0.0]])
    assignment, _ = assign_with_ambiguity(
        scores, ["NK", "T"], exclusion=np.array([True])
    )
    assert assignment[0] == AMBIGUOUS


def test_assignment_needs_at_least_two_candidates():
    with pytest.raises(ValueError, match="at least two candidate labels"):
        assign_with_ambiguity(np.ones((3, 1)), ["NK"])


def test_foundation_model_must_beat_every_baseline():
    baselines = {"hvg": 0.64, "pca": 0.60, "harmony": 0.66, "scvi": 0.68}
    weak = must_beat_baselines("scGPT", 0.61, baselines, n_held_out_patients=10)
    strong = must_beat_baselines("Nicheformer", 0.74, baselines, n_held_out_patients=10)
    assert not weak.admitted and "does not beat" in weak.reason
    assert strong.admitted and strong.best_baseline == "scvi"


def test_beating_only_some_baselines_is_a_rejection():
    """Beating HVG but not scVI is the documented zero-shot failure mode."""
    v = must_beat_baselines(
        "scGPT", 0.65, {"hvg": 0.64, "pca": 0.60, "harmony": 0.66, "scvi": 0.68},
        n_held_out_patients=10,
    )
    assert not v.admitted


def test_too_few_held_out_patients_is_a_rejection():
    v = must_beat_baselines(
        "Nicheformer", 0.95, {"hvg": 0.1, "pca": 0.1, "harmony": 0.1, "scvi": 0.1},
        n_held_out_patients=2,
    )
    assert not v.admitted and "held-out patient" in v.reason


def test_missing_required_baseline_blocks_admission():
    v = must_beat_baselines("scGPT", 0.99, {"hvg": 0.5}, n_held_out_patients=10)
    assert not v.admitted and "required baseline" in v.reason


# -------------------------------------------------------------- deconvolution


def test_simulated_mixtures_have_exact_known_proportions():
    rng = np.random.default_rng(0)
    counts = rng.poisson(3, size=(50, 300))
    labels = np.repeat(["Tumor", "CAF", "NK"], 100)
    sim, truth, types = simulate_mixtures(
        counts, labels, n_spots=100, cells_per_spot=8, rng=rng
    )
    assert sim.shape == (50, 100)
    assert truth.shape == (100, 3)
    np.testing.assert_allclose(truth.sum(axis=1), 1.0)
    assert types == ["CAF", "NK", "Tumor"]


def test_rare_boost_lifts_rare_type_representation():
    rng = np.random.default_rng(0)
    counts = rng.poisson(3, size=(30, 400))
    labels = np.array(["Tumor"] * 380 + ["NK"] * 20)   # NK at 5%
    _, natural, types = simulate_mixtures(counts, labels, n_spots=200, rng=np.random.default_rng(1))
    _, boosted, _ = simulate_mixtures(
        counts, labels, n_spots=200, rare_boost=1.0, rng=np.random.default_rng(1)
    )
    nk = types.index("NK")
    assert boosted[:, nk].mean() > 5 * natural[:, nk].mean()


def test_reliability_gate_passes_accurate_and_fails_noisy_types():
    rng = np.random.default_rng(0)
    truth = rng.dirichlet([4, 4, 0.6], size=300)
    est = truth.copy()
    est[:, 2] = rng.dirichlet([1] * 300) * 3        # NK column is pure noise
    report = evaluate_deconvolution(truth, est, ["Tumor", "CAF", "NK"], min_r=0.7)
    assert "Tumor" in report.reliable and "CAF" in report.reliable
    assert "NK" in report.unreliable
    assert "may not carry a conclusion" in [r.reason for r in report.rows
                                            if r.cell_type == "NK"][0]


def test_types_too_rare_to_evaluate_fail_with_that_reason():
    truth = np.zeros((100, 2))
    truth[:, 0] = 1.0
    truth[:3, 1] = 0.5
    report = evaluate_deconvolution(truth, truth.copy(), ["Tumor", "NK"], min_spots_present=20)
    nk = [r for r in report.rows if r.cell_type == "NK"][0]
    assert not nk.reliable and "too rare" in nk.reason


def test_consensus_flags_types_the_methods_disagree_on():
    rng = np.random.default_rng(0)
    spots = [f"s{i}" for i in range(80)]
    types = ["Tumor", "CAF", "NK"]
    shared = pd.DataFrame(rng.dirichlet([5, 3, 0.5], 80), index=spots, columns=types)
    methods = {}
    for m in ("cell2location", "rctd", "card"):
        f = shared.copy()
        f["NK"] = rng.random(80)          # each method invents its own NK signal
        methods[m] = f
    result = consensus_weights(methods)
    assert "Tumor" in result.concordant_types and "CAF" in result.concordant_types
    assert "NK" not in result.concordant_types


def test_consensus_of_one_method_is_refused():
    frame = pd.DataFrame({"Tumor": [0.5, 0.5]})
    with pytest.raises(ValueError, match="at least two methods"):
        consensus_weights({"rctd": frame})


def test_mode_agreement_flags_similarity_carried_weight():
    weights = pd.DataFrame({"Tumor": [0.6] * 50, "NK": [0.1] * 50})
    calls = {f"s{i}": ["Tumor"] for i in range(50)}   # doublet mode never calls NK
    out = mode_agreement(weights, calls)
    nk = out[out["cell_type"] == "NK"].iloc[0]
    assert nk["similarity_carried"]
    assert not out[out["cell_type"] == "Tumor"].iloc[0]["similarity_carried"]
