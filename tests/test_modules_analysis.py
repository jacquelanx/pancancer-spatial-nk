"""Tests for the phase 4-6 analysis layers."""

import numpy as np
import pandas as pd
import pytest

from panspatial.alignment import check_alignment, select_aligner
from panspatial.comodules import module_preservation
from panspatial.domains import BOUNDARY_UNCERTAIN, consensus_domains, select_method
from panspatial.dynamics import corroborate_direction, view_family
from panspatial.grn import conserved_regulons, regulon_niche_enrichment
from panspatial.histology import per_organ_performance
from panspatial.interactions import (
    contact_frequency, robust_rank_aggregate, spatial_constraint,
)
from panspatial.malignant import (
    malignant_fraction_table, match_meta_programs, meta_program_table, select_cnv_caller,
)
from panspatial.modules import ModuleError
from panspatial.niches import composition_matrix, match_niches, recurrence_test
from panspatial.predict import patient_level_cv
from panspatial.stats.leakage import LeakageError
from panspatial.translate import TargetCandidate, prioritise_targets, score_signature


# ------------------------------------------------------------------- domains


@pytest.mark.parametrize("platform,expected", [
    ("Visium", "graphst"), ("ST (legacy)", "bayesspace"), ("Xenium", "banksy"),
    ("CosMx", "banksy"), ("Slide-seq", "bayesspace"),
])
def test_domain_method_follows_the_platform_benchmark(platform, expected):
    assert select_method(platform).primary == expected


def test_domain_selection_requires_a_platform():
    with pytest.raises(ModuleError, match="needs a platform"):
        select_method(None)


def test_consensus_is_invariant_to_arbitrary_cluster_names():
    """Cluster indices are arbitrary; comparing them directly would measure nothing."""
    a = np.array(["A"] * 40 + ["B"] * 40)
    b = np.array([9] * 40 + [3] * 40)
    c = np.array(["red"] * 40 + ["blue"] * 40)
    result = consensus_domains({"graphst": a, "stagate": b, "banksy": c})
    assert result.mean_agreement == pytest.approx(1.0)
    assert result.uncertain.sum() == 0


def test_spots_where_methods_split_evenly_are_boundary_uncertain():
    a = np.array(["A"] * 50 + ["B"] * 50)
    b = np.array(["A"] * 50 + ["B"] * 50)
    c = np.array(["A"] * 25 + ["B"] * 75)     # disagrees on spots 25-49
    result = consensus_domains({"m1": a, "m2": b, "m3": c}, min_agreement=0.7)
    assert result.uncertain[25:50].all()
    assert (result.labels[25:50] == BOUNDARY_UNCERTAIN).all()
    assert not result.uncertain[:25].any()


def test_consensus_of_one_method_is_refused():
    with pytest.raises(ValueError, match="at least two labelings"):
        consensus_domains({"only": np.array(["A", "B"])})


# -------------------------------------------------------------------- niches


def test_composition_matrix_rows_sum_to_one():
    labels = np.array(["n1"] * 10 + ["n2"] * 10)
    types = np.array(["Tumor"] * 6 + ["CAF"] * 4 + ["NK"] * 10)
    mat, niches, tnames = composition_matrix(labels, types)
    np.testing.assert_allclose(mat.sum(axis=1), 1.0)
    assert niches == ["n1", "n2"] and set(tnames) == {"CAF", "NK", "Tumor"}


def test_dissimilar_local_niches_stay_unmatched():
    reference = np.array([[0.9, 0.05, 0.05], [0.05, 0.9, 0.05]])
    sample = {"S1": (np.array([[0.33, 0.33, 0.34]]), ["odd"])}
    matches = match_niches(reference, ["core", "stroma"], sample, min_similarity=0.9)
    assert not matches[0].matched
    assert matches[0].reference_niche is None


def test_recurrence_requires_both_patients_and_cancer_types():
    reference = np.eye(3)
    samples, patient_of, cancer_of = {}, {}, {}
    for i in range(12):
        # niche 0 everywhere; niche 2 only in LUAD samples
        comp = reference[[0, 1]] if i % 4 else reference[[0, 1, 2]]
        names = [f"n{j}" for j in range(len(comp))]
        s = f"S{i}"
        samples[s] = (comp, names)
        patient_of[s] = f"P{i}"
        cancer_of[s] = "LUAD" if i % 4 == 0 else ["BRCA", "CRC", "GBM"][i % 3]
    matches = match_niches(reference, ["a", "b", "c"], samples, min_similarity=0.8)
    frame = recurrence_test(matches, patient_of, cancer_of, min_cancer_types=3)
    a = frame[frame["niche"] == "a"].iloc[0]
    c = frame[frame["niche"] == "c"].iloc[0]
    assert a["conserved"]
    assert not c["conserved"] and "restricted" in c["verdict"]


def test_recurrence_counts_patients_not_sections():
    """Three sections from one patient are one observation."""
    reference = np.eye(2)
    samples = {f"S{i}": (reference, ["n0", "n1"]) for i in range(6)}
    patient_of = {f"S{i}": "P0" for i in range(6)}       # all one patient
    cancer_of = {f"S{i}": "LUAD" for i in range(6)}
    frame = recurrence_test(match_niches(reference, ["a", "b"], samples), patient_of, cancer_of)
    assert (frame["n_patients_total"] == 1).all()
    assert not frame["conserved"].any()


def test_null_prevalence_must_be_a_probability():
    reference = np.eye(2)
    samples = {"S0": (reference, ["n0", "n1"])}
    with pytest.raises(ValueError, match="null_prevalence"):
        recurrence_test(
            match_niches(reference, ["a", "b"], samples), {"S0": "P0"}, {"S0": "LUAD"},
            null_prevalence=1.5,
        )


# ----------------------------------------------------------------- malignant


@pytest.mark.parametrize("allele,depth,expected", [
    (True, 10000, "numbat"),
    (True, 1200, "copykat"),      # Numbat degrades sharply at low depth
    (False, 8000, "copykat"),     # best expression-only tool
    (False, 800, "copykat"),
])
def test_cnv_caller_choice_follows_allele_and_depth(allele, depth, expected):
    assert select_cnv_caller(has_allele_info=allele, median_umi_per_cell=depth).primary == expected


def test_very_low_depth_attaches_a_caveat():
    choice = select_cnv_caller(has_allele_info=False, median_umi_per_cell=600)
    assert choice.depth_caveat and "very low depth" in choice.depth_caveat


def test_malignant_fraction_reports_depth_comparability():
    calls = np.array(["malignant"] * 30 + ["normal"] * 70)
    samples = np.array(["S1"] * 50 + ["S2"] * 50)
    frame = malignant_fraction_table(calls, samples, {"S1": 8000, "S2": 500}, "copykat")
    assert frame.set_index("sample_id").loc["S1", "depth_comparable"]
    assert not frame.set_index("sample_id").loc["S2", "depth_comparable"]


def test_meta_programs_need_several_patients():
    programs = {"P1": {"p": ["A", "B", "C", "D"]}, "P2": {"p": ["A", "B", "C", "E"]}}
    assert match_meta_programs(programs, min_patients=3) == []
    assert len(match_meta_programs(programs, min_patients=2)) == 1


def test_programs_within_one_patient_are_never_merged():
    programs = {"P1": {"a": ["X", "Y", "Z"], "b": ["X", "Y", "Z"]}}
    assert match_meta_programs(programs, min_patients=1) != []
    metas = match_meta_programs(programs, min_patients=1)
    assert all(m.n_patients == 1 for m in metas)
    assert len(metas) == 2, "identical programs from the same patient must stay separate"


def test_single_cancer_type_program_is_not_pan_cancer():
    programs = {f"P{i}": {"p": ["A", "B", "C", "D", "E"]} for i in range(4)}
    ct = {f"P{i}": "LUAD" for i in range(4)}
    metas = match_meta_programs(programs, cancer_type_of=ct, min_patients=3)
    frame = meta_program_table(metas, min_cancer_types=3)
    assert not frame["pan_cancer"].any()
    assert "not a pan-cancer program" in frame.iloc[0]["verdict"]


# -------------------------------------------------------------- interactions


def test_separated_populations_have_no_contacts():
    coords = np.vstack([np.random.default_rng(0).normal(0, 10, (40, 2)),
                        np.random.default_rng(1).normal(500, 10, (40, 2))])
    types = np.array(["NK"] * 40 + ["CAF"] * 40)
    frame = contact_frequency(coords, types, max_distance_um=30)
    cross = frame[(frame["sender"] == "NK") & (frame["receiver"] == "CAF")].iloc[0]
    assert cross["n_contacts"] == 0 and not cross["in_contact"]


def test_spatial_constraint_flags_expression_only_calls():
    interactions = pd.DataFrame({"sender": ["NK", "NK"], "receiver": ["CAF", "NK"],
                                 "score": [0.9, 0.5]})
    contacts = pd.DataFrame({"sender": ["NK", "NK"], "receiver": ["CAF", "NK"],
                             "n_contacts": [0, 120]})
    out = spatial_constraint(interactions, contacts)
    assert not out.iloc[0]["spatially_supported"], "high score with no contact must be flagged"
    assert out.iloc[1]["spatially_supported"]


def test_rank_aggregation_prefers_agreement_over_a_single_top_hit():
    ranks = {
        "cellphonedb": {"agreed": 2, "solo": 1, "filler1": 3, "filler2": 4, "filler3": 5},
        "cellchat":    {"agreed": 1, "solo": 5, "filler1": 2, "filler2": 3, "filler3": 4},
        "nichenet":    {"agreed": 2, "solo": 5, "filler1": 1, "filler2": 3, "filler3": 4},
    }
    frame = robust_rank_aggregate(ranks).set_index("interaction")
    assert frame.loc["agreed", "consensus_rank"] < frame.loc["solo", "consensus_rank"]


def test_rank_aggregation_needs_two_methods():
    with pytest.raises(ValueError, match="at least two methods"):
        robust_rank_aggregate({"only": {"a": 1}})


# ----------------------------------------------------------------------- GRN


def _regulon_fixture(seed=0, effect=1.2):
    rng = np.random.default_rng(seed)
    n_patients, per_patient = 8, 200
    patient = np.repeat([f"P{i}" for i in range(n_patients)], per_patient)
    niche = np.tile(np.array(["target"] * 100 + ["other"] * 100), n_patients)
    real = rng.normal(0, 1, len(patient)) + (niche == "target") * effect
    noise = rng.normal(0, 1, len(patient))
    return pd.DataFrame({"TF_real": real, "TF_noise": noise}), niche, patient


def test_regulon_enrichment_separates_signal_from_noise():
    activity, niche, patient = _regulon_fixture()
    frame = regulon_niche_enrichment(activity, niche, patient, target_niche="target")
    real = frame[frame["regulon"] == "TF_real"].iloc[0]
    noise = frame[frame["regulon"] == "TF_noise"].iloc[0]
    assert real["hedges_g"] > 0.8 and real["fdr"] < 0.01
    assert abs(noise["hedges_g"]) < 0.3 and noise["fdr"] > 0.05
    assert real["n_patients"] == 8, "the n must be patients, not cells"


def test_conservation_requires_a_consistent_direction():
    up = pd.DataFrame({"regulon": ["TF1"], "hedges_g": [1.0], "p_value": [0.001], "fdr": [0.001]})
    down = pd.DataFrame({"regulon": ["TF1"], "hedges_g": [-1.0], "p_value": [0.001], "fdr": [0.001]})
    out = conserved_regulons({"LUAD": up, "BRCA": up, "CRC": down})
    assert not out.iloc[0]["conserved"]
    assert "direction flips" in out.iloc[0]["verdict"]


def test_conservation_needs_enough_cancer_types():
    up = pd.DataFrame({"regulon": ["TF1"], "hedges_g": [1.0], "p_value": [0.001], "fdr": [0.001]})
    out = conserved_regulons({"LUAD": up, "BRCA": up}, min_cancer_types=3)
    assert not out.iloc[0]["conserved"] and "cancer-type-specific" in out.iloc[0]["verdict"]


# ------------------------------------------------------------------- modules


def test_preservation_separates_preserved_from_cohort_specific_modules():
    genes = [f"g{i}" for i in range(150)]

    def cohort(seed, coherent):
        rng = np.random.default_rng(seed)
        X = rng.normal(0, 1, (150, 50))
        X[:20] += rng.normal(0, 1, 50) * 2.5                  # module A: always coherent
        if coherent:
            X[20:40] += rng.normal(0, 1, 50) * 2.5            # module B: reference only
        return pd.DataFrame(X, index=genes)

    frame = module_preservation(
        cohort(1, True), cohort(2, False),
        {"A": genes[:20], "B": genes[20:40], "noise": genes[100:120]},
        n_permutations=60, rng=np.random.default_rng(3),
    ).set_index("module")
    assert frame.loc["A", "preserved"]
    assert not frame.loc["B", "preserved"]
    assert not frame.loc["noise", "preserved"]
    assert "not report this module as pan-cancer" in frame.loc["B", "interpretation"]


# ------------------------------------------------------------------ dynamics


def test_two_velocity_tools_are_one_line_of_evidence():
    r = corroborate_direction("A->B", {"scvelo": "forward", "unitvelo": "forward"})
    assert not r.corroborated
    assert "one view" in r.reason


def test_independent_families_corroborate():
    r = corroborate_direction("A->B", {"scvelo": "forward", "palantir": "forward"})
    assert r.corroborated and r.direction == "forward"
    assert set(r.families_agreeing) == {"velocity", "pseudotime"}


def test_conflicting_views_report_nothing():
    r = corroborate_direction("A->B", {"scvelo": "forward", "palantir": "reverse",
                                       "moscot": "reverse"})
    assert not r.corroborated and r.direction is None and "disagree" in r.reason


def test_view_families_are_recognised():
    assert view_family("scVelo") == "velocity"
    assert view_family("monocle3") == "pseudotime"
    assert view_family("timepoint") == "realtime"


def test_invalid_direction_is_rejected():
    with pytest.raises(ValueError, match="forward\\|reverse\\|none"):
        corroborate_direction("A->B", {"scvelo": "upwards"})


# ----------------------------------------------------------------- alignment


@pytest.mark.parametrize("same_tech,full,distort,expected", [
    (True, True, False, "paste"),
    (True, False, False, "paste2"),
    (True, True, True, "stalign"),
    (False, True, False, "slat"),
])
def test_aligner_choice_follows_section_properties(same_tech, full, distort, expected):
    choice = select_aligner(
        same_technology=same_tech, full_overlap=full, nonlinear_distortion=distort
    )
    assert choice.method == expected


def test_partial_overlap_estimates_rather_than_assumes():
    assert select_aligner(same_technology=True, full_overlap=False).estimate_overlap


def test_bad_alignment_is_reported_as_unacceptable():
    source = np.zeros((20, 2))
    target = np.full((20, 2), 400.0)
    stats = check_alignment(source, target, list(range(20)), max_residual_um=100)
    assert not stats["acceptable"] and stats["median_residual_um"] > 100


def test_good_alignment_passes():
    rng = np.random.default_rng(0)
    source = rng.normal(0, 100, (30, 2))
    target = source + rng.normal(0, 5, (30, 2))
    assert check_alignment(source, target, list(range(30)))["acceptable"]


# ----------------------------------------------------------------- histology


def test_histology_names_the_organs_where_prediction_fails():
    rng = np.random.default_rng(0)
    n = 200
    truth = rng.normal(0, 1, (n, 10))
    pred = truth.copy()
    organs = np.array(["lung"] * 100 + ["brain"] * 100)
    pred[100:] = rng.normal(0, 1, (100, 10))          # brain prediction is noise
    patients = np.array([f"P{i%10}" for i in range(n)])
    frame = per_organ_performance(pred, truth, organs, patients).set_index("organ")
    assert frame.loc["lung", "usable"]
    assert not frame.loc["brain", "usable"]
    assert "exclude this organ" in frame.loc["brain", "verdict"]


def test_organs_with_too_few_patients_are_not_evaluable():
    truth = np.random.default_rng(0).normal(0, 1, (20, 5))
    frame = per_organ_performance(
        truth, truth, ["kidney"] * 20, ["P0"] * 20
    ).set_index("organ")
    assert not frame.loc["kidney", "usable"] and "not evaluable" in frame.loc["kidney", "verdict"]


# ------------------------------------------------------------------- predict


def _cv_fixture(n_patients=20, per_patient=25, seed=0):
    rng = np.random.default_rng(seed)
    patients = np.repeat([f"P{i:02d}" for i in range(n_patients)], per_patient)
    signal = np.repeat(rng.normal(0, 1, n_patients), per_patient)
    X = np.column_stack([signal + rng.normal(0, 0.4, len(patients)),
                         rng.normal(0, 1, len(patients))])
    y = signal + rng.normal(0, 0.3, len(patients))
    return X, y, patients


def _ridge(Xtr, ytr, Xte):
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-9        # fitted inside the fold
    A = np.c_[(Xtr - mu) / sd, np.ones(len(Xtr))]
    w = np.linalg.lstsq(A, ytr, rcond=None)[0]
    return np.c_[(Xte - mu) / sd, np.ones(len(Xte))] @ w


def _r2(yt, yp):
    ss_res = float(np.sum((yt - yp) ** 2))
    ss_tot = float(np.sum((yt - yt.mean()) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot else float("nan")


def test_patient_level_cv_runs_a_leakage_check_on_every_fold():
    X, y, patients = _cv_fixture()
    report = patient_level_cv(
        X, y, patients, fit_predict=_ridge, score=_r2, metric="r2", n_splits=5, seed=1
    )
    assert len(report.folds) == 5
    assert report.leakage_checks == 5
    assert all(f.n_test_patients >= 1 for f in report.folds)
    assert report.mean > 0.3


def test_buffer_removes_boundary_training_points_and_asserts_again():
    X, y, patients = _cv_fixture(n_patients=10, per_patient=20)
    coords = np.column_stack([np.arange(len(X)) * 50.0, np.zeros(len(X))])
    sections = patients.copy()
    report = patient_level_cv(
        X, y, patients, fit_predict=_ridge, score=_r2, metric="r2",
        sections=sections, coords=coords, buffer_um=100.0, n_splits=5, seed=2,
    )
    assert report.leakage_checks == 10, "one assertion before the buffer, one after"


def test_single_patient_cohort_is_refused():
    X, y, patients = _cv_fixture(n_patients=1, per_patient=50)
    with pytest.raises(LeakageError, match="one patient"):
        patient_level_cv(X, y, patients, fit_predict=_ridge, score=_r2)


# ----------------------------------------------------------------- translate


def test_signature_scoring_is_robust_to_one_dominant_gene():
    expr = pd.DataFrame(
        {"s1": [1.0, 1.0, 1000.0], "s2": [2.0, 2.0, 1.0], "s3": [3.0, 3.0, 500.0]},
        index=["A", "B", "HUGE"],
    )
    scores = score_signature(expr, ["A", "B", "HUGE"])
    assert np.isfinite(scores).all() and len(scores) == 3
    assert scores.argmax() != 0, "z-scoring must stop one gene dominating"


def test_signature_needs_at_least_one_present_gene():
    expr = pd.DataFrame({"s1": [1.0]}, index=["A"])
    with pytest.raises(ValueError, match="none of the signature genes"):
        score_signature(expr, ["X", "Y"])


def test_target_without_a_falsifier_is_refused():
    bad = TargetCandidate(gene="TGFB1", niche="CAF", spatial_effect=0.9, conservation=0.8)
    with pytest.raises(ValueError, match="falsifier"):
        prioritise_targets([bad])


def test_missing_evidence_is_penalised_not_ignored():
    complete = TargetCandidate(
        "A", "n", 0.8, 0.8, dependency=0.8, survival=0.8, tractability=0.8,
        falsifier="blockade in co-culture fails to restore cytotoxicity",
    )
    sparse = TargetCandidate(
        "B", "n", 0.8, 0.8, falsifier="same assay",
    )
    frame = prioritise_targets([complete, sparse]).set_index("gene")
    assert frame.loc["A", "priority_score"] > frame.loc["B", "priority_score"]
    assert frame.loc["B", "missing_evidence"] != ""
    assert frame.loc["A", "rank"] == 1


def test_weights_must_sum_to_one():
    c = TargetCandidate("A", "n", 0.5, 0.5, falsifier="x")
    with pytest.raises(ValueError, match="sum to 1"):
        prioritise_targets([c], weights={"spatial_effect": 0.5, "conservation": 0.2})


def test_scores_outside_the_unit_interval_are_refused():
    c = TargetCandidate("A", "n", 1.8, 0.5, falsifier="x")
    with pytest.raises(ValueError, match="scaled to \\[0, 1\\]"):
        prioritise_targets([c])


def test_unadjusted_survival_analysis_is_refused():
    """An unadjusted HR for a microenvironment signature mostly restates stage."""
    from panspatial.translate import survival_association

    with pytest.raises(ValueError, match="required and must include"):
        survival_association([0.1] * 50, [10] * 50, [1] * 50)


def test_survival_requires_the_named_covariates():
    from panspatial.translate import survival_association

    cov = pd.DataFrame({"age": np.arange(50)})       # no stage column
    with pytest.raises(ValueError, match="required covariate"):
        survival_association([0.1] * 50, [10] * 50, [1] * 50, cov)


def test_survival_declines_on_too_few_events():
    from panspatial.translate import survival_association

    cov = pd.DataFrame({"stage": ["I"] * 50})
    with pytest.raises(ValueError, match="not interpretable"):
        survival_association(
            np.linspace(0, 1, 50), np.linspace(1, 100, 50), [1] * 3 + [0] * 47, cov
        )
