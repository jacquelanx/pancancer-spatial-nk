"""Tests for GEO metadata classification and structure verification (no network)."""

import pytest

from panspatial.ingest.geo_sra import (
    GeoRecord,
    classify_platform_sc,
    classify_platform_st,
    classify_treatment,
    classify_tumor_type,
    st_resolution,
    supplementary_url,
    verify_file_structure,
)


@pytest.mark.parametrize(
    "text,expected",
    [
        ("Spatial profiling with 10x Visium", "Visium"),
        ("Visium HD 2um binned data", "Visium HD"),
        ("Xenium in situ 5K panel", "Xenium"),
        ("NanoString CosMx SMI", "CosMx"),
        ("MERSCOPE imaging of tumors", "MERFISH"),
        ("Slide-seqV2 pucks", "Slide-seq"),
        ("Stereo-seq chips", "Stereo-seq"),
    ],
)
def test_st_platforms_are_recognized(text, expected):
    assert classify_platform_st(text)[0] == expected


def test_visium_hd_wins_over_plain_visium():
    """Order matters: the specific pattern must be tried before the general one."""
    assert classify_platform_st("10x Visium HD spatial gene expression")[0] == "Visium HD"


@pytest.mark.parametrize(
    "text",
    [
        "bulk RNA-seq, not Visium",
        "non-spatial transcriptomic profiling",
        "scRNA-seq without Xenium validation",
    ],
)
def test_negated_mentions_do_not_produce_a_platform(text):
    assert classify_platform_st(text)[0] is None


def test_classification_returns_auditable_evidence():
    label, evidence = classify_platform_st("Paired scRNA-seq and 10x Visium of PDAC")
    assert label == "Visium"
    assert "Visium" in evidence


@pytest.mark.parametrize(
    "text,expected",
    [
        ("Chromium Single Cell 3' v3", "10x v3"),
        ("10x Genomics v2 chemistry", "10x v2"),
        ("Smart-seq2 plate-based", "Smart-seq2"),
        ("processed with CellRanger", "10x (version unstated)"),
    ],
)
def test_sc_platforms_are_recognized(text, expected):
    assert classify_platform_sc(text)[0] == expected


@pytest.mark.parametrize(
    "text,expected",
    [
        ("lung adenocarcinoma cohort", "LUAD"),
        ("triple negative breast cancer (TNBC)", "BRCA"),
        ("colorectal liver metastases", "CRC"),
        ("glioblastoma multiforme", "GBM"),
        ("head and neck squamous cell carcinoma", "HNSCC"),
        ("hepatocellular carcinoma", "HCC"),
    ],
)
def test_tumor_types_are_recognized(text, expected):
    assert classify_tumor_type(text)[0] == expected


@pytest.mark.parametrize(
    "text,expected",
    [
        ("samples collected after anti-PD-1 therapy", "post-ICB"),
        ("neoadjuvant chemotherapy cohort", "post-chemo"),
        ("treatment-naive primary tumors", "treatment-naive"),
        ("resected primary tumors", "unstated"),
    ],
)
def test_treatment_status_is_recognized(text, expected):
    assert classify_treatment(text)[0] == expected


def test_platform_resolution_drives_the_single_cell_flag():
    assert st_resolution("Visium") == (55.0, False)
    assert st_resolution("Xenium")[1] is True
    assert st_resolution("CosMx")[1] is True
    assert st_resolution(None) == (None, False)
    assert st_resolution("SomeNewPlatform") == (None, False)


def test_supplementary_url_bucketing():
    assert supplementary_url("GSE181919").endswith("/GSE181nnn/GSE181919/suppl/")
    assert supplementary_url("GSE999").endswith("/GSEnnn/GSE999/suppl/")
    assert supplementary_url("not-an-accession") == ""


def test_structure_verification_finds_visium_layout():
    report = verify_file_structure(
        [
            "GSE1_filtered_feature_bc_matrix.h5",
            "GSE1_tissue_positions_list.csv.gz",
            "GSE1_tissue_hires_image.png",
            "GSE1_scalefactors_json.json",
        ]
    )
    assert report.is_analysis_ready
    assert report.present["image"]
    assert report.missing == []


def test_structure_verification_records_what_is_missing():
    report = verify_file_structure(["GSE2_counts.csv.gz", "GSE2_metadata.txt"])
    assert not report.is_analysis_ready
    assert report.missing == ["spatial_coords"]
    assert report.has_expression


def test_record_marks_a_two_modality_series_as_paired():
    record = GeoRecord.from_summary(
        {
            "Accession": "GSE123456",
            "title": "Paired scRNA-seq (Chromium v3) and Xenium of treatment-naive PDAC",
            "summary": "Twelve patients profiled before therapy.",
            "n_samples": 24,
            "suppFile": "MTX,CSV",
        }
    )
    assert record.Pairing == "paired"
    assert record.Platform_ST == "Xenium"
    assert record.Platform_scRNA == "10x v3"
    assert record.Tumor_Type == "PDAC"
    assert record.Treatment_Status == "treatment-naive"
    assert record.ST_is_single_cell
    assert record.Sample_Count == 24
    assert "Xenium" in record.Evidence


def test_record_flags_spot_resolution_for_downstream_gating():
    record = GeoRecord.from_summary(
        {
            "Accession": "GSE654321",
            "title": "10x Visium of colorectal cancer with matched scRNA-seq",
            "n_samples": 8,
            "suppFile": "TAR",
        }
    )
    assert record.ST_is_single_cell is False
    assert record.ST_resolution_um == 55.0
