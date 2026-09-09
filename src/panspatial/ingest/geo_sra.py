"""GEO/SRA harvesting for paired scRNA-seq + spatial transcriptomics cancer datasets.

The classification layer is deliberately pure and regex-driven, not LLM-driven: the same
GSE must classify identically on every re-run, and each label carries the substring that
produced it so a curator can audit it. Negations ("not Visium", "non-spatial") are excluded
explicitly, which a bare substring test gets wrong.

Network access is confined to :class:`EntrezClient`, which throttles to NCBI's published
limits (3 req/s anonymous, 10 req/s with ``NCBI_API_KEY``) and caches responses on disk, so
re-running the manifest build does not re-hit the API.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
import re
import time
from dataclasses import asdict, dataclass, field
from enum import Enum
from pathlib import Path
from typing import Iterable, Sequence

log = logging.getLogger("panspatial.ingest")

# --------------------------------------------------------------------- classification

# (label, pattern). Order matters: the first match wins, so put specific before general.
ST_PLATFORM_PATTERNS: list[tuple[str, str]] = [
    ("Visium HD", r"visium\s*hd"),
    ("Visium", r"\bvisium\b|10x\s+spatial|spatial\s+gene\s+expression"),
    ("Xenium", r"\bxenium\b"),
    ("CosMx", r"\bcosmx\b|nanostring\s+spatial\s+molecular"),
    ("MERFISH", r"\bmerfish\b|merscope"),
    ("Slide-seq", r"slide-?seq(?:\s*v?2)?"),
    ("GeoMx", r"\bgeomx\b|digital\s+spatial\s+profil"),
    ("STARmap", r"\bstarmap\b"),
    ("seqFISH", r"seq-?fish\+?"),
    ("Stereo-seq", r"stereo-?seq"),
    ("ST (legacy)", r"spatial\s+transcriptomic"),
]

SC_PLATFORM_PATTERNS: list[tuple[str, str]] = [
    ("10x v4", r"10x.{0,20}\bv4\b|chromium.{0,20}\bv4\b|gem-?x"),
    ("10x v3", r"10x.{0,20}\bv3\b|chromium.{0,20}\bv3\b|single\s+cell\s+3'\s*v3"),
    ("10x v2", r"10x.{0,20}\bv2\b|chromium.{0,20}\bv2\b"),
    ("10x (version unstated)", r"\b10x\b|chromium|cellranger"),
    ("Smart-seq2", r"smart-?seq\s*2"),
    ("Smart-seq3", r"smart-?seq\s*3"),
    ("BD Rhapsody", r"bd\s+rhapsody"),
    ("Drop-seq", r"drop-?seq"),
    ("inDrop", r"in-?drop"),
]

TUMOR_PATTERNS: list[tuple[str, str]] = [
    ("LUAD", r"\bluad\b|lung\s+adenocarcinoma"),
    ("LUSC", r"\blusc\b|lung\s+squamous"),
    ("NSCLC", r"\bnsclc\b|non-?small\s+cell\s+lung"),
    ("BRCA", r"\bbrca\b|breast\s+(?:cancer|carcinoma|tumou?r)|\btnbc\b"),
    ("CRC", r"\bcrc\b|colorectal|colon\s+(?:cancer|adenocarcinoma)|\bcoad\b|\bread\b"),
    ("GBM", r"\bgbm\b|glioblastoma|glioma"),
    ("HNSCC", r"\bhnscc\b|head\s+and\s+neck"),
    ("PDAC", r"\bpdac\b|pancreatic\s+(?:ductal\s+)?adenocarcinoma|\bpaad\b"),
    ("HCC", r"\bhcc\b|hepatocellular"),
    ("SKCM", r"\bskcm\b|melanoma"),
    ("RCC", r"\brcc\b|renal\s+cell\s+carcinoma|\bkirc\b"),
    ("STAD", r"\bstad\b|gastric\s+(?:cancer|adenocarcinoma)"),
    ("OV", r"\bov\b|ovarian\s+(?:cancer|carcinoma)"),
    ("PRAD", r"\bprad\b|prostate\s+(?:cancer|adenocarcinoma)"),
    ("BLCA", r"\bblca\b|bladder\s+(?:cancer|carcinoma)|urothelial"),
    ("ESCA", r"\besca\b|esophageal|oesophageal"),
    ("CESC", r"\bcesc\b|cervical\s+(?:cancer|carcinoma)"),
    ("UCEC", r"\bucec\b|endometrial"),
]

TREATMENT_PATTERNS: list[tuple[str, str]] = [
    ("post-ICB", r"anti-?pd-?[l1]|pembrolizumab|nivolumab|atezolizumab|durvalumab|"
                 r"ipilimumab|immune\s+checkpoint|\bicb\b|\bici\b"),
    ("post-chemo", r"neoadjuvant|chemotherap|post-?treatment|after\s+treatment"),
    ("treatment-naive", r"treatment[-\s]?na[iï]ve|untreated|therapy[-\s]?na[iï]ve|"
                        r"chemotherapy[-\s]?na[iï]ve"),
]

# platform -> (nominal resolution in um, captures single cells)
ST_RESOLUTION: dict[str, tuple[float | None, bool]] = {
    "Visium": (55.0, False),
    "Visium HD": (2.0, False),  # 2 um bins; sub-cellular bins, still not segmented cells
    "Xenium": (None, True),
    "CosMx": (None, True),
    "MERFISH": (None, True),
    "seqFISH": (None, True),
    "STARmap": (None, True),
    "Slide-seq": (10.0, False),
    "Stereo-seq": (0.5, False),
    "GeoMx": (None, False),  # region-of-interest, not a grid
    "ST (legacy)": (100.0, False),
}

# A match preceded by any of these is a negation and must not count.
_NEGATION = re.compile(r"(?:not|non|no|without|excluding|lack(?:ing|s)?\s+of|rather\s+than)[\s-]+$",
                       re.IGNORECASE)


class PairingLevel(str, Enum):
    """How strong the link between the modalities is. These are not interchangeable."""

    PAIRED = "paired"          # both modalities in the same series
    MATCHED = "matched"        # same patient identifiers recoverable across series
    UNPAIRED = "unpaired"      # one modality only, or no recoverable link


def _match(text: str, patterns: Sequence[tuple[str, str]]) -> tuple[str | None, str | None]:
    """First non-negated pattern match; returns (label, matched evidence)."""
    if not text:
        return None, None
    for label, pattern in patterns:
        for m in re.finditer(pattern, text, flags=re.IGNORECASE):
            if _NEGATION.search(text[max(0, m.start() - 24) : m.start()]):
                log.debug("skipping negated match %r for %s", m.group(0), label)
                continue
            start = max(0, m.start() - 32)
            return label, text[start : m.end() + 32].strip()
    return None, None


def classify_platform_st(text: str) -> tuple[str | None, str | None]:
    """Spatial platform label plus the substring that produced it."""
    return _match(text, ST_PLATFORM_PATTERNS)


def classify_platform_sc(text: str) -> tuple[str | None, str | None]:
    """Single-cell platform label plus the substring that produced it."""
    return _match(text, SC_PLATFORM_PATTERNS)


def classify_tumor_type(text: str) -> tuple[str | None, str | None]:
    return _match(text, TUMOR_PATTERNS)


def classify_treatment(text: str) -> tuple[str, str | None]:
    label, evidence = _match(text, TREATMENT_PATTERNS)
    return (label or "unstated"), evidence


def st_resolution(platform: str | None) -> tuple[float | None, bool]:
    """(nominal resolution um, is_single_cell) for a spatial platform.

    ``is_single_cell`` gates every NK-specific spatial claim downstream. NK cells are 0.5-3%
    of a tumor; on a 55 um Visium spot holding 1-10 cells they are recoverable only as a
    deconvolution weight, and that weight is not separable from CD8 T cells at realistic
    depth. Spot platforms are usable for niche context, not for NK-resolved conclusions.
    """
    if platform is None:
        return None, False
    return ST_RESOLUTION.get(platform, (None, False))


# --------------------------------------------------------------- structure verification

_FILE_SIGNALS: dict[str, str] = {
    "matrix": r"matrix\.mtx|_matrix|counts?\.(?:csv|tsv|txt|h5|h5ad)|filtered_feature_bc_matrix|\.h5ad",
    "barcodes": r"barcodes?\.(?:tsv|csv|txt)",
    "features": r"(?:features|genes)\.(?:tsv|csv|txt)",
    "spatial_coords": r"tissue_positions|spatial|coordinates?|centroids|cell_boundaries|_pos\.",
    "image": r"\.(?:tif|tiff|png|jpg|jpeg)(?:\.gz)?$|hires_image|lowres_image",
    "scalefactors": r"scalefactors",
}


@dataclass
class StructureReport:
    """What a series' supplementary directory actually contains."""

    present: dict[str, bool]
    n_files: int
    missing: list[str] = field(default_factory=list)

    @property
    def has_expression(self) -> bool:
        return self.present.get("matrix", False)

    @property
    def has_spatial(self) -> bool:
        return self.present.get("spatial_coords", False)

    @property
    def is_analysis_ready(self) -> bool:
        return self.has_expression and self.has_spatial


def verify_file_structure(filenames: Iterable[str]) -> StructureReport:
    """Check a supplementary file listing for the components an analysis needs.

    Missing components are recorded, not used to silently drop the series -- a study whose
    counts live in a single ``.h5ad`` legitimately has no ``barcodes.tsv``, and a study with
    processed-only deposits is still worth flagging for a data request.
    """
    names = [n.lower() for n in filenames]
    joined = "\n".join(names)
    present = {key: bool(re.search(pattern, joined)) for key, pattern in _FILE_SIGNALS.items()}
    missing = [k for k in ("matrix", "spatial_coords") if not present[k]]
    return StructureReport(present=present, n_files=len(names), missing=missing)


# ------------------------------------------------------------------------- records


@dataclass
class GeoRecord:
    Dataset_ID: str
    Title: str
    Tumor_Type: str | None
    Platform_scRNA: str | None
    Platform_ST: str | None
    ST_resolution_um: float | None
    ST_is_single_cell: bool
    Pairing: str
    Treatment_Status: str
    Sample_Count: int
    Patient_Count: int | None
    Has_Open_Access: bool
    Download_URL: str
    Structure_OK: bool
    Structure_Missing: str
    Evidence: str

    @classmethod
    def from_summary(
        cls,
        summary: dict,
        *,
        structure: StructureReport | None = None,
        patient_count: int | None = None,
    ) -> "GeoRecord":
        """Build a record from a GEO DocSummary (``Entrez.esummary`` over db='gds')."""
        accession = str(summary.get("Accession", "")).strip()
        title = str(summary.get("title", "")).strip()
        text = " ".join(
            str(summary.get(k, "")) for k in ("title", "summary", "taxon", "gdsType", "suppFile")
        )

        st_platform, st_evidence = classify_platform_st(text)
        sc_platform, sc_evidence = classify_platform_sc(text)
        tumor, tumor_evidence = classify_tumor_type(text)
        treatment, treatment_evidence = classify_treatment(text)
        resolution, single_cell = st_resolution(st_platform)

        if st_platform and sc_platform:
            pairing = PairingLevel.PAIRED
        elif st_platform or sc_platform:
            pairing = PairingLevel.UNPAIRED
        else:
            pairing = PairingLevel.UNPAIRED

        structure = structure or StructureReport(present={}, n_files=0, missing=["unchecked"])
        evidence = " | ".join(
            f"{tag}: {snippet}"
            for tag, snippet in (
                ("ST", st_evidence),
                ("scRNA", sc_evidence),
                ("tumor", tumor_evidence),
                ("treatment", treatment_evidence),
            )
            if snippet
        )
        return cls(
            Dataset_ID=accession,
            Title=title,
            Tumor_Type=tumor,
            Platform_scRNA=sc_platform,
            Platform_ST=st_platform,
            ST_resolution_um=resolution,
            ST_is_single_cell=single_cell,
            Pairing=pairing.value,
            Treatment_Status=treatment,
            Sample_Count=int(summary.get("n_samples", 0) or 0),
            Patient_Count=patient_count,
            Has_Open_Access=bool(summary.get("suppFile")),
            Download_URL=supplementary_url(accession),
            Structure_OK=structure.is_analysis_ready,
            Structure_Missing=",".join(structure.missing),
            Evidence=evidence,
        )


def supplementary_url(accession: str) -> str:
    """GEO supplementary directory URL for a GSE accession (HTTPS, not FTP)."""
    m = re.fullmatch(r"GSE(\d+)", accession.strip(), flags=re.IGNORECASE)
    if not m:
        return ""
    digits = m.group(1)
    stub = f"GSE{digits[:-3]}nnn" if len(digits) > 3 else "GSEnnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/GSE{digits}/suppl/"


# ------------------------------------------------------------------- network access


class EntrezClient:
    """Throttled, cached wrapper around Bio.Entrez.

    NCBI allows 3 requests/s anonymously and 10/s with an API key; exceeding it gets the
    caller blocked, which in a long cohort sweep means losing the whole run.
    """

    def __init__(
        self,
        email: str | None = None,
        api_key: str | None = None,
        cache_dir: str | Path = ".cache/entrez",
    ) -> None:
        try:
            from Bio import Entrez
        except ImportError as exc:  # pragma: no cover - environment dependent
            raise RuntimeError("biopython is required for GEO queries: pip install biopython") from exc

        self.email = email or os.environ.get("NCBI_EMAIL", "")
        if not self.email:
            raise ValueError(
                "NCBI requires a contact email for E-utilities; set NCBI_EMAIL or pass email="
            )
        self.api_key = api_key or os.environ.get("NCBI_API_KEY")
        Entrez.email = self.email
        if self.api_key:
            Entrez.api_key = self.api_key
        self._entrez = Entrez
        self._min_interval = 1 / 10 if self.api_key else 1 / 3
        self._last_call = 0.0
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        log.info(
            "Entrez ready (%s key) -> %.0f req/s, cache=%s",
            "with" if self.api_key else "no",
            1 / self._min_interval,
            self.cache_dir,
        )

    def _throttle(self) -> None:
        elapsed = time.monotonic() - self._last_call
        if elapsed < self._min_interval:
            time.sleep(self._min_interval - elapsed)
        self._last_call = time.monotonic()

    def _cache_path(self, key: str) -> Path:
        return self.cache_dir / f"{hashlib.sha256(key.encode()).hexdigest()[:20]}.json"

    def _cached(self, key: str, fetch, *, max_attempts: int = 4):
        path = self._cache_path(key)
        if path.is_file():
            log.debug("cache hit: %s", key)
            return json.loads(path.read_text())
        for attempt in range(max_attempts):
            self._throttle()
            try:
                result = fetch()
                path.write_text(json.dumps(result, default=str))
                return result
            except Exception as exc:  # Entrez raises HTTPError / URLError / RuntimeError
                if attempt == max_attempts - 1:
                    log.error("giving up on %s after %d attempts: %s", key, max_attempts, exc)
                    raise
                delay = 2**attempt
                log.warning("Entrez call failed (%s); retrying in %ds", exc, delay)
                time.sleep(delay)

    def search_series(self, query: str, retmax: int = 500) -> list[str]:
        """Return GDS UIDs for a GEO DataSets query."""

        def _fetch():
            handle = self._entrez.esearch(db="gds", term=query, retmax=retmax)
            try:
                return dict(self._entrez.read(handle))
            finally:
                handle.close()

        return list(self._cached(f"esearch:{query}:{retmax}", _fetch).get("IdList", []))

    def summaries(self, uids: Sequence[str], batch_size: int = 100) -> list[dict]:
        """Fetch DocSummaries for GDS UIDs, in batches."""
        out: list[dict] = []
        for start in range(0, len(uids), batch_size):
            batch = list(uids[start : start + batch_size])

            def _fetch(batch=batch):
                handle = self._entrez.esummary(db="gds", id=",".join(batch))
                try:
                    return [dict(rec) for rec in self._entrez.read(handle)]
                finally:
                    handle.close()

            out.extend(self._cached(f"esummary:{','.join(batch)}", _fetch))
            log.info("fetched summaries %d-%d of %d", start + 1, start + len(batch), len(uids))
        return out


def default_query(tumor_terms: Sequence[str]) -> str:
    """GEO query for human cancer series carrying spatial and/or single-cell data."""
    tumors = " OR ".join(f'"{t}"[All Fields]' for t in tumor_terms)
    modality = (
        '("spatial transcriptomics"[All Fields] OR "Visium"[All Fields] OR "Xenium"[All Fields] '
        'OR "MERFISH"[All Fields] OR "CosMx"[All Fields] OR "Slide-seq"[All Fields])'
    )
    return (
        f'("Homo sapiens"[Organism]) AND (GSE[Entry Type]) AND ({tumors}) AND {modality}'
    )


def build_manifest(
    tumor_terms: Sequence[str],
    *,
    client: EntrezClient | None = None,
    retmax: int = 500,
    query: str | None = None,
):
    """Search GEO and return a pandas DataFrame of candidate datasets.

    Nothing is filtered out here. Series that fail structure verification stay in the table
    with the reason recorded, because "processed data only" is a curation decision, not a
    parsing failure.
    """
    import pandas as pd

    client = client or EntrezClient()
    q = query or default_query(tumor_terms)
    log.info("GEO query: %s", q)
    uids = client.search_series(q, retmax=retmax)
    log.info("%d series matched", len(uids))
    records = [asdict(GeoRecord.from_summary(s)) for s in client.summaries(uids)]
    frame = pd.DataFrame(records)
    if not frame.empty:
        frame = frame.sort_values(
            ["ST_is_single_cell", "Pairing", "Sample_Count"], ascending=[False, True, False]
        ).reset_index(drop=True)
    log.info(
        "manifest: %d rows, %d paired, %d single-cell-resolution ST",
        len(frame),
        int((frame.get("Pairing") == "paired").sum()) if not frame.empty else 0,
        int(frame.get("ST_is_single_cell", pd.Series(dtype=bool)).sum()) if not frame.empty else 0,
    )
    return frame
