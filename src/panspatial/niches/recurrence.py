"""Matching niches across samples, and testing whether they recur.

A niche found in one tumour is an observation. A niche recurring across patients and cancer
types is a finding. The step between them is the one most often skipped: niches are called
per section, so "the same niche" in two samples means two clusters whose cellular
compositions are similar, and that correspondence has to be established before any
cross-cohort statement is possible.

Matching here is by composition cosine similarity with a Hungarian assignment, and
recurrence is reported at patient and cancer-type level with an exact binomial test against
a stated null.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

import numpy as np

log = logging.getLogger("panspatial.niches.recurrence")


def composition_matrix(
    labels: Sequence[Any], cell_types: Sequence[Any], types: Sequence[str] | None = None
) -> tuple[np.ndarray, list[str], list[str]]:
    """Cell-type composition of each niche in one section.

    Returns ``(matrix (n_niches, n_types), niche_names, type_names)`` with rows summing to 1.
    """
    labels = np.asarray(labels, dtype=object)
    cell_types = np.asarray(cell_types, dtype=object)
    if labels.shape != cell_types.shape:
        raise ValueError("labels and cell_types must describe the same cells")
    niches = sorted({str(x) for x in labels})
    tnames = sorted({str(x) for x in cell_types}) if types is None else list(types)
    t_index = {t: i for i, t in enumerate(tnames)}

    mat = np.zeros((len(niches), len(tnames)), dtype=float)
    for i, nic in enumerate(niches):
        mask = labels.astype(str) == nic
        for ct in cell_types[mask]:
            j = t_index.get(str(ct))
            if j is not None:
                mat[i, j] += 1.0
    totals = mat.sum(axis=1, keepdims=True)
    mat = np.divide(mat, totals, out=np.zeros_like(mat), where=totals > 0)
    return mat, niches, tnames


def _cosine(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    an = a / np.clip(np.linalg.norm(a, axis=1, keepdims=True), 1e-12, None)
    bn = b / np.clip(np.linalg.norm(b, axis=1, keepdims=True), 1e-12, None)
    return an @ bn.T


@dataclass
class NicheMatch:
    sample: str
    local_niche: str
    reference_niche: str | None
    similarity: float
    matched: bool


def match_niches(
    reference: np.ndarray,
    reference_names: Sequence[str],
    sample_compositions: Mapping[str, tuple[np.ndarray, Sequence[str]]],
    *,
    min_similarity: float = 0.7,
) -> list[NicheMatch]:
    """Assign each sample's niches to reference niches by composition similarity.

    Assignment is one-to-one within a sample (Hungarian), so two local niches cannot both
    claim the same reference niche. A local niche whose best match falls below
    ``min_similarity`` stays unmatched — sample-specific structure is a real result, and
    forcing it onto a reference niche inflates apparent recurrence.
    """
    from scipy.optimize import linear_sum_assignment

    reference = np.asarray(reference, dtype=float)
    if reference.ndim != 2:
        raise ValueError(f"reference must be (n_niches, n_types), got {reference.shape}")
    if len(reference_names) != reference.shape[0]:
        raise ValueError("reference_names must have one entry per reference niche")

    out: list[NicheMatch] = []
    for sample, (comp, names) in sample_compositions.items():
        comp = np.asarray(comp, dtype=float)
        if comp.shape[1] != reference.shape[1]:
            raise ValueError(
                f"{sample}: composition has {comp.shape[1]} cell types but the reference "
                f"has {reference.shape[1]}; harmonise the cell-type vocabulary first"
            )
        sim = _cosine(comp, reference)
        rows, cols = linear_sum_assignment(-sim)
        assigned = {int(r): int(c) for r, c in zip(rows, cols)}
        for i, local in enumerate(names):
            j = assigned.get(i)
            s = float(sim[i, j]) if j is not None else float(sim[i].max())
            ok = j is not None and s >= min_similarity
            out.append(NicheMatch(
                sample=str(sample), local_niche=str(local),
                reference_niche=str(reference_names[j]) if ok else None,
                similarity=round(s, 4), matched=ok,
            ))
    n_matched = sum(m.matched for m in out)
    log.info("matched %d/%d local niches to the reference at cosine >= %.2f",
             n_matched, len(out), min_similarity)
    return out


@dataclass
class RecurrenceRow:
    niche: str
    n_patients_present: int
    n_patients_total: int
    n_cancer_types_present: int
    n_cancer_types_total: int
    patient_fraction: float
    p_value: float
    conserved: bool
    verdict: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "niche": self.niche,
            "n_patients_present": self.n_patients_present,
            "n_patients_total": self.n_patients_total,
            "patient_fraction": round(self.patient_fraction, 4),
            "n_cancer_types_present": self.n_cancer_types_present,
            "n_cancer_types_total": self.n_cancer_types_total,
            "p_value": self.p_value,
            "conserved": self.conserved,
            "verdict": self.verdict,
        }


def recurrence_test(
    matches: Sequence[NicheMatch],
    patient_of: Mapping[str, str],
    cancer_type_of: Mapping[str, str],
    *,
    null_prevalence: float = 0.25,
    min_cancer_types: int = 3,
    min_patient_fraction: float = 0.5,
) -> Any:
    """Test each reference niche for recurrence at patient and cancer-type level.

    The unit is the patient, never the section: a patient contributing four sections
    contributes one observation. ``null_prevalence`` is the prevalence a niche would have to
    beat to be interesting and must be stated — a niche present in 30% of patients is only
    "recurrent" relative to an explicit expectation.

    A niche is called conserved only if it clears both the patient-fraction threshold and
    appears in at least ``min_cancer_types`` cancer types, so an effect concentrated in one
    tumour type cannot be reported as pan-cancer.
    """
    import pandas as pd
    from scipy.stats import binomtest

    if not 0.0 < null_prevalence < 1.0:
        raise ValueError(f"null_prevalence must be in (0, 1), got {null_prevalence}")

    all_patients = {patient_of[m.sample] for m in matches if m.sample in patient_of}
    all_cancers = {cancer_type_of[m.sample] for m in matches if m.sample in cancer_type_of}
    if not all_patients:
        raise ValueError("no sample in `matches` could be mapped to a patient")

    by_niche: dict[str, tuple[set[str], set[str]]] = {}
    for m in matches:
        if not m.matched or m.reference_niche is None:
            continue
        pats, cans = by_niche.setdefault(m.reference_niche, (set(), set()))
        if m.sample in patient_of:
            pats.add(patient_of[m.sample])
        if m.sample in cancer_type_of:
            cans.add(cancer_type_of[m.sample])

    rows = []
    for niche, (pats, cans) in sorted(by_niche.items()):
        k, n = len(pats), len(all_patients)
        frac = k / n
        p = binomtest(k, n, null_prevalence, alternative="greater").pvalue
        enough_types = len(cans) >= min_cancer_types
        conserved = bool(frac >= min_patient_fraction and enough_types and p < 0.05)
        if conserved:
            verdict = (f"conserved: present in {k}/{n} patients across {len(cans)} cancer types")
        elif not enough_types:
            verdict = (f"cancer-type-restricted: only {len(cans)} cancer type(s) "
                       f"(< {min_cancer_types}); report as such, not as pan-cancer")
        else:
            verdict = (f"not recurrent: {k}/{n} patients ({frac:.0%}) does not clear "
                       f"{min_patient_fraction:.0%} at p < 0.05 against a {null_prevalence:.0%} null")
        rows.append(RecurrenceRow(
            niche, k, n, len(cans), len(all_cancers), frac, float(p), conserved, verdict,
        ))
    frame = pd.DataFrame([r.to_dict() for r in rows])
    if not frame.empty:
        frame = frame.sort_values("patient_fraction", ascending=False).reset_index(drop=True)
        frame.attrs["null_prevalence"] = null_prevalence
        log.info("%d/%d niches called conserved", int(frame["conserved"].sum()), len(frame))
    return frame
