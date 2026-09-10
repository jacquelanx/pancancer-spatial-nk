"""Matching per-patient malignant programs into recurrent meta-programs.

cNMF run on one patient's malignant cells returns programs specific to that patient. The
pan-cancer question is which programs recur, and answering it means matching gene
signatures across patients — by overlap, since the numeric factors are not comparable.

The matching is deliberately conservative: a meta-program must be supported by several
patients and, to be called pan-cancer, by several cancer types.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

import numpy as np

log = logging.getLogger("panspatial.malignant.programs")


def jaccard(a: set[str], b: set[str]) -> float:
    if not a and not b:
        return 0.0
    return len(a & b) / len(a | b)


@dataclass
class MetaProgram:
    name: str
    members: list[tuple[str, str]] = field(default_factory=list)   # (patient, program)
    core_genes: list[str] = field(default_factory=list)
    cancer_types: set[str] = field(default_factory=set)

    @property
    def n_patients(self) -> int:
        return len({p for p, _ in self.members})

    def to_dict(self) -> dict[str, Any]:
        return {
            "meta_program": self.name,
            "n_patients": self.n_patients,
            "n_programs": len(self.members),
            "n_cancer_types": len(self.cancer_types),
            "cancer_types": ",".join(sorted(self.cancer_types)),
            "core_genes": ",".join(self.core_genes[:30]),
            "n_core_genes": len(self.core_genes),
        }


def match_meta_programs(
    programs: Mapping[str, Mapping[str, Sequence[str]]],
    *,
    cancer_type_of: Mapping[str, str] | None = None,
    min_jaccard: float = 0.15,
    min_patients: int = 3,
    core_fraction: float = 0.5,
) -> list[MetaProgram]:
    """Cluster per-patient gene programs into recurrent meta-programs by gene overlap.

    Args:
        programs: ``{patient: {program_name: [genes]}}`` — typically the top-N genes of each
            cNMF factor.
        min_jaccard: overlap required to join two programs. Single-linkage agglomeration is
            used, which is permissive; the patient and cancer-type thresholds downstream do
            the filtering.
        core_fraction: a gene is a core gene of a meta-program if it appears in at least
            this fraction of the member programs.

    Programs from the same patient are never merged with each other: within one patient,
    cNMF factors are already distinct by construction, so merging them would inflate the
    apparent size of a meta-program.
    """
    if not 0.0 < min_jaccard <= 1.0:
        raise ValueError(f"min_jaccard must be in (0, 1], got {min_jaccard}")

    items: list[tuple[str, str, set[str]]] = [
        (patient, name, {str(g).upper() for g in genes})
        for patient, progs in programs.items()
        for name, genes in progs.items()
        if genes
    ]
    if not items:
        return []

    parent = list(range(len(items)))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i: int, j: int) -> None:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[max(ri, rj)] = min(ri, rj)

    for i in range(len(items)):
        for j in range(i + 1, len(items)):
            if items[i][0] == items[j][0]:
                continue  # never merge two programs from the same patient
            if jaccard(items[i][2], items[j][2]) >= min_jaccard:
                union(i, j)

    clusters: dict[int, list[int]] = {}
    for i in range(len(items)):
        clusters.setdefault(find(i), []).append(i)

    metas: list[MetaProgram] = []
    for members in sorted(clusters.values(), key=lambda m: -len({items[i][0] for i in m})):
        patients = {items[i][0] for i in members}
        if len(patients) < min_patients:
            continue
        counts: dict[str, int] = {}
        for i in members:
            for g in items[i][2]:
                counts[g] = counts.get(g, 0) + 1
        threshold = core_fraction * len(members)
        core = sorted([g for g, c in counts.items() if c >= threshold],
                      key=lambda g: (-counts[g], g))
        meta = MetaProgram(
            name=f"MP{len(metas) + 1}",
            members=[(items[i][0], items[i][1]) for i in members],
            core_genes=core,
            cancer_types={cancer_type_of[items[i][0]] for i in members
                          if cancer_type_of and items[i][0] in cancer_type_of},
        )
        metas.append(meta)

    log.info(
        "%d meta-programs from %d patient programs (>= %d patients, Jaccard >= %.2f)",
        len(metas), len(items), min_patients, min_jaccard,
    )
    return metas


def meta_program_table(metas: Sequence[MetaProgram], *, min_cancer_types: int = 3) -> Any:
    """Tidy frame with an explicit pan-cancer verdict per meta-program."""
    import pandas as pd

    rows = []
    for m in metas:
        d = m.to_dict()
        d["pan_cancer"] = len(m.cancer_types) >= min_cancer_types
        d["verdict"] = (
            f"recurrent across {len(m.cancer_types)} cancer types"
            if d["pan_cancer"]
            else f"restricted to {len(m.cancer_types) or 'unannotated'} cancer type(s); "
                 "not a pan-cancer program"
        )
        rows.append(d)
    return pd.DataFrame(rows)
