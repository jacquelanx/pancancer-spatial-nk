"""Cell-type annotation with an explicit ambiguity class.

Every marker panel has cells that score for two things at once. Forcing them into one label
manufactures clean structure; the populations that sit in the gaps — NKT and γδ T between NK
and T, ILC1 without a cytotoxic program, CD16⁺ monocytes sharing FCGR3A with NK — are
reported here rather than absorbed.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

import numpy as np

from panspatial.modules.base import (
    AnalysisModule,
    Evidence,
    ModuleContext,
    ModuleError,
    ModuleResult,
    PlatformClass,
    ReplicationUnit,
    register,
    require,
)

log = logging.getLogger("panspatial.annotate")

UNASSIGNED = "unassigned"
AMBIGUOUS = "ambiguous"


@dataclass(frozen=True)
class MarkerPanel:
    """Positive and negative markers for one cell type, with the populations excluded."""

    label: str
    positive: tuple[str, ...]
    negative: tuple[str, ...] = ()
    discriminating: tuple[str, ...] = ()   # the markers that actually separate it
    excludes: tuple[str, ...] = ()         # populations the negatives are there to remove

    def genes(self) -> set[str]:
        return set(self.positive) | set(self.negative) | set(self.discriminating)


# A starting ontology. Deliberately small and explicit; extend per cohort.
DEFAULT_PANELS: tuple[MarkerPanel, ...] = (
    MarkerPanel("T_CD8", ("CD3D", "CD3E", "CD8A", "GZMK"), ("KLRF1", "NCR1"), ("CD8A", "CD8B")),
    MarkerPanel("T_CD4", ("CD3D", "CD3E", "CD4", "IL7R"), ("CD8A", "KLRF1"), ("CD4",)),
    MarkerPanel("Treg", ("FOXP3", "IL2RA", "CTLA4", "IKZF2"), (), ("FOXP3",)),
    MarkerPanel("NK", ("NCAM1", "KLRD1", "NKG7", "GNLY", "PRF1"), ("CD3D", "CD3E", "TRAC", "LYZ"),
                ("KLRF1", "NCR1"), ("NKT", "gamma-delta T", "ILC1", "CD16+ monocyte")),
    MarkerPanel("B", ("MS4A1", "CD79A", "CD19"), (), ("MS4A1", "CD79A")),
    MarkerPanel("Plasma", ("MZB1", "JCHAIN", "XBP1"), ("MS4A1",), ("MZB1",)),
    MarkerPanel("Macrophage", ("CD68", "CD163", "MRC1", "LYZ"), ("CD3E",), ("CD68", "CD163")),
    MarkerPanel("cDC1", ("CLEC9A", "XCR1", "BATF3"), (), ("CLEC9A", "XCR1")),
    MarkerPanel("Neutrophil", ("FCGR3B", "CSF3R", "S100A8"), (), ("FCGR3B",)),
    MarkerPanel("mCAF", ("COL1A1", "ACTA2", "TAGLN", "POSTN"), ("PTPRC",), ("POSTN", "ACTA2")),
    MarkerPanel("iCAF", ("IL6", "CXCL12", "PDGFRA", "HAS1"), ("PTPRC", "ACTA2"), ("IL6", "CXCL12")),
    MarkerPanel("Endothelial", ("PECAM1", "VWF", "CDH5"), ("PTPRC",), ("PECAM1", "VWF")),
    MarkerPanel("Pericyte", ("RGS5", "NOTCH3", "PDGFRB"), (), ("RGS5",)),
    MarkerPanel("Epithelial", ("EPCAM", "KRT8", "KRT18"), ("PTPRC",), ("EPCAM",)),
)


@dataclass
class AmbiguityReport:
    """The counts that belong in a supplementary table, not in a footnote."""

    n_cells: int
    counts: dict[str, int] = field(default_factory=dict)
    n_ambiguous: int = 0
    n_unassigned: int = 0
    missing_markers: dict[str, list[str]] = field(default_factory=dict)
    panels_without_discriminating: list[str] = field(default_factory=list)

    @property
    def fraction_resolved(self) -> float:
        unresolved = self.n_ambiguous + self.n_unassigned
        return 1.0 - unresolved / self.n_cells if self.n_cells else 0.0

    def summary(self) -> str:
        top = ", ".join(f"{k}={v}" for k, v in sorted(
            self.counts.items(), key=lambda kv: -kv[1])[:6])
        return (
            f"{self.n_cells} cells, {100*self.fraction_resolved:.1f}% resolved "
            f"(ambiguous={self.n_ambiguous}, unassigned={self.n_unassigned}); {top}"
        )


def assign_with_ambiguity(
    scores: np.ndarray,
    labels: Sequence[str],
    *,
    min_margin: float = 0.1,
    min_score: float = 0.0,
    exclusion: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Assign each cell to its top-scoring label, or to an explicit non-answer.

    A cell becomes ``unassigned`` when no label clears ``min_score``, and ``ambiguous`` when
    the top two labels are within ``min_margin`` or an exclusion program also fires. Both
    are real answers; forcing a call is not.

    Returns ``(assignment, margin)``.
    """
    scores = np.asarray(scores, dtype=float)
    if scores.ndim != 2:
        raise ValueError(f"scores must be (n_cells, n_labels), got {scores.shape}")
    if scores.shape[1] != len(labels):
        raise ValueError(
            f"scores has {scores.shape[1]} columns but {len(labels)} labels were given"
        )
    if scores.shape[1] < 2:
        raise ValueError("need at least two candidate labels to measure a margin")

    order = np.argsort(-scores, axis=1)
    rows = np.arange(len(scores))
    top, second = order[:, 0], order[:, 1]
    best = scores[rows, top]
    margin = best - scores[rows, second]

    assignment = np.asarray(labels, dtype=object)[top].astype(object)
    assignment[margin < min_margin] = AMBIGUOUS
    assignment[best <= min_score] = UNASSIGNED
    if exclusion is not None:
        excluded = np.asarray(exclusion, dtype=bool)
        if excluded.shape[0] != len(scores):
            raise ValueError("exclusion must have one entry per cell")
        assignment[excluded & (assignment != UNASSIGNED)] = AMBIGUOUS
    return assignment, margin


def score_panels(
    adata: Any,
    panels: Sequence[MarkerPanel] = DEFAULT_PANELS,
    *,
    seed: int = 0,
) -> tuple[np.ndarray, list[str], dict[str, list[str]], list[str]]:
    """Score every panel over an AnnData; report which markers were absent.

    Returns ``(scores, labels, missing_by_label, panels_lacking_discriminating_markers)``.
    """
    sc = require("scanpy", extra="spatial", purpose="marker scoring")
    present_index = {g.upper() for g in adata.var_names}

    cols, labels, missing, weak = [], [], {}, []
    for panel in panels:
        found = [g for g in panel.positive if g.upper() in present_index]
        absent = sorted(set(panel.positive) - set(found))
        if absent:
            missing[panel.label] = absent
        if panel.discriminating and not any(
            g.upper() in present_index for g in panel.discriminating
        ):
            weak.append(panel.label)
            log.warning(
                "%s: none of its discriminating markers (%s) are in the panel; calls will "
                "rest on markers shared with %s and are provisional",
                panel.label, ", ".join(panel.discriminating),
                ", ".join(panel.excludes) or "other populations",
            )
        key = f"_score_{panel.label}"
        if found:
            sc.tl.score_genes(adata, found, score_name=key, random_state=seed)
            cols.append(np.asarray(adata.obs[key]))
        else:
            cols.append(np.zeros(adata.n_obs))
        labels.append(panel.label)
    return np.column_stack(cols), labels, missing, weak


def annotate(
    adata: Any,
    panels: Sequence[MarkerPanel] = DEFAULT_PANELS,
    *,
    min_margin: float = 0.1,
    seed: int = 0,
    key_added: str = "cell_type",
) -> AmbiguityReport:
    """Annotate an AnnData in place and return the ambiguity report."""
    scores, labels, missing, weak = score_panels(adata, panels, seed=seed)
    assignment, margin = assign_with_ambiguity(scores, labels, min_margin=min_margin)
    adata.obs[key_added] = assignment
    adata.obs[f"{key_added}_margin"] = margin

    counts = {l: int(np.sum(assignment == l)) for l in labels}
    report = AmbiguityReport(
        n_cells=int(adata.n_obs), counts=counts,
        n_ambiguous=int(np.sum(assignment == AMBIGUOUS)),
        n_unassigned=int(np.sum(assignment == UNASSIGNED)),
        missing_markers=missing, panels_without_discriminating=weak,
    )
    log.info(report.summary())
    return report


@register
class Annotation(AnalysisModule):
    name = "annotate"
    version = "0.1.0"
    phase = "2"
    description = "Marker-panel and reference-based cell typing with an explicit ambiguity class."
    evidence = Evidence.INFERRED
    replication_unit = ReplicationUnit.CELL
    supported_platforms = frozenset(
        {PlatformClass.SINGLE_CELL, PlatformClass.DISSOCIATED}
    )
    optional_dependencies = ("scanpy", "celltypist")

    def run(self, ctx: ModuleContext, *, adata: Any = None, **inputs: Any) -> ModuleResult:
        import pandas as pd

        if adata is None:
            raise ModuleError("annotate.run requires adata=")
        report = annotate(
            adata, inputs.get("panels", DEFAULT_PANELS),
            min_margin=ctx.param("min_margin", 0.1), seed=ctx.seed,
            key_added=ctx.param("key_added", "cell_type"),
        )
        rows = [{"label": k, "n_cells": v,
                 "fraction": round(v / max(report.n_cells, 1), 5),
                 "missing_markers": ",".join(report.missing_markers.get(k, [])),
                 "discriminating_markers_present": k not in report.panels_without_discriminating}
                for k, v in report.counts.items()]
        rows += [
            {"label": AMBIGUOUS, "n_cells": report.n_ambiguous,
             "fraction": round(report.n_ambiguous / max(report.n_cells, 1), 5),
             "missing_markers": "", "discriminating_markers_present": True},
            {"label": UNASSIGNED, "n_cells": report.n_unassigned,
             "fraction": round(report.n_unassigned / max(report.n_cells, 1), 5),
             "missing_markers": "", "discriminating_markers_present": True},
        ]
        notes = ["Ambiguous and unassigned are reported, not redistributed."]
        if report.panels_without_discriminating:
            notes.append(
                "Provisional calls (no discriminating marker in panel): "
                + ", ".join(report.panels_without_discriminating)
            )
        return self.result(
            ctx, pd.DataFrame(rows), notes=notes,
            fraction_resolved=report.fraction_resolved,
            reliable_entities=[l for l in report.counts
                               if l not in report.panels_without_discriminating],
        )
