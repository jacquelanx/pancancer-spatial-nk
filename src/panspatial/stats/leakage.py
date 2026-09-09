"""Cross-validation splitting for spatial models, and the leakage assertions to go with it.

A model that predicts an NK state from its spatial neighborhood will score ~0.95 AUC and
mean nothing if the splits are made over spots. Two spots 100 um apart share the same
tissue microenvironment, the same patient, the same batch, and on Visium sometimes the same
cells; one in train and one in test is memorization, not generalization.

The rule enforced here: **split on patients, then add a spatial buffer.** The patient split
handles biological and batch non-independence; the buffer handles the residual case where a
section is split internally on purpose (e.g. region-held-out designs within one patient).
"""

from __future__ import annotations

import logging
from typing import Sequence

import numpy as np
from scipy.spatial import cKDTree

log = logging.getLogger("panspatial.stats.leakage")


class LeakageError(AssertionError):
    """Raised when a train/test split shares information across the boundary."""


def patient_blocked_splits(
    patient_ids: Sequence[str],
    *,
    n_splits: int = 5,
    rng: np.random.Generator | None = None,
    stratify_by: Sequence[str] | None = None,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Grouped k-fold where whole patients move together between folds.

    Args:
        patient_ids: one label per observation (cell or spot). Observations from a patient
            are never split across folds.
        n_splits: number of folds; capped at the number of distinct patients.
        rng: seeded generator.
        stratify_by: optional per-observation label (typically cancer type) held constant
            within patient. Folds are balanced across its levels so no fold is missing a
            cancer type -- which otherwise makes fold scores incomparable.

    Returns:
        List of (train_idx, test_idx) index arrays.
    """
    patient_ids = np.asarray(patient_ids)
    if patient_ids.ndim != 1:
        raise ValueError("patient_ids must be 1-D, one entry per observation")
    rng = rng or np.random.default_rng(0)

    patients = np.unique(patient_ids)
    if len(patients) < 2:
        raise LeakageError(
            f"only {len(patients)} patient(s); cross-validation over one patient measures "
            "memorization of that patient, not generalization"
        )
    n_splits = min(n_splits, len(patients))
    if n_splits < 2:
        raise ValueError("n_splits must be at least 2")

    if stratify_by is not None:
        strat = np.asarray(stratify_by)
        if len(strat) != len(patient_ids):
            raise ValueError("stratify_by must be one entry per observation")
        patient_stratum = {}
        for p in patients:
            levels = np.unique(strat[patient_ids == p])
            if len(levels) > 1:
                raise ValueError(
                    f"patient {p!r} spans multiple strata {list(levels)}; stratify_by must be "
                    "constant within patient"
                )
            patient_stratum[p] = levels[0]
        fold_of: dict[str, int] = {}
        for level in np.unique(list(patient_stratum.values())):
            members = np.array([p for p in patients if patient_stratum[p] == level])
            rng.shuffle(members)
            for i, p in enumerate(members):
                fold_of[p] = i % n_splits
    else:
        shuffled = patients.copy()
        rng.shuffle(shuffled)
        fold_of = {p: i % n_splits for i, p in enumerate(shuffled)}

    assignment = np.array([fold_of[p] for p in patient_ids])
    splits = []
    for fold in range(n_splits):
        test = np.flatnonzero(assignment == fold)
        train = np.flatnonzero(assignment != fold)
        if len(test) == 0 or len(train) == 0:
            continue
        splits.append((train, test))
    return splits


def spatial_buffer_filter(
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    coords: np.ndarray,
    section_ids: Sequence[str],
    *,
    buffer_um: float,
) -> np.ndarray:
    """Drop training observations within ``buffer_um`` of any test observation, per section.

    Only within-section proximity matters -- coordinates from different sections share no
    frame of reference, so a distance between them is meaningless. Returns the filtered
    training indices.
    """
    coords = np.asarray(coords, dtype=float)
    section_ids = np.asarray(section_ids)
    if coords.shape[0] != len(section_ids):
        raise ValueError("coords and section_ids must describe the same observations")
    if buffer_um <= 0:
        return np.asarray(train_idx)

    train_idx = np.asarray(train_idx)
    test_idx = np.asarray(test_idx)
    drop = np.zeros(len(train_idx), dtype=bool)
    for section in np.unique(section_ids[test_idx]):
        test_here = test_idx[section_ids[test_idx] == section]
        train_mask = section_ids[train_idx] == section
        if not train_mask.any():
            continue
        tree = cKDTree(coords[test_here])
        distances, _ = tree.query(coords[train_idx[train_mask]], k=1)
        drop[train_mask] |= distances < buffer_um

    if drop.any():
        log.info(
            "spatial buffer (%.0f um) removed %d/%d training observations",
            buffer_um,
            drop.sum(),
            len(train_idx),
        )
    return train_idx[~drop]


def assert_no_leakage(
    train_idx: np.ndarray,
    test_idx: np.ndarray,
    patient_ids: Sequence[str],
    *,
    coords: np.ndarray | None = None,
    section_ids: Sequence[str] | None = None,
    buffer_um: float = 0.0,
    allow_shared_patients: bool = False,
) -> None:
    """Raise :class:`LeakageError` if a split leaks. Call this in the CV loop, not once.

    Checks, in order of severity: index overlap, patient overlap, and -- when coordinates
    are supplied -- training points closer than ``buffer_um`` to a test point in the same
    section.

    Args:
        allow_shared_patients: set only for a deliberate region-held-out design, where one
            section is split spatially on purpose. It disables the patient check, so a
            positive ``buffer_um`` becomes mandatory -- otherwise nothing is guarding the
            split at all. Such a design measures within-patient spatial generalization,
            which is a weaker claim than generalization to a new patient; say so in the
            manuscript.
    """
    if allow_shared_patients and buffer_um <= 0:
        raise ValueError(
            "allow_shared_patients=True without a spatial buffer disables every leakage "
            "guard; pass buffer_um > 0 (typically several times the platform's spot pitch)"
        )
    train_idx = np.asarray(train_idx)
    test_idx = np.asarray(test_idx)
    patient_ids = np.asarray(patient_ids)

    overlap = np.intersect1d(train_idx, test_idx)
    if overlap.size:
        raise LeakageError(f"{overlap.size} observation(s) appear in both train and test")

    shared = np.intersect1d(patient_ids[train_idx], patient_ids[test_idx])
    if shared.size and not allow_shared_patients:
        raise LeakageError(
            f"{shared.size} patient(s) in both train and test: {list(shared[:5])}. "
            "Split on patients, not on cells or spots."
        )

    if coords is not None and buffer_um > 0:
        if section_ids is None:
            raise ValueError("section_ids is required for the spatial buffer check")
        remaining = spatial_buffer_filter(
            train_idx, test_idx, coords, section_ids, buffer_um=buffer_um
        )
        if len(remaining) != len(train_idx):
            raise LeakageError(
                f"{len(train_idx) - len(remaining)} training observation(s) lie within "
                f"{buffer_um} um of a test observation in the same section. Apply "
                "spatial_buffer_filter to the training set before fitting."
            )
