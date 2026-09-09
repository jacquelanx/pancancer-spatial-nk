"""Spatial null models, cross-sample meta-analysis, and multiple-testing correction.

Two failure modes dominate spatial-omics statistics, and both are addressed here.

**Inflated significance from naive permutation.** Cell-type labels in tissue are strongly
spatially autocorrelated. Shuffling labels over the tissue destroys that autocorrelation,
so the null distribution is far tighter than reality and almost any co-localization test
comes out significant. :func:`torus_shift_permutation_test` instead translates the focal
point pattern rigidly on a torus: the internal structure of the pattern is preserved and
only its registration against the anchor pattern is randomized.

**Pseudoreplication.** With 10^5 cells per section, a cell-level test has an n large enough
to make any effect significant, but the cells within a patient are not independent draws.
The correct unit is the patient. :func:`random_effects_meta` combines per-sample effects
with a DerSimonian-Laird random-effects model and reports between-sample heterogeneity and
the leave-one-out range, which is what tells you whether a "pan-cancer" effect is real or
is one cohort in a trench coat.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Callable, Sequence

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import norm

log = logging.getLogger("panspatial.stats")

Statistic = Callable[[np.ndarray, np.ndarray], float]


# --------------------------------------------------------------------------- statistics


def mean_nn_distance(focal_xy: np.ndarray, anchor_xy: np.ndarray) -> float:
    """Mean distance from each focal point to its nearest anchor point.

    Coordinates must be in physical units (um). Mixing pixel coordinates from one platform
    with micron coordinates from another silently produces meaningless distances, which is
    why the callers in this package convert at load time.
    """
    focal_xy = np.asarray(focal_xy, dtype=float)
    anchor_xy = np.asarray(anchor_xy, dtype=float)
    if focal_xy.ndim != 2 or focal_xy.shape[1] != 2:
        raise ValueError(f"focal_xy must be (n, 2), got {focal_xy.shape}")
    if anchor_xy.ndim != 2 or anchor_xy.shape[1] != 2:
        raise ValueError(f"anchor_xy must be (n, 2), got {anchor_xy.shape}")
    if len(focal_xy) == 0 or len(anchor_xy) == 0:
        return float("nan")
    distances, _ = cKDTree(anchor_xy).query(focal_xy, k=1)
    return float(np.mean(distances))


def fraction_within(radius_um: float) -> Statistic:
    """Statistic factory: fraction of focal points with an anchor within ``radius_um``."""

    def _stat(focal_xy: np.ndarray, anchor_xy: np.ndarray) -> float:
        if len(focal_xy) == 0 or len(anchor_xy) == 0:
            return float("nan")
        distances, _ = cKDTree(np.asarray(anchor_xy, float)).query(np.asarray(focal_xy, float), k=1)
        return float(np.mean(distances <= radius_um))

    return _stat


# --------------------------------------------------------------------- permutation test


@dataclass(frozen=True)
class PermutationResult:
    observed: float
    null_mean: float
    null_sd: float
    effect: float
    z: float
    p_value: float
    n_perms: int
    n_focal: int
    n_anchor: int
    alternative: str

    def __str__(self) -> str:  # pragma: no cover - display only
        return (
            f"obs={self.observed:.2f} null={self.null_mean:.2f}+/-{self.null_sd:.2f} "
            f"z={self.z:.2f} p={self.p_value:.2g} (n_focal={self.n_focal}, "
            f"{self.n_perms} torus shifts, alternative={self.alternative})"
        )


def torus_shift_permutation_test(
    focal_xy: np.ndarray,
    anchor_xy: np.ndarray,
    *,
    statistic: Statistic = mean_nn_distance,
    alternative: str = "less",
    n_perms: int = 999,
    inside_tissue: Callable[[np.ndarray], np.ndarray] | None = None,
    min_retained: float = 0.5,
    bounds: tuple[float, float, float, float] | None = None,
    rng: np.random.Generator | None = None,
) -> PermutationResult:
    """Test focal-to-anchor spatial association against an autocorrelation-preserving null.

    The focal pattern is rigidly translated on the torus of the section's bounding box.
    This keeps the focal pattern's own clustering intact -- only its alignment with the
    anchor pattern is randomized -- so a significant result means "these two patterns are
    registered with each other", not merely "these cells are clustered".

    Args:
        focal_xy: (n, 2) coordinates of the focal population, in um, for ONE section.
        anchor_xy: (m, 2) coordinates of the reference structure, in um, same section.
        statistic: summary of the focal-anchor relationship. Lower must mean "closer" for
            ``alternative="less"``.
        alternative: ``"less"`` tests attraction (observed statistic smaller than null),
            ``"greater"`` tests exclusion, ``"two-sided"`` tests either.
        inside_tissue: predicate mapping (k, 2) coordinates to a boolean mask. Supply the
            section's tissue mask: without it, shifted points land in empty slide area and
            the null is biased towards large distances, which manufactures "attraction".
        min_retained: resample a shift if it leaves less than this fraction of focal points
            inside the tissue mask.
        bounds: (xmin, xmax, ymin, ymax) torus extent; defaults to the union bounding box.
        rng: seeded generator. Always pass one -- an unseeded permutation test is not
            reproducible and will not survive review.

    Returns:
        :class:`PermutationResult` with the observed statistic, the null distribution
        summary, and a permutation p-value with the standard +1 correction (never zero).
    """
    if alternative not in {"less", "greater", "two-sided"}:
        raise ValueError(f"alternative must be less|greater|two-sided, got {alternative!r}")
    focal_xy = np.asarray(focal_xy, dtype=float)
    anchor_xy = np.asarray(anchor_xy, dtype=float)
    if len(focal_xy) < 3:
        raise ValueError(
            f"only {len(focal_xy)} focal points in this section; a per-section spatial test "
            "on a handful of cells is noise. Drop the section and record it as dropped."
        )
    if rng is None:
        log.warning("torus_shift_permutation_test called without a seeded rng; using entropy")
        rng = np.random.default_rng()

    observed = statistic(focal_xy, anchor_xy)
    if not np.isfinite(observed):
        raise ValueError("observed statistic is not finite; check for empty populations")

    if bounds is None:
        stacked = np.vstack([focal_xy, anchor_xy])
        xmin, ymin = stacked.min(axis=0)
        xmax, ymax = stacked.max(axis=0)
    else:
        xmin, xmax, ymin, ymax = bounds
    width, height = xmax - xmin, ymax - ymin
    if width <= 0 or height <= 0:
        raise ValueError("degenerate bounding box: all points are collinear or identical")

    origin = np.array([xmin, ymin])
    extent = np.array([width, height])
    relative = focal_xy - origin

    null = np.empty(n_perms, dtype=float)
    rejected = 0
    max_rejects = 20 * n_perms
    i = 0
    while i < n_perms:
        shifted = (relative + rng.uniform(0.0, 1.0, size=2) * extent) % extent + origin
        if inside_tissue is not None:
            keep = np.asarray(inside_tissue(shifted), dtype=bool)
            if keep.mean() < min_retained:
                rejected += 1
                if rejected > max_rejects:
                    raise RuntimeError(
                        "torus shifts keep falling outside the tissue mask; the section is "
                        "probably a thin or non-convex fragment. Use a block bootstrap "
                        "restricted to the mask instead of a torus shift."
                    )
                continue
            shifted = shifted[keep]
        value = statistic(shifted, anchor_xy)
        if not np.isfinite(value):
            rejected += 1
            continue
        null[i] = value
        i += 1

    null_mean = float(np.mean(null))
    null_sd = float(np.std(null, ddof=1))
    if alternative == "less":
        p = (1.0 + np.sum(null <= observed)) / (n_perms + 1.0)
    elif alternative == "greater":
        p = (1.0 + np.sum(null >= observed)) / (n_perms + 1.0)
    else:
        centered = np.abs(null - null_mean)
        p = (1.0 + np.sum(centered >= abs(observed - null_mean))) / (n_perms + 1.0)

    return PermutationResult(
        observed=observed,
        null_mean=null_mean,
        null_sd=null_sd,
        effect=observed - null_mean,
        z=float((observed - null_mean) / null_sd) if null_sd > 0 else float("nan"),
        p_value=float(p),
        n_perms=n_perms,
        n_focal=len(focal_xy),
        n_anchor=len(anchor_xy),
        alternative=alternative,
    )


# ------------------------------------------------------------------------ meta-analysis


@dataclass(frozen=True)
class MetaAnalysisResult:
    estimate: float
    se: float
    ci_low: float
    ci_high: float
    z: float
    p_value: float
    tau2: float
    i2: float
    q: float
    k: int
    loo_min: float
    loo_max: float
    loo_driver: str | None

    @property
    def is_driven_by_one_cohort(self) -> bool:
        """True when dropping one unit halves the pooled estimate or flips its sign.

        A "pan-cancer" effect that fails this check is a single-cohort effect reported
        with a pan-cancer n, and must be described as such.
        """
        if not np.isfinite(self.estimate) or self.estimate == 0:
            return False
        extremes = np.array([self.loo_min, self.loo_max])
        closest_to_zero = float(extremes[np.argmin(np.abs(extremes))])
        return (closest_to_zero / self.estimate) < 0.5

    def __str__(self) -> str:  # pragma: no cover - display only
        return (
            f"{self.estimate:.3f} [{self.ci_low:.3f}, {self.ci_high:.3f}] "
            f"p={self.p_value:.2g}, k={self.k}, I2={self.i2:.0%}, "
            f"LOO range [{self.loo_min:.3f}, {self.loo_max:.3f}]"
        )


def random_effects_meta(
    effects: Sequence[float],
    variances: Sequence[float],
    labels: Sequence[str] | None = None,
) -> MetaAnalysisResult:
    """DerSimonian-Laird random-effects pooling of per-sample (or per-cohort) effects.

    This is the layer that makes a pan-cancer claim legitimate: each patient or cohort
    contributes one effect with its own uncertainty, between-study heterogeneity is
    estimated rather than assumed away, and the leave-one-out range exposes an effect that
    exists only because of one cohort.

    Args:
        effects: per-unit effect estimates (e.g. per-patient log ratio of observed to null
            mean distance). One value per patient or cohort -- never per cell.
        variances: sampling variance of each effect. Must be positive.
        labels: unit names, used to report which unit drives the pooled estimate.
    """
    y = np.asarray(effects, dtype=float)
    v = np.asarray(variances, dtype=float)
    if y.shape != v.shape:
        raise ValueError(f"effects {y.shape} and variances {v.shape} must have the same shape")
    finite = np.isfinite(y) & np.isfinite(v) & (v > 0)
    if finite.sum() != len(y):
        log.warning("dropping %d unit(s) with non-finite or non-positive variance", (~finite).sum())
    names = list(labels) if labels is not None else [f"unit{i}" for i in range(len(y))]
    y, v = y[finite], v[finite]
    names = [n for n, keep in zip(names, finite) if keep]

    k = len(y)
    if k < 2:
        raise ValueError(
            f"random-effects pooling needs at least 2 units, got {k}. With one sample there "
            "is no evidence about between-sample variability, so no pan-cancer claim."
        )

    pooled = _dl_pool(y, v)
    estimate, se, tau2, q = pooled

    loo = []
    for i in range(k):
        keep = np.ones(k, dtype=bool)
        keep[i] = False
        loo.append(_dl_pool(y[keep], v[keep])[0])
    loo_arr = np.asarray(loo)
    driver = names[int(np.argmax(np.abs(loo_arr - estimate)))] if k > 1 else None

    z = estimate / se if se > 0 else float("nan")
    return MetaAnalysisResult(
        estimate=float(estimate),
        se=float(se),
        ci_low=float(estimate - 1.959963985 * se),
        ci_high=float(estimate + 1.959963985 * se),
        z=float(z),
        p_value=float(2 * norm.sf(abs(z))) if np.isfinite(z) else float("nan"),
        tau2=float(tau2),
        i2=float(max(0.0, (q - (k - 1)) / q)) if q > 0 else 0.0,
        q=float(q),
        k=k,
        loo_min=float(loo_arr.min()),
        loo_max=float(loo_arr.max()),
        loo_driver=driver,
    )


def _dl_pool(y: np.ndarray, v: np.ndarray) -> tuple[float, float, float, float]:
    """Return (estimate, se, tau2, Q) for a DerSimonian-Laird random-effects fit."""
    w = 1.0 / v
    fixed = float(np.sum(w * y) / np.sum(w))
    q = float(np.sum(w * (y - fixed) ** 2))
    k = len(y)
    c = float(np.sum(w) - np.sum(w**2) / np.sum(w))
    tau2 = max(0.0, (q - (k - 1)) / c) if c > 0 else 0.0
    w_star = 1.0 / (v + tau2)
    estimate = float(np.sum(w_star * y) / np.sum(w_star))
    se = float(np.sqrt(1.0 / np.sum(w_star)))
    return estimate, se, tau2, q


# ------------------------------------------------------------------- multiple testing


def benjamini_hochberg(p_values: Sequence[float], family_size: int | None = None) -> np.ndarray:
    """BH-FDR adjusted p-values.

    ``family_size`` makes the correction family explicit and lets you correct against a
    family larger than the vector you happen to be holding -- the common way a pan-cancer
    screen understates its own multiplicity is by correcting each cancer type separately
    and reporting the union.
    """
    p = np.asarray(p_values, dtype=float)
    if p.ndim != 1:
        raise ValueError(f"p_values must be 1-D, got shape {p.shape}")
    if np.any((p < 0) | (p > 1)):
        raise ValueError("p_values must lie in [0, 1]")
    n = len(p)
    m = n if family_size is None else int(family_size)
    if m < n:
        raise ValueError(f"family_size ({m}) cannot be smaller than the number of tests ({n})")

    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * m / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(adjusted, 0.0, 1.0)
    return out
