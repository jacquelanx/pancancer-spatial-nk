"""CLI: pool per-section niche tests to patients, then meta-analyze across the cohort.

    python -m panspatial.nk.run_meta --inputs a.tsv b.tsv --anchor mCAF --out meta.json

This is the only step permitted to emit a cross-patient number, and the only place the
family-wide FDR correction is applied.
"""

from __future__ import annotations

import argparse
import json
import logging
from dataclasses import asdict
from pathlib import Path

import pandas as pd

from panspatial.nk.nk_spatial import cohort_meta_analysis
from panspatial.stats.spatial_stats import benjamini_hochberg

log = logging.getLogger("panspatial.nk.run_meta")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", nargs="+", required=True, help="per-sample niche TSVs")
    ap.add_argument("--anchor", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument(
        "--family-size",
        type=int,
        required=True,
        help="total number of tests in the study-wide family (anchors x cancer types), so "
        "the correction reflects everything tested, not just this anchor",
    )
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-7s %(message)s")

    frames = [pd.read_csv(p, sep="\t") for p in args.inputs]
    frames = [f for f in frames if len(f)]
    if not frames:
        raise SystemExit(
            f"no sections survived filtering for anchor {args.anchor!r}; nothing to pool"
        )
    combined = pd.concat(frames, ignore_index=True)
    combined.attrs["focal_label"] = "NK"
    combined.attrs["anchor_label"] = args.anchor

    result = cohort_meta_analysis(combined)
    section_fdr = benjamini_hochberg(
        combined["p_value"].to_numpy(), family_size=max(args.family_size, len(combined))
    )

    payload = {
        "anchor": args.anchor,
        "focal": "NK",
        "n_sections": int(len(combined)),
        "n_patients": int(result.k),
        "family_size": int(args.family_size),
        "meta": asdict(result),
        "interpretation": _interpret(result),
        "per_section_fdr_min": float(section_fdr.min()),
        "platforms": sorted(combined.get("platform", pd.Series(dtype=str)).unique().tolist()),
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(payload, indent=2))
    log.info("wrote %s | %s", args.out, result)
    return 0


def _interpret(result) -> str:
    """A one-line reading, phrased at the strength the evidence actually supports."""
    if result.p_value >= 0.05:
        return "no consistent spatial association across patients"
    direction = "closer to" if result.estimate < 0 else "farther from"
    strength = "single-cohort" if result.is_driven_by_one_cohort else "conserved"
    return (
        f"NK cells are {direction} this structure than an autocorrelation-preserving null "
        f"({strength} across {result.k} patients, I2={result.i2:.0%}). This is spatial "
        f"co-localization, not evidence of signalling."
    )


if __name__ == "__main__":
    raise SystemExit(main())
