"""CLI: assemble the verified-results table.

This file is the contract behind ``--results-verified``. A manuscript prompt may quote a
number only if it appears here, which is what keeps a plausible-sounding invented effect
size out of a draft.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import pandas as pd

log = logging.getLogger("panspatial.report")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--meta", nargs="+", required=True, help="niche meta-analysis JSONs")
    ap.add_argument("--reliability", nargs="*", default=[], help="deconvolution reliability CSVs")
    ap.add_argument("--stability", nargs="*", default=[], help="NK state stability JSONs")
    ap.add_argument("--summary", required=True)
    ap.add_argument("--table", required=True)
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-7s %(message)s")

    rows = []
    for path in args.meta:
        payload = json.loads(Path(path).read_text())
        meta = payload["meta"]
        rows.append(
            {
                "comparison": f"NK vs {payload['anchor']}",
                "estimate_log_ratio": round(meta["estimate"], 4),
                "ci_low": round(meta["ci_low"], 4),
                "ci_high": round(meta["ci_high"], 4),
                "p_value": meta["p_value"],
                "n_patients": meta["k"],
                "n_sections": payload["n_sections"],
                "i2": round(meta["i2"], 3),
                "loo_min": round(meta["loo_min"], 4),
                "loo_max": round(meta["loo_max"], 4),
                "single_cohort_driven": bool(
                    (min(abs(meta["loo_min"]), abs(meta["loo_max"])) / meta["estimate"] < 0.5)
                    if meta["estimate"]
                    else False
                ),
                "family_size": payload["family_size"],
                "platforms": ",".join(payload["platforms"]),
                "interpretation": payload["interpretation"],
            }
        )
    table = pd.DataFrame(rows).sort_values("p_value").reset_index(drop=True)
    Path(args.table).parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.table, sep="\t", index=False)

    unreliable: list[str] = []
    for path in args.reliability:
        frame = pd.read_csv(path)
        unreliable += frame.loc[~frame["reliable"], "cell_type"].tolist()

    stabilities = [json.loads(Path(p).read_text()).get("bootstrap_ari", 0.0) for p in args.stability]

    lines = [
        "# Verified results",
        "",
        "Every number a manuscript prompt is permitted to quote appears in "
        f"`{Path(args.table).name}`. Anything absent must be written as [TBD].",
        "",
        f"- Comparisons tested: {len(table)}",
        f"- Significant at p < 0.05 before family correction: {int((table['p_value'] < 0.05).sum())}",
        f"- Significant AND driven by a single cohort (report as single-cohort, not "
        f"conserved): "
        f"{int((table['single_cohort_driven'] & (table['p_value'] < 0.05)).sum())}",
    ]
    if stabilities:
        lines.append(
            f"- NK state cluster stability (bootstrap ARI): "
            f"min {min(stabilities):.2f}, median {pd.Series(stabilities).median():.2f}"
        )
        if min(stabilities) < 0.6:
            lines.append(
                "  - At least one sample has unstable NK states; describe those as a "
                "continuum of programs, not discrete states."
            )
    if unreliable:
        lines.append(
            f"- Cell types that failed the deconvolution benchmark and may not carry a "
            f"conclusion: {', '.join(sorted(set(unreliable)))}"
        )
    lines += ["", "## Comparisons", "", _markdown_table(table)]

    Path(args.summary).write_text("\n".join(lines) + "\n")
    log.info("wrote %s and %s", args.summary, args.table)
    return 0


def _markdown_table(frame: pd.DataFrame) -> str:
    """Minimal markdown table, so the report has no optional-dependency footgun."""
    header = "| " + " | ".join(frame.columns) + " |"
    rule = "| " + " | ".join("---" for _ in frame.columns) + " |"
    rows = [
        "| " + " | ".join(_fmt(v) for v in row) + " |"
        for row in frame.itertuples(index=False, name=None)
    ]
    return "\n".join([header, rule, *rows])


def _fmt(value) -> str:
    if isinstance(value, float):
        return f"{value:.3g}"
    return str(value).replace("|", "\\|")


if __name__ == "__main__":
    raise SystemExit(main())
