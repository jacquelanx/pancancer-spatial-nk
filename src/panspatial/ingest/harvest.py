"""CLI: build the GEO candidate manifest.  python -m panspatial.ingest.harvest --help"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

from panspatial.ingest.geo_sra import build_manifest


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tumor-types", required=True, help="comma-separated tumor abbreviations")
    ap.add_argument("--out", required=True, help="output TSV")
    ap.add_argument("--retmax", type=int, default=500)
    ap.add_argument("--cache-dir", default=".cache/entrez")
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-7s %(message)s")
    tumors = [t.strip() for t in args.tumor_types.split(",") if t.strip()]
    frame = build_manifest(tumors, retmax=args.retmax)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out, sep="\t", index=False)
    print(f"wrote {out} ({len(frame)} candidate series)")
    if not frame.empty:
        usable = frame[frame["ST_is_single_cell"] & (frame["Pairing"] == "paired")]
        print(
            f"  {len(usable)} paired series on single-cell-resolution platforms "
            f"(the only ones eligible for NK-resolved spatial claims)"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
