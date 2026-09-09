"""CLI: gate NK cells from an AnnData and write the identity audit report."""

from __future__ import annotations

import argparse
import json
import logging
from dataclasses import asdict
from pathlib import Path

from panspatial.nk.nk_spatial import extract_nk_cells


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--sample-key", default="sample_id")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-7s %(message)s")
    import anndata as ad

    adata = ad.read_h5ad(args.input)
    nk, report = extract_nk_cells(adata, sample_key=args.sample_key, seed=args.seed)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    nk.write_h5ad(args.out)
    Path(args.report).write_text(json.dumps(asdict(report), indent=2))
    print(report.summary())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
