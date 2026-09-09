"""CLI: score NK functional states and report cluster stability."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--stability", required=True)
    ap.add_argument("--resolution", type=float, default=0.6)
    ap.add_argument("--min-margin", type=float, default=0.1)
    ap.add_argument("--n-bootstrap", type=int, default=20)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-7s %(message)s")
    import anndata as ad

    from panspatial.nk.nk_spatial import score_nk_states

    adata = ad.read_h5ad(args.input)
    adata, stats = score_nk_states(
        adata,
        resolution=args.resolution,
        min_margin=args.min_margin,
        seed=args.seed,
        n_bootstrap=args.n_bootstrap,
    )
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out)
    Path(args.stability).write_text(json.dumps(stats, indent=2))
    print(json.dumps(stats, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
