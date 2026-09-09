"""CLI: per-section torus-shift niche test for one sample.

    python -m panspatial.nk.run_niche --input nk_states.h5ad --anchor mCAF --out out.tsv
"""

from __future__ import annotations

import argparse
import logging

import numpy as np

from panspatial.nk.nk_spatial import per_section_niche_test, require_single_cell_resolution

log = logging.getLogger("panspatial.nk.run_niche")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, help="AnnData with spatial coords and cell labels")
    ap.add_argument("--anchor", required=True, help="anchor cell class, e.g. mCAF")
    ap.add_argument("--out", required=True)
    ap.add_argument("--focal", default="NK")
    ap.add_argument("--label-key", default="cell_type")
    ap.add_argument("--section-key", default="section_id")
    ap.add_argument("--patient-key", default="patient")
    ap.add_argument("--platform-key", default="platform_st")
    ap.add_argument("--coords-key", default="spatial_um")
    ap.add_argument("--n-perms", type=int, default=999)
    ap.add_argument("--min-focal", type=int, default=30)
    ap.add_argument("--alternative", default="less", choices=["less", "greater", "two-sided"])
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-7s %(message)s")
    import anndata as ad

    adata = ad.read_h5ad(args.input)

    platform = str(adata.uns.get(args.platform_key, adata.obs.get(args.platform_key, ["unknown"])[0]))
    require_single_cell_resolution(platform)

    if args.coords_key not in adata.obsm:
        raise SystemExit(
            f"{args.input} has no obsm[{args.coords_key!r}]. Coordinates must be converted to "
            "microns at load time -- pixel coordinates make every distance meaningless."
        )
    coords = np.asarray(adata.obsm[args.coords_key], dtype=float)[:, :2]

    frame = per_section_niche_test(
        coords,
        adata.obs[args.label_key].to_numpy(),
        adata.obs[args.section_key].to_numpy(),
        adata.obs[args.patient_key].to_numpy(),
        focal_label=args.focal,
        anchor_label=args.anchor,
        n_perms=args.n_perms,
        alternative=args.alternative,
        min_focal=args.min_focal,
        seed=args.seed,
    )
    frame["platform"] = platform
    frame["anchor"] = args.anchor
    frame.to_csv(args.out, sep="\t", index=False)
    log.info("wrote %s (%d sections, %d skipped)", args.out, len(frame),
             len(frame.attrs["skipped_sections"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
