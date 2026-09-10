"""Generic module runner: one entrypoint that dispatches to any registered module.

    python -m panspatial.modules.run --job jobs/niches_LUAD.json

A job file is the module name, the context, and the inputs. Keeping this generic means the
workflow gains a rule per analysis layer without a bespoke script per layer, and every run
writes the same three artifacts: the tidy table, a summary with provenance, and the notes
that must travel with the result.

Job file shape::

    {
      "module": "niches",
      "context": {"sample_id": "LUAD_cohort", "platform": "Xenium", "seed": 42,
                  "params": {"min_patient_fraction": 0.5}},
      "inputs":  {"sample_compositions": "results/niches/compositions.pkl",
                  "reference": "results/niches/reference.npy"}
    }

Inputs are resolved by file extension: ``.h5ad`` via anndata, ``.tsv``/``.csv``/``.parquet``
via pandas, ``.npy`` via numpy, ``.json`` as parsed JSON, ``.pkl`` via pickle. A value that
is not an existing path is passed through unchanged, so scalars can be given inline.
"""

from __future__ import annotations

import argparse
import json
import logging
import pickle
import sys
from pathlib import Path
from typing import Any

from panspatial.modules import ModuleContext, ModuleError, get_module

log = logging.getLogger("panspatial.modules.run")


def load_input(value: Any) -> Any:
    """Resolve one job input: a path is loaded, anything else passes through."""
    if not isinstance(value, str):
        return value
    path = Path(value)
    if not path.is_file():
        return value

    suffix = path.suffix.lower()
    if suffix == ".h5ad":
        import anndata
        return anndata.read_h5ad(path)
    if suffix in {".tsv", ".txt"}:
        import pandas as pd
        return pd.read_csv(path, sep="\t")
    if suffix == ".csv":
        import pandas as pd
        return pd.read_csv(path)
    if suffix == ".parquet":
        import pandas as pd
        return pd.read_parquet(path)
    if suffix == ".npy":
        import numpy as np
        return np.load(path, allow_pickle=False)
    if suffix == ".json":
        return json.loads(path.read_text())
    if suffix in {".pkl", ".pickle"}:
        # Only ever load pickles this pipeline wrote; they are intermediates, not inputs
        # from outside the run.
        with path.open("rb") as fh:
            return pickle.load(fh)
    raise ModuleError(f"unsupported input type {suffix!r} for {path}")


def run_job(spec: dict[str, Any], outdir: Path) -> dict[str, Any]:
    name = spec.get("module")
    if not name:
        raise ModuleError("job file must name a `module`")
    module = get_module(name)

    ctx_spec = dict(spec.get("context", {}))
    ctx = ModuleContext(
        sample_id=ctx_spec.get("sample_id", "unnamed"),
        platform=ctx_spec.get("platform"),
        patient_id=ctx_spec.get("patient_id"),
        cancer_type=ctx_spec.get("cancer_type"),
        seed=int(ctx_spec.get("seed", 42)),
        params=dict(ctx_spec.get("params", {})),
        paths=dict(ctx_spec.get("paths", {})),
    )
    inputs = {k: load_input(v) for k, v in (spec.get("inputs") or {}).items()}

    log.info("running %s v%s on %s (%s)", module.name, module.version,
             ctx.sample_id, ctx.platform or "no platform")
    result = module.run(ctx, **inputs)

    outdir.mkdir(parents=True, exist_ok=True)
    stem = f"{module.name}_{ctx.sample_id}"
    table_path = outdir / f"{stem}.tsv"
    summary_path = outdir / f"{stem}.summary.json"

    if result.table is not None and hasattr(result.table, "to_csv"):
        result.table.to_csv(table_path, sep="\t", index=False)
    summary = result.summary()
    summary["table"] = str(table_path)
    summary_path.write_text(json.dumps(summary, indent=2, default=str))

    log.info("wrote %s (%d rows) and %s", table_path, result.n_rows, summary_path)
    for note in result.notes:
        log.info("note: %s", note)
    return summary


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--job", required=True, help="job JSON file")
    ap.add_argument("--outdir", default="results/modules")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)-7s %(name)s | %(message)s",
    )
    try:
        summary = run_job(json.loads(Path(args.job).read_text()), Path(args.outdir))
    except ModuleError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(json.dumps({k: summary[k] for k in
                      ("module", "rows", "evidence", "replication_unit", "table")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
