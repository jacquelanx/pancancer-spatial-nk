"""Inspect the analysis module registry.

    python -m panspatial.modules.cli list
    python -m panspatial.modules.cli list --phase 5
    python -m panspatial.modules.cli show niches
    python -m panspatial.modules.cli check          # which optional deps are installed
    python -m panspatial.modules.cli rules          # the governing rule per module
"""

from __future__ import annotations

import argparse
import importlib.util
import sys

from panspatial.modules import (
    ModuleError, describe_registry, get_module, modules_for_phase, registered_modules,
)

# The one-line rule each module enforces, shown by `rules`.
GOVERNING_RULES = {
    "qc": "Thresholds are per sample; every cell lost is recorded",
    "integrate": "Refuses global correction when batch is confounded with cancer type",
    "annotate": "Ambiguous and unassigned are real answers, not redistributed",
    "deconv": "Consensus across methods; disagreement is reported, not averaged away",
    "domains": "Method chosen per platform; disagreement becomes boundary-uncertain",
    "niches": "Recurrence tested over patients, across cancer types",
    "malignant": "Caller chosen from allele availability and depth; depth travels with the estimate",
    "nk_spatial": "Single-cell platforms only; torus-shift null preserves autocorrelation",
    "interactions": "Sender and receiver must be in physical contact; co-localisation is not communication",
    "grn": "Per-patient enrichment pooled by meta-analysis; regulons are hypotheses",
    "comodules": "Zsummary < 2 means cohort-specific, never pan-cancer",
    "dynamics": "Two independent method families must agree; velocity alone never suffices",
    "alignment": "Overlap estimated, not assumed; residual QC gates 3D claims",
    "histology": "Imputed evidence: extends cohorts, never validates",
    "predict": "Patient-level folds with a spatial buffer, asserted every fold",
    "translate": "No target without a falsifier; survival must be stage-adjusted",
}


def _installed(dep: str) -> bool:
    root = dep.split()[0].split("[")[0].replace("-", "_")
    aliases = {"scvi_tools": "scvi", "scikit_learn": "sklearn", "scikit_survival": "sksurv",
               "torch_geometric": "torch_geometric", "paste_bio": "paste", "tangram_sc": "tangram",
               "scib_metrics": "scib_metrics", "hdWGCNA_(R)": None}
    name = aliases.get(root, root)
    if name is None:
        return False
    try:
        return importlib.util.find_spec(name) is not None
    except (ImportError, ValueError):
        return False


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="panspatial-modules", description=__doc__)
    sub = p.add_subparsers(dest="command", required=True)

    p_list = sub.add_parser("list", help="list registered analysis modules")
    p_list.add_argument("--phase")

    p_show = sub.add_parser("show", help="describe one module")
    p_show.add_argument("module")

    sub.add_parser("check", help="report which optional dependencies are installed")
    sub.add_parser("rules", help="the governing rule each module enforces")

    args = p.parse_args(argv)

    try:
        if args.command == "list":
            frame = describe_registry()
            if args.phase:
                frame = frame[frame["phase"] == str(args.phase)]
                if frame.empty:
                    print(f"no modules in phase {args.phase}", file=sys.stderr)
                    return 1
            print(f"{len(frame)} module(s)\n")
            for _, r in frame.sort_values(["phase", "module"]).iterrows():
                print(f"  {r['module']:<14} phase {r['phase']}  "
                      f"{r['evidence']:<9} n={r['replication_unit']:<8} [{r['platforms']}]")
                print(f"  {'':<14} {r['description']}")
            return 0

        if args.command == "show":
            m = get_module(args.module)
            print(f"{m.name} v{m.version}  (phase {m.phase})")
            print(f"  {m.description}\n")
            print(f"  evidence          {m.evidence.value}")
            print(f"  replication unit  {m.replication_unit.value}")
            print(f"  platforms         {', '.join(sorted(p.value for p in m.supported_platforms))}")
            print(f"  governing rule    {GOVERNING_RULES.get(m.name, '-')}")
            if m.optional_dependencies:
                print("  dependencies")
                for d in m.optional_dependencies:
                    print(f"    {'OK ' if _installed(d) else '-- '} {d}")
            return 0

        if args.command == "check":
            deps: dict[str, list[str]] = {}
            for name, cls in registered_modules().items():
                for d in cls.optional_dependencies:
                    deps.setdefault(d, []).append(name)
            missing = 0
            for dep in sorted(deps):
                ok = _installed(dep)
                missing += not ok
                print(f"  {'OK ' if ok else 'MISSING'} {dep:<22} used by: {', '.join(deps[dep])}")
            print(f"\n{len(deps) - missing}/{len(deps)} optional dependencies installed.")
            if missing:
                print("Install groups with: pip install 'panspatial[spatial,translate,azure]'")
            return 0

        if args.command == "rules":
            for name in sorted(registered_modules()):
                print(f"  {name:<14} {GOVERNING_RULES.get(name, '-')}")
            return 0

    except ModuleError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
