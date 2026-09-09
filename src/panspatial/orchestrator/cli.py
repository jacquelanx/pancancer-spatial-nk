"""CLI for the prompt pipeline.

    python -m panspatial.orchestrator.cli list
    python -m panspatial.orchestrator.cli show phase3.nk_spatial
    python -m panspatial.orchestrator.cli run phase3.nk_spatial \
        --var adata_path=data/pan_cancer.h5ad --var platform=Xenium --dry-run
    python -m panspatial.orchestrator.cli run phase6.audit \
        --var target_paths=src/panspatial/nk/nk_spatial.py \
        --var-file source_bundle=src/panspatial/nk/nk_spatial.py --out audit.md

``--var-file`` reads a file into a variable, which is how source bundles and results tables
get into a prompt without shell quoting hazards.
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from panspatial.orchestrator.client import LLMError, OrchestratorClient
from panspatial.orchestrator.registry import PromptError, PromptRegistry


def _collect_vars(args: argparse.Namespace) -> dict[str, str]:
    values: dict[str, str] = {}
    for item in args.var or []:
        key, sep, val = item.partition("=")
        if not sep:
            raise SystemExit(f"--var expects key=value, got {item!r}")
        values[key.strip()] = val
    for item in args.var_file or []:
        key, sep, path = item.partition("=")
        if not sep:
            raise SystemExit(f"--var-file expects key=path, got {item!r}")
        target = Path(path)
        if not target.is_file():
            raise SystemExit(f"--var-file {key}: no such file: {path}")
        values[key.strip()] = target.read_text(encoding="utf-8")
    return values


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="panspatial-prompt", description=__doc__)
    parser.add_argument("--prompts", default="prompts", help="prompt root (default: prompts)")
    parser.add_argument("--runs", default="runs", help="run-record directory (default: runs)")
    parser.add_argument("-v", "--verbose", action="store_true")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("list", help="list registered prompts")

    p_show = sub.add_parser("show", help="print a template and its variables")
    p_show.add_argument("template_id")

    p_run = sub.add_parser("run", help="render and execute a prompt")
    p_run.add_argument("template_id")
    p_run.add_argument("--var", action="append", metavar="KEY=VALUE")
    p_run.add_argument("--var-file", action="append", metavar="KEY=PATH")
    p_run.add_argument("--model", help="override the template's model")
    p_run.add_argument("--effort", choices=["low", "medium", "high", "xhigh", "max"])
    p_run.add_argument("--schema", help="JSON Schema file for structured output")
    p_run.add_argument("--out", help="write the response text here")
    p_run.add_argument("--dry-run", action="store_true", help="render only; no API call")
    p_run.add_argument(
        "--results-verified",
        action="store_true",
        help="assert that result values came from an executed analysis (required by "
        "manuscript templates)",
    )

    args = parser.parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)-7s %(name)s | %(message)s",
    )

    try:
        registry = PromptRegistry(args.prompts)

        if args.command == "list":
            for tid in registry.ids():
                t = registry.get(tid)
                gate = " [gated: requires verified results]" if t.requires_verified_results else ""
                print(f"{tid:26s} v{t.version}  phase {t.phase:3s} {t.model} ({t.effort}){gate}")
                print(f"{'':26s} {t.meta.get('description', '')}")
            return 0

        if args.command == "show":
            t = registry.get(args.template_id)
            print(f"# {t.id} v{t.version}  ({t.path})")
            for name, spec in t.variables.items():
                flag = "required" if spec.required else f"default={spec.default!r}"
                print(f"#   {name:20s} [{flag}] {spec.description}")
            print()
            print(t.body)
            return 0

        values = _collect_vars(args)
        prompt = registry.render(
            args.template_id,
            values,
            results_verified=args.results_verified,
            model=args.model,
            effort=args.effort,
        )
        schema = json.loads(Path(args.schema).read_text()) if args.schema else None
        if schema is None and prompt.structured_output:
            candidate = Path(args.prompts) / prompt.structured_output
            if candidate.is_file():
                schema = json.loads(candidate.read_text())

        client = OrchestratorClient(args.runs, dry_run=args.dry_run)
        result = client.run(prompt, schema=schema, write_record=not args.dry_run)

        if args.dry_run:
            print(f"--- system anchor ({len(prompt.system)} chars, cached prefix) ---")
            print(f"--- {prompt.template_id} v{prompt.template_version} fp={prompt.fingerprint} ---")
            print(prompt.user)
            return 0

        if args.out:
            Path(args.out).write_text(result.text, encoding="utf-8")
            print(f"wrote {args.out} ({len(result.text)} chars, ${result.cost_usd:.4f})")
        else:
            print(result.text)
        return 0

    except (PromptError, LLMError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
