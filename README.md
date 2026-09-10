# Pan-cancer spatial + single-cell NK atlas — LLM-orchestrated pipeline

An implementation of the six-phase prompt architecture: a versioned prompt registry driving
an LLM (Anthropic or Azure OpenAI), and the analysis code the prompts describe, wired
together by Snakemake.

The organising idea is that the prompts are **code artifacts**, not chat messages. Each one
is a file with frontmatter, a version, declared variables, and a provenance fingerprint.
Nothing is sent to a model that was not rendered from a template in `prompts/`.

```
prompts/          versioned templates; system/00_meta.md is the anchor on every call
src/panspatial/
  orchestrator/   registry (strict binding, gating, fingerprints), backends, client, CLI
  modules/        the analysis-module contract, registry and generic runner
  ingest/         GEO/SRA harvesting and platform classification    [phase 2]
  qc/             per-sample adaptive QC, coordinate normalisation  [phase 2]
  integrate/      confounding-aware batch integration               [phase 2]
  annotate/       cell typing + the foundation-model admission gate [phase 2]
  deconv/         consensus deconvolution + reliability gate        [phase 4]
  domains/        platform-aware domain detection + consensus       [phase 4]
  niches/         niche matching and patient-level recurrence       [phase 4]
  malignant/      CNV caller selection, cNMF meta-programs          [phase 4]
  nk/             NK identity, states, spatial niche testing        [phase 4]
  interactions/   spatially constrained CCI + rank consensus        [phase 5]
  grn/            regulon niche enrichment and conservation         [phase 5]
  comodules/      co-expression modules + cross-cohort preservation [phase 5]
  dynamics/       fate mapping with multi-view corroboration        [phase 5]
  alignment/      multi-section alignment + residual QC             [phase 5]
  histology/      H&E expression prediction, flagged as imputation  [phase 5]
  predict/        outcome models under leakage-safe CV              [phase 6]
  translate/      signatures, adjusted survival, target ranking     [phase 6]
  stats/          spatial nulls, meta-analysis, CV leakage guards   [cross-cutting]
  report/         the verified-results table that gates drafting    [phase 5]
R/                Seurat v5 + Harmony + RCTD with a reliability benchmark [phase 2]
workflow/         Snakefile, config, per-rule resources             [phase 6]
docs/CAVEATS.md   the constraints that decide whether this is publishable
runs/             one JSON record per LLM call: prompt, inputs, usage, response
```

## Quick start

```bash
pip install -e ".[dev]"
PYTHONPATH=src python -m pytest tests/ -q
```

List and inspect the prompt library:

```bash
python -m panspatial.orchestrator.cli list
python -m panspatial.orchestrator.cli show phase3.nk_spatial
```

Render a prompt without calling the API — useful for reviewing what a phase will actually
ask for:

```bash
python -m panspatial.orchestrator.cli run phase3.nk_spatial \
    --var adata_path=data/pan_cancer.h5ad --var platform=Xenium --dry-run
```

Run it for real:

```bash
python -m panspatial.orchestrator.cli run phase6.audit \
    --var target_paths=src/panspatial/nk/nk_spatial.py \
    --var-file source_bundle=src/panspatial/nk/nk_spatial.py \
    --out results/audit/nk_spatial.md
```

## The analysis modules

Sixteen analysis layers, each implementing a common contract
([modules/base.py](src/panspatial/modules/base.py)). A module cannot produce a result
without declaring four things, and the framework checks them:

- **Evidence level** — `measured`, `inferred` or `imputed`. `assert_can_validate` refuses a
  validation whose evidence is weaker than the claim, so histology-predicted expression
  cannot validate the spatial model it came from.
- **Replication unit** — anything at cell or spot level is marked as such, so the
  aggregation layer knows it must pool before a cross-patient claim.
- **Platform support** — a rare-cell claim on a 55 µm spot is refused at the point of
  production, not caught in review.
- **Provenance** — parameters, seed, and external tool versions, hashed into a fingerprint.

```bash
python -m panspatial.modules.cli list          # every layer, phase, evidence class
python -m panspatial.modules.cli rules         # the one rule each layer enforces
python -m panspatial.modules.cli show niches   # detail + dependency availability
python -m panspatial.modules.cli check         # which optional deps are installed
```

Any module runs through one generic entrypoint driven by a job file, so adding a layer to
the workflow adds a job file rather than a rule:

```bash
python -m panspatial.modules.run --job workflow/jobs/niches/pan_cohort.json
```

| Layer | Governing rule enforced in code |
|---|---|
| `qc` | Thresholds are per sample; every cell lost is recorded |
| `integrate` | Refuses global correction when batch is confounded with cancer type |
| `annotate` | Ambiguous and unassigned are real answers, not redistributed |
| `deconv` | Consensus across methods; disagreement reported, not averaged away |
| `domains` | Method chosen per platform; disagreement becomes boundary-uncertain |
| `niches` | Recurrence tested over patients, across cancer types |
| `malignant` | Caller chosen from allele availability and depth; depth travels with the estimate |
| `nk_spatial` | Single-cell platforms only; torus-shift null preserves autocorrelation |
| `interactions` | Sender and receiver must be in contact; co-localisation is not communication |
| `grn` | Per-patient enrichment pooled by meta-analysis; regulons are hypotheses |
| `comodules` | Zsummary < 2 means cohort-specific, never pan-cancer |
| `dynamics` | Two independent method families must agree; velocity alone never suffices |
| `alignment` | Overlap estimated, not assumed; residual QC gates 3D claims |
| `histology` | Imputed evidence: extends cohorts, never validates |
| `predict` | Patient-level folds with a spatial buffer, asserted every fold |
| `translate` | No target without a falsifier; survival must be stage-adjusted |

Method choices are tied to published benchmarks, including where those benchmarks disagree
— the deconvolution layer runs four tools because no two benchmarks agree on a single
winner, and the annotation layer refuses a foundation model that cannot beat scVI, Harmony
and a plain HVG baseline on held-out patients.

## Providers

Two backends, selected by `--backend {anthropic,azure}`, `PANSPATIAL_BACKEND`, or
auto-detection from the environment. Copy `.env.example` to `.env` and fill it in; the CLI
loads it automatically and real environment variables always win over the file.

```bash
cp .env.example .env          # then fill in
pip install -e ".[azure]"     # or ".[llm]" for Anthropic
python -m panspatial.orchestrator.cli list        # prints the active backend
```

**Anthropic** needs `ANTHROPIC_API_KEY` or `ant auth login`.

**Azure OpenAI** needs `AZURE_OPENAI_API_KEY`, `AZURE_OPENAI_ENDPOINT`,
`AZURE_OPENAI_API_VERSION`, and `AZURE_OPENAI_DEPLOYMENTS`. Three Azure-specific things the
backend handles so callers never see them:

*Deployment routing.* Azure routes on deployment name, not model name. Prompt templates name
a logical model (`claude-opus-5`); `AZURE_OPENAI_DEPLOYMENTS` maps it to whatever your
deployment is called, so the same templates run on either provider unedited:

```
AZURE_OPENAI_DEPLOYMENTS=claude-opus-5=my-gpt5-deployment,fast=my-4o-mini-deployment
```

An unmapped model falls back to the first entry. A bare name (`AZURE_OPENAI_DEPLOYMENTS=my-deployment`)
maps to itself.

*Reasoning vs chat parameters.* Reasoning deployments want `max_completion_tokens` and accept
`reasoning_effort`; chat deployments want `max_tokens` and reject it. Which is which cannot
be inferred from an arbitrary deployment name, so capabilities are probed once per deployment
from the API's own 400 response and cached for the process. Set
`AZURE_OPENAI_REASONING_DEPLOYMENTS` / `AZURE_OPENAI_CHAT_DEPLOYMENTS` to skip the probe. The
five-level `effort` in templates maps onto Azure's three (`xhigh`/`max` → `high`).

*Caching and cost.* Azure prompt caching is automatic and prefix-based above ~1024 tokens —
there is no `cache_control`, so the frozen system anchor is what earns the hit. Azure pricing
varies by region, tier, and agreement, so run records report tokens and leave cost
`unpriced` rather than printing an Anthropic list price that would be wrong.

Both backends normalize the two failure modes that must never pass silently: truncation
(`finish_reason == "length"` / `stop_reason == "max_tokens"`) and policy refusal (Azure
content filter / Anthropic refusal). Both raise; neither returns partial text that reads
complete.

## What the orchestrator layer guarantees

**Strict variable binding.** An unbound `{{variable}}` raises instead of being sent
literally. Unresolved placeholders in a manuscript prompt are how invented numbers get into
a draft.

**A frozen, cached system anchor.** `prompts/system/00_meta.md` is byte-identical on every
call and carries the statistical constraints, so they are always in force — and, being the
cached prefix, nearly free. The client warns if a call shows no cache activity, which means
something per-run has leaked into `system`.

**Provenance.** Every call writes `runs/<timestamp>_<template>_<fingerprint>.json` with the
prompt, the bound inputs, token usage, and the response. The fingerprint is a hash of
(system + rendered prompt + model + effort).

**Gating on executed results.** Templates that write result-bearing prose carry
`requires_verified_results: true` and refuse to render without `--results-verified`. The
permitted numbers live in `results/report/verified_results.tsv`, built from executed
outputs. See [docs/CAVEATS.md §9](docs/CAVEATS.md).

**Truncation is an error.** A response that stops at `max_tokens` raises rather than
returning a half-finished audit that reads complete.

All four guarantees are provider-independent — they live above the backend layer, so
switching to Azure does not weaken the manuscript gate.

## What the analysis layer enforces

These are the parts that determine whether the study survives review. Each is enforced in
code and documented in [docs/CAVEATS.md](docs/CAVEATS.md):

| Constraint | Enforcement |
|---|---|
| NK-resolved spatial claims need single-cell-resolution platforms | `nk.require_single_cell_resolution` |
| The replication unit is the patient, not the cell or spot | `nk.per_section_niche_test` → `nk.cohort_meta_analysis` |
| Spatial nulls preserve autocorrelation | `stats.torus_shift_permutation_test` |
| FDR is corrected across the whole family tested | `stats.benjamini_hochberg(family_size=...)` |
| CV splits are patient-level with a spatial buffer | `stats.patient_blocked_splits`, `assert_no_leakage` |
| Deconvolution weights are benchmarked before use | `R/harmonize_deconvolve.R::benchmark_deconvolution` |
| Ambiguous NK calls are reported, not absorbed | `nk.nk_identity_report` |

On simulated tissue with no true association, the torus-shift null rejects at 8% (2/24) at
α = 0.05 while a naive label shuffle rejects at 100% (24/24). That difference is the
difference between a finding and an artifact; the test lives in
`tests/test_spatial_stats.py`.

## Running the workflow

```bash
snakemake -s workflow/Snakefile --use-conda --cores 16 -n     # dry run
snakemake -s workflow/Snakefile --use-conda --profile slurm   # cluster
```

`workflow/config/config.yaml` holds the seed, QC and spatial parameters, and the cohort and
sample tables. Populate `samples:` from the harvested manifest plus curation; the
`platform_st` field decides whether a sample takes the NK-resolved path or the
niche-context path.

## Status

Executed and tested in this repository: the orchestrator (both backends), the module
contract and all sixteen analysis layers' decision logic, ingest classification, the NK
spatial layer, and the statistics — **233 tests, all passing**, including simulated-tissue
runs end to end through meta-analysis, report assembly, and the generic module runner.

What is tested is the logic this repository owns: the gates, the consensus and
concordance calculations, the reliability and preservation statistics, the recurrence and
corroboration rules, and the leakage assertions. External heavy tools (cell2location, RCTD,
GraphST, NicheCompass, pySCENIC, CellRank, PASTE) are called through adapters that are
import-guarded and unexercised here — install the relevant extra and smoke-test one sample
before committing a cohort.

Not exercised against a live Azure endpoint — no credentials here. The first real call will
confirm the deployment mapping and capability probe against your actual deployments; run one
`--dry-run --show-target` first to check routing without spending tokens.

Not executed here: `R/harmonize_deconvolve.R` (no R toolchain in this environment) and the
AnnData-backed entrypoints (`run_extract`, `run_states`, `run_niche`), which need
scanpy/squidpy and real data. They are written against the documented APIs; run the
Snakemake dry run and a single-sample smoke test before committing a cohort to them.
