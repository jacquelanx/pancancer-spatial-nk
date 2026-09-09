# Pan-cancer spatial + single-cell NK atlas — LLM-orchestrated pipeline

An implementation of the six-phase prompt architecture: a versioned prompt registry driving
Claude, and the analysis code the prompts describe, wired together by Snakemake.

The organising idea is that the prompts are **code artifacts**, not chat messages. Each one
is a file with frontmatter, a version, declared variables, and a provenance fingerprint.
Nothing is sent to a model that was not rendered from a template in `prompts/`.

```
prompts/          versioned templates; system/00_meta.md is the anchor on every call
src/panspatial/
  orchestrator/   registry (strict binding, gating, fingerprints), client, CLI
  ingest/         GEO/SRA harvesting and platform classification    [phase 2]
  nk/             NK identity, states, spatial niche testing        [phase 3]
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

Run it for real (needs `ANTHROPIC_API_KEY` or `ant auth login`):

```bash
python -m panspatial.orchestrator.cli run phase6.audit \
    --var target_paths=src/panspatial/nk/nk_spatial.py \
    --var-file source_bundle=src/panspatial/nk/nk_spatial.py \
    --out results/audit/nk_spatial.md
```

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

Executed and tested in this repository: the orchestrator, ingest classification, NK spatial
layer, and statistics — 93 tests, all passing, including simulated-tissue end-to-end runs
through meta-analysis and report assembly.

Not executed here: `R/harmonize_deconvolve.R` (no R toolchain in this environment) and the
AnnData-backed entrypoints (`run_extract`, `run_states`, `run_niche`), which need
scanpy/squidpy and real data. They are written against the documented APIs; run the
Snakemake dry run and a single-sample smoke test before committing a cohort to them.
