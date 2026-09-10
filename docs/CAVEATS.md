# Constraints that decide whether this study is publishable

Every item here is enforced somewhere in the code, not left to discipline. The enforcement
point is named so it can be audited or deliberately overridden with a recorded reason.

---

## 1. NK cells are below the resolution of spot-based platforms

NK cells are 0.5–3% of cells in most solid tumors. A 55 µm Visium spot contains 1–10 cells.
An "NK cell" on Visium is therefore a deconvolution weight — a similarity score between the
spot's expression and a reference NK profile — and the nearest competing reference profile
is CD8 effector T cells, whose cytotoxic program is largely shared. At realistic sequencing
depth those two weights are not separable, and the resulting "NK spatial map" reproduces
the CD8 map.

**Consequence.** A spatial claim resolved to NK cells requires Xenium, CosMx, or MERFISH.
Visium and Slide-seq contribute niche context — stromal architecture, TLS position, CAF
density — not NK-resolved conclusions.

**Enforced by** `panspatial.nk.require_single_cell_resolution`, which raises on spot
platforms, and by the `ST_is_single_cell` column that
`panspatial.ingest.geo_sra` writes into the manifest.

**Also required.** `R/harmonize_deconvolve.R::benchmark_deconvolution` builds simulated spot
mixtures of known composition and reports per-cell-type Pearson r. Cell types below the
threshold are flagged unreliable and barred from carrying a conclusion. Run it; do not skip
it because the pipeline is slow.

---

## 2. The replication unit is the patient

A pan-cancer atlas has ~10⁶ cells and ~10² patients. A test whose n is cells will return
p < 10⁻³⁰⁰ for effects that do not replicate in a second cohort, because cells within a
patient are not independent draws — they share genotype, treatment history, batch, and
section.

**Consequence.** Compute the statistic per section, pool sections to patients, then
meta-analyze across patients. The n you may quote is the number of patients.

**Enforced by** `panspatial.nk.per_section_niche_test` (one row per section, carrying its
patient) and `panspatial.nk.cohort_meta_analysis`, which averages within patient before
pooling and refuses fewer than two units.

**Report alongside every pooled effect:** I², and the leave-one-out range. An effect whose
LOO range collapses when one cohort is dropped is a single-cohort finding reported with a
pan-cancer n. `MetaAnalysisResult.is_driven_by_one_cohort` flags this.

---

## 3. Spatial nulls must preserve autocorrelation

Cell-type labels in tissue are strongly spatially autocorrelated. Shuffling labels over
positions destroys that structure, producing a null distribution far tighter than reality.

Measured on simulated tissue (two independently clustered patterns, no true association,
α = 0.05, 24 replicates — `tests/test_spatial_stats.py`):

| Null model | False-positive rate |
|---|---|
| Torus shift (preserves autocorrelation) | 2/24 (8%) |
| Label shuffle over pooled positions | 24/24 (100%) |

The naive test rejects on every single true-null replicate.

**Enforced by** `panspatial.stats.torus_shift_permutation_test`, which rigidly translates
the focal pattern on the section's torus so only its registration against the anchor is
randomized. Pass the section's tissue mask via `inside_tissue`: without it, shifted points
land in empty slide area, the null inflates, and the test manufactures "attraction".

For thin or non-convex sections the torus shift breaks down; the function raises rather
than returning a quietly wrong answer, and a mask-restricted block bootstrap is the
alternative.

---

## 4. Batch is confounded with cancer type

In a pan-cancer cohort assembled from public data, "batch" and "cancer type" are the same
variable. Integration that removes batch removes the biology.

**Consequence.** Correct within cancer type over patient. Compare across cancer types on
shared cell-type labels and per-sample summary statistics — replication across cohorts —
not in a shared latent space.

**Enforced by** `R/harmonize_deconvolve.R::integrate_within_cancer_type`.

---

## 5. Multiple testing is corrected across everything tested

A screen over cell-type pairs × cancer types × anchors that corrects each cancer type
separately and reports the union has not corrected at all.

**Enforced by** `panspatial.stats.benjamini_hochberg(..., family_size=...)`, which requires
the family size explicitly and refuses a family smaller than the vector passed, and by the
`--family-size` argument threaded through `panspatial.nk.run_meta` from the Snakemake rule.

---

## 6. Spatial ML splits leak by default

Two spots 100 µm apart share microenvironment, patient, and batch. Split over spots and the
model memorizes; AUC 0.95 means nothing.

**Consequence.** Split on patients. For a deliberate region-held-out design within one
section, add a spatial buffer of several times the platform's spot pitch.

**Enforced by** `panspatial.stats.patient_blocked_splits`, `spatial_buffer_filter`, and
`assert_no_leakage` — the last of which must be called inside the CV loop, not once. Its
`allow_shared_patients` escape hatch requires a positive buffer, because relaxing the
patient check with no buffer leaves nothing guarding the split.

Feature selection, HVG choice, scaling, and integration must be fitted inside the training
fold. Fitting them on the full cohort before splitting leaks the test set into the model.

---

## 7. NK identity is a discrimination problem

`FCGR3A` is on CD16⁺ monocytes. `NCAM1` is on neurons and some tumors. NK and CD8 effector
T transcriptomes are close, and the separating evidence — absence of `CD3D/E/G`, `TRAC` —
is a dropout-prone negative at low depth. `KLRF1` and `NCR1` carry the positive
discrimination.

NKT cells, γδ T cells, and ILC1 sit in the gap. They are reported as ambiguous, not
absorbed into "NK".

**Enforced by** `panspatial.nk.nk_identity_report`, which counts what was dropped and why,
warns when neither `KLRF1` nor `NCR1` is in the panel, and writes the per-sample NK fraction
that belongs in the supplement.

**States, not clusters.** `score_nk_states` leaves cells whose top two program scores are
within `min_margin` as `unassigned`, and reports bootstrap ARI. A state count that is not
stable to resampling is a continuum, and should be described as one.

---

## 8. Co-localization is not communication

A ligand–receptor score is co-expression in adjacent neighborhoods. It is a hypothesis about
signalling. `HLA-E`–`KLRC1` proximity is not evidence that the axis is active, and neither
is a significant `squidpy.gr.ligrec` result.

The verb must match the evidence: *associated with* for observational spatial correlation;
*drives* or *mediates* only where a perturbation experiment exists.

**Enforced by** the system anchor's evidence-discipline block and by the phase 5 prompts,
which require every claim to be traced to a line in the verified-results table.

(Note: `squidpy.gr.ligand_receptor_score` does not exist. The function is
`squidpy.gr.ligrec`.)

---

## 9. Manuscript text may only quote executed results

The failure mode that ends a project: an LLM asked to draft an abstract "based on our
finding that X" produces a fluent paragraph containing a plausible hazard ratio that was
never computed, and it survives into a submission.

**Enforced by** three layers:

1. Templates that write result-bearing prose set `requires_verified_results: true` and
   refuse to render without `--results-verified`.
2. `results/report/verified_results.tsv`, built by `panspatial.report.build` from executed
   outputs, is the only permitted source of numbers. Anything absent is written `[TBD]`.
3. The system anchor forbids inventing values and requires an `INSUFFICIENT EVIDENCE`
   section when the input cannot support the requested claim.

Every call writes a run record to `runs/` with the prompt fingerprint, the bound inputs, the
provider, the resolved model or deployment, and the response, so any sentence in the draft
can be traced back to what produced it.

All three layers sit above the provider backend, so they apply identically on Anthropic and
Azure OpenAI. Switching providers changes which model wrote the text; it does not change
what the text is allowed to claim. Record the provider and model in the methods — they are
in the run record for exactly this reason.
