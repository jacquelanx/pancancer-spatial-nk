# Module job files

One JSON file per analysis run. The `run_module` Snakemake rule discovers them by path:
`workflow/jobs/{module}/{job}.json` produces `results/modules/{module}/{job}.tsv`.

```json
{
  "module": "niches",
  "context": {
    "sample_id": "pan_cohort",
    "platform": "Xenium",
    "seed": 42,
    "params": {"min_patient_fraction": 0.5, "min_cancer_types": 3}
  },
  "inputs": {
    "sample_compositions": "results/niches/compositions.pkl",
    "reference": "results/niches/reference.npy",
    "reference_names": ["core", "stroma", "immune"],
    "patient_of": "results/manifest/patient_of.json",
    "cancer_type_of": "results/manifest/cancer_type_of.json"
  }
}
```

Inputs are resolved by extension (`.h5ad`, `.tsv`, `.csv`, `.parquet`, `.npy`, `.json`,
`.pkl`); anything that is not an existing path is passed through as a literal, so scalars
and lists can be given inline.

Adding an analysis to the workflow means adding a job file, not a rule. Run one directly:

```bash
python -m panspatial.modules.run --job workflow/jobs/niches/pan_cohort.json
```

`python -m panspatial.modules.cli show <module>` lists the inputs a module expects and
which of its optional dependencies are installed.
