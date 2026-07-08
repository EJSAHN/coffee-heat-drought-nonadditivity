# Results directory

Generated outputs are written here by the analysis workflow and are excluded from git.

Primary outputs include:

```text
results/deseq2/
results/tables/
results/next/
results/analysis_tables/
results/qc_validation/
```

The `qc_validation` directory is produced only when `RUN_PIPELINE.ps1 -Step qc-validation` is executed.
