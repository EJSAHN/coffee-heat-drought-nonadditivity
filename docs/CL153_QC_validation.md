# CL153 QC validation analyses

This repository includes three CL153-focused validation analyses that evaluate whether the primary non-additive gene counts are sensitive to mapping-rate variation or reference framework.

## Analyses

1. **Low-mapping exclusion analysis**  
   Refits the CL153 Salmon/tximport DESeq2 model after excluding the five CL153 libraries with Salmon mapping rates below 50%.

2. **Stage-local interaction analysis**  
   Fits CL153 models using only T25 plus one target stage at a time. The T25+T37 model excludes all T42 and REC14 libraries from dispersion estimation.

3. **STAR genome-alignment validation**  
   Reprocesses all 24 CL153 libraries with STAR against the Ensembl Plants release 62 AUK_PRJEB4211_v1 genome and compares STAR gene-count non-additive totals with the primary Salmon/tximport totals.

## Run

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
.\RUN_PIPELINE.ps1 -Step setup -Threads 8
.\RUN_PIPELINE.ps1 -Step qc-validation -Threads 8
```

The STAR genome-alignment validation is computationally heavier than the table-generation steps. It writes STAR temporary files under the WSL home directory because Windows-mounted drives do not support the FIFO files used by STAR.

## Outputs

```text
results/qc_validation/cl153_low_mapping_exclusion/
results/qc_validation/cl153_stage_local/
results/qc_validation/cl153_star_genome_alignment/
```

The generated Excel workbooks and TSV files are intended for QC reporting and reproducible validation.
