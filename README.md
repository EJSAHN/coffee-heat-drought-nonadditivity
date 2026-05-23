# Coffee heat-drought non-additivity reanalysis

## Repository description

Reproducible Windows/WSL workflow for reconstructing public *Coffea arabica* Icatu and *Coffea canephora* CL153 RNA-seq data, quantifying transcript abundance with Salmon, modeling heat-drought non-additivity with DESeq2, and generating analysis-ready tables.

## Recommended repository name

`coffee-heat-drought-nonadditivity`

## Short GitHub description

Reproducible reanalysis pipeline for coffee heat-drought RNA-seq non-additivity, recovery memory, functional modules, and analysis-ready tables.

## Overview

This workflow reconstructs a balanced public RNA-seq design from five BioProjects:

- PRJNA787748
- PRJNA1087442
- PRJNA1087679
- PRJNA1088119
- PRJNA1135679

The final curated design contains:

- *C. arabica* cv. Icatu and *C. canephora* cv. CL153
- WW and SWD water regimes
- T25, T37, T42, and REC14 stages
- 3 biological replicates per genotype-water-stage combination
- 48 retained RNA-seq libraries

The workflow writes all generated files under `results/`.

## Directory layout

```text
scripts/
  powershell/   Windows and WSL runners
  bash/         Linux workflow steps executed inside WSL
  python/       Metadata parsing, metadata curation, tx2gene, Salmon command writing
  R/            DESeq2, annotation, QC tables, candidate prioritization, table synthesis
config/         Curated sample metadata after pipeline execution
manuscript/     Methods and reference draft text
results/        Generated outputs; not tracked by git
```

## Quick start

Run commands from the cloned repository root. The PowerShell scripts automatically use the repository root by default; pass `-ProjectRoot` only if you want to run from another location.

### 1. Bootstrap project metadata, references, and optional FASTQ download

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
.\scripts\powershell\00_bootstrap_coffee_multistress.ps1
.\scripts\powershell\04_download_rich_metadata.ps1
```

Optional custom project root example:

```powershell
.\RUN_PIPELINE.ps1 -Step all -ProjectRoot "E:\coffee-heat-drought-nonadditivity"
```

Test FASTQ download:

```powershell
.\scripts\powershell\02_download_fastq.ps1 -Mode test -MaxRuns 2
```

Full FASTQ download:

```powershell
.\scripts\powershell\02_download_fastq.ps1 -Mode all
```

### 2. Setup the WSL environment

```powershell
.\RUN_PIPELINE.ps1 -Step setup -Threads 8
```

### 3. Build transcriptome indices and tx2gene tables

```powershell
.\RUN_PIPELINE.ps1 -Step index -Threads 8
```

### 4. Generate and curate sample metadata

```powershell
.\RUN_PIPELINE.ps1 -Step metadata -Threads 8
```

This step creates `config/sample_metadata_auto.tsv` and `config/sample_metadata_curated.tsv`. If `data/raw/sra_metadata/metadata_triage_rich.csv` is available, the coffee-specific curation script retains the 48 target samples and excludes unrelated runs.

### 5. Run FastQC, MultiQC, and Salmon quantification

```powershell
.\RUN_PIPELINE.ps1 -Step quant -Threads 8
```

The quantification step is resume-aware. Existing FastQC outputs and existing Salmon `quant.sf` files are skipped.

### 6. Run DESeq2 non-additivity modeling

```powershell
.\RUN_PIPELINE.ps1 -Step deseq2 -Threads 8
```

### 7. Generate annotation/module summaries, QC sensitivity tables, and candidate synthesis

```powershell
.\RUN_PIPELINE.ps1 -Step next -Threads 8
```

### 8. Build synthesized analysis tables

The `next` step writes annotation, module, QC, candidate-prioritization tables, and a consolidated Excel workbook. The `all` step runs through this table-synthesis stage. This repository intentionally does not include manuscript figure-generation code.

## Main outputs

```text
results/deseq2/<analysis_group>/
results/tables/*_nonadditivity_scores.tsv
results/tables/*_top50_interaction_candidates.tsv
results/next/tables/
results/paper_pack_v1/
results/paper_pack_v1/coffee_analysis_tables.xlsx
```

## Manuscript methods

A draft Materials and Methods section and the corresponding reference list are provided in:

```text
manuscript/materials_and_methods.md
manuscript/references_cited.md
```

## Notes

- Raw FASTQ files, reference files, and generated results are intentionally excluded from git.
- PowerShell scripts default to the cloned repository root and also accept `-ProjectRoot` for custom locations.
- The micromamba root defaults to `$HOME/micromamba_roots/coffee_multistress` inside WSL.
