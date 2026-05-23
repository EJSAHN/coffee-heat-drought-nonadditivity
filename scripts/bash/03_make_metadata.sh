#!/usr/bin/env bash
set -euo pipefail
source scripts/bash/env_common.sh

LOG="logs/pipeline_metadata_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

mkdir -p config data/processed/metadata
run_env python scripts/python/01_make_sample_metadata_auto.py \
  --manifest data/raw/sra_metadata/fastq_download_manifest.tsv \
  --ena data/raw/sra_metadata/coffee_all_ena_runs.tsv \
  --out-auto config/sample_metadata_auto.tsv \
  --out-curated config/sample_metadata_curated.tsv


if [[ -s data/raw/sra_metadata/metadata_triage_rich.csv ]]; then
  run_env python scripts/python/04_curate_sample_metadata_coffee.py \
    --metadata config/sample_metadata_curated.tsv \
    --rich data/raw/sra_metadata/metadata_triage_rich.csv \
    --out config/sample_metadata_curated.tsv
else
  echo "Rich metadata not found; run scripts/powershell/04_download_rich_metadata.ps1, then rerun metadata."
fi

echo "== DONE metadata"

