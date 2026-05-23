#!/usr/bin/env bash
set -euo pipefail
source scripts/bash/env_common.sh

LOG="logs/pipeline_deseq2_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

mkdir -p results/deseq2 results/tables
run_env Rscript scripts/R/04_deseq2_nonadditive.R

echo "== DONE DESeq2 non-additivity"

