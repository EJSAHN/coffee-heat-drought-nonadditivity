#!/usr/bin/env bash
set -euo pipefail
source scripts/bash/env_common.sh

LOG="logs/pipeline_index_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "== Project: $PROJECT_ROOT"
echo "== Log: $LOG"

mkdir -p data/ref/salmon_transcripts data/ref/salmon_index data/processed/annotations

ARABICA_RNA_GZ="data/ref/arabica_Cara_1_0_NCBI/GCF_003713225.1_Cara_1.0_rna.fna.gz"
CANEPHORA_CDNA_GZ="data/ref/canephora_AUK_PRJEB4211_v1_EnsemblPlants/Coffea_canephora.AUK_PRJEB4211_v1.cdna.all.fa.gz"
ARABICA_FA="data/ref/salmon_transcripts/arabica_Cara_1_0_rna.fa"
CANEPHORA_FA="data/ref/salmon_transcripts/canephora_AUK_PRJEB4211_v1_cdna.fa"

[[ -s "$ARABICA_RNA_GZ" ]] || { echo "Missing $ARABICA_RNA_GZ"; exit 1; }
[[ -s "$CANEPHORA_CDNA_GZ" ]] || { echo "Missing $CANEPHORA_CDNA_GZ"; exit 1; }

if [[ ! -s "$ARABICA_FA" ]]; then
  echo "== Decompressing Arabica transcript FASTA"
  gzip -dc "$ARABICA_RNA_GZ" > "$ARABICA_FA"
fi
if [[ ! -s "$CANEPHORA_FA" ]]; then
  echo "== Decompressing Canephora cDNA FASTA"
  gzip -dc "$CANEPHORA_CDNA_GZ" > "$CANEPHORA_FA"
fi

if [[ ! -s data/processed/annotations/tx2gene_arabica.tsv ]]; then
  echo "== Building tx2gene for Arabica"
  run_env python scripts/python/02_make_tx2gene.py --fasta "$ARABICA_FA" --species arabica --out data/processed/annotations/tx2gene_arabica.tsv
fi
if [[ ! -s data/processed/annotations/tx2gene_canephora.tsv ]]; then
  echo "== Building tx2gene for Canephora"
  run_env python scripts/python/02_make_tx2gene.py --fasta "$CANEPHORA_FA" --species canephora --out data/processed/annotations/tx2gene_canephora.tsv
fi

if [[ ! -d data/ref/salmon_index/arabica_Cara_1_0 ]]; then
  echo "== Salmon index: Arabica"
  run_env salmon index -t "$ARABICA_FA" -i data/ref/salmon_index/arabica_Cara_1_0 -k 31
else
  echo "== Salmon index exists: Arabica"
fi

if [[ ! -d data/ref/salmon_index/canephora_AUK_PRJEB4211_v1 ]]; then
  echo "== Salmon index: Canephora"
  run_env salmon index -t "$CANEPHORA_FA" -i data/ref/salmon_index/canephora_AUK_PRJEB4211_v1 -k 31
else
  echo "== Salmon index exists: Canephora"
fi

echo "== DONE indices and tx2gene"

