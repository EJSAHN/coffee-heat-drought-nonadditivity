#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(pwd)"
source scripts/bash/env_common.sh

THREADS="${THREADS:-8}"
MAMBA="${MAMBA:-./tools/micromamba/micromamba}"
ENV_NAME="${ENV_NAME:-coffee-rnaseq}"

mkdir -p results/qc results/multiqc results/salmon logs

echo "== Resume-aware FastQC + Salmon"
echo "== Project: $PROJECT_ROOT"
echo "== Threads: $THREADS"


META="config/sample_metadata_curated.tsv"
if [[ ! -s "$META" ]]; then
  META="config/sample_metadata_auto.tsv"
fi
[[ -s "$META" ]] || { echo "Missing sample metadata. Run metadata step first."; exit 1; }

run_env python scripts/python/03_write_salmon_commands.py   --metadata "$META"   --threads "$THREADS"   --out logs/salmon_commands.sh

echo "== FastQC resume mode"
mapfile -t FQS < <(find data/raw/fastq -type f -name "*.fastq.gz" | sort)
echo "FASTQ files found for QC: ${#FQS[@]}"

for fq in "${FQS[@]}"; do
    base="$(basename "$fq" .fastq.gz)"
    fq_dir="$(dirname "$fq")"

    existing="$(find results/qc data/raw/fastq -type f \( -name "${base}*fastqc.zip" -o -name "${base}*fastqc.html" \) -size +0c 2>/dev/null | head -n 1 || true)"

    if [[ -n "$existing" ]]; then
        echo "FastQC skip: $fq  [found: $existing]"
        continue
    fi

    echo "FastQC run:  $fq"
    "$MAMBA" run -n "$ENV_NAME" fastqc -t 1 -o results/qc "$fq"
done

echo "== MultiQC after FastQC"
"$MAMBA" run -n "$ENV_NAME" multiqc -f -o results/multiqc results/qc || true

if [[ ! -s logs/salmon_commands.sh ]]; then
    echo "ERROR: logs/salmon_commands.sh not found."
    echo "Run original Quant once until it writes logs/salmon_commands.sh."
    exit 1
fi

echo "== Salmon quant resume mode"
done_count=0
run_count=0

while IFS= read -r cmd; do
    [[ -z "$cmd" ]] && continue
    [[ "$cmd" =~ ^# ]] && continue

    out=""
    if [[ "$cmd" =~ [[:space:]]-o[[:space:]]+([^[:space:]]+) ]]; then
        out="${BASH_REMATCH[1]}"
    elif [[ "$cmd" =~ [[:space:]]--output[[:space:]]+([^[:space:]]+) ]]; then
        out="${BASH_REMATCH[1]}"
    fi

    if [[ -n "$out" && -s "$out/quant.sf" ]]; then
        echo "Salmon skip: $out"
        done_count=$((done_count + 1))
    else
        echo "Salmon run:  ${out:-unknown_output}"
        "$MAMBA" run -n "$ENV_NAME" bash -lc "$cmd"
        run_count=$((run_count + 1))
    fi
done < logs/salmon_commands.sh

echo "== Salmon summary: skipped=$done_count newly_run=$run_count"

echo "== Final MultiQC"
"$MAMBA" run -n "$ENV_NAME" multiqc -f -o results/multiqc results/qc results/salmon || true

echo "== DONE resume-aware Quant"