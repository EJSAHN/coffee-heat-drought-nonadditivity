#!/usr/bin/env bash
set -euo pipefail

mkdir -p tools/micromamba logs
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba_roots/coffee_multistress}"
MAMBA="$(pwd)/tools/micromamba/micromamba"
ENV_NAME="coffee-rnaseq"
LOG="logs/pipeline_setup_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

echo "== Project: $(pwd)"
echo "== Log: $LOG"

if [[ ! -x "$MAMBA" ]]; then
  echo "== Installing micromamba locally"
  tmpdir="$(mktemp -d)"
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest -o "$tmpdir/micromamba.tar.bz2"
  tar -xjf "$tmpdir/micromamba.tar.bz2" -C "$tmpdir"
  cp "$tmpdir/bin/micromamba" "$MAMBA"
  chmod +x "$MAMBA"
  rm -rf "$tmpdir"
else
  echo "== micromamba already present: $MAMBA"
fi

if "$MAMBA" env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
  echo "== Environment already exists: $ENV_NAME"
else
  echo "== Creating environment: $ENV_NAME"
  "$MAMBA" create -y -n "$ENV_NAME" \
    -c conda-forge -c bioconda \
    python=3.11 pandas biopython \
    salmon=1.10.3 fastqc multiqc pigz \
    r-base r-readr r-dplyr r-tidyr r-stringr r-tibble r-openxlsx \
    bioconductor-deseq2 bioconductor-tximport r-jsonlite
fi

echo "== Tool versions"
"$MAMBA" run -n "$ENV_NAME" salmon --version || true
"$MAMBA" run -n "$ENV_NAME" python --version
"$MAMBA" run -n "$ENV_NAME" R --version | head -n 1

echo "== DONE setup"

