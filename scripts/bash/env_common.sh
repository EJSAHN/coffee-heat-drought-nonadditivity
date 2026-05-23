#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(pwd)"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba_roots/coffee_multistress}"
MAMBA="$PROJECT_ROOT/tools/micromamba/micromamba"
ENV_NAME="coffee-rnaseq"

run_env() {
  "$MAMBA" run -n "$ENV_NAME" "$@"
}

