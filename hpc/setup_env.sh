#!/usr/bin/env bash
# Create the benchmark environment on the HPC. Two routes, pick whichever the cluster
# supports; both end by running install_packages.R.
#
#   A) conda/mamba (preferred, self-contained):   bash hpc/setup_env.sh conda
#   B) module-provided R (>= 4.3) + user library:  bash hpc/setup_env.sh module [R-module-name]
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(cd "$HERE/.." && pwd)"
MODE="${1:-conda}"
if [ "$MODE" = "conda" ]; then
  if command -v mamba >/dev/null; then CONDA=mamba; elif command -v conda >/dev/null; then CONDA=conda; else echo "no conda/mamba found; use: bash hpc/setup_env.sh module"; exit 1; fi
  $CONDA env create -f "$HERE/environment.yml" -n pursue-bench || $CONDA env update -f "$HERE/environment.yml" -n pursue-bench
  echo ">> activate with:  conda activate pursue-bench"
  eval "$($CONDA shell.bash hook)"; conda activate pursue-bench
  Rscript "$HERE/install_packages.R"
else
  MOD="${2:-R}"; module load "$MOD" || true
  export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/pursue-bench-lib}"; mkdir -p "$R_LIBS_USER"
  echo ">> R user library: $R_LIBS_USER (export R_LIBS_USER before running jobs)"
  Rscript "$HERE/install_packages.R"
fi
echo ">> now: bash hpc/download_templates.sh && Rscript hpc/smoke_test.R"
