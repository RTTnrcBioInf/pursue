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
  # Keep the user's existing library on the path: R_LIBS_USER is a colon-separated list, and
  # replacing it outright hides packages already installed on this cluster (that is what hid
  # MicrobiomeStat / corncob / LDM / LOCOM in the 2026-09-11 smoke test). New installs go to
  # the first entry; everything already installed stays visible.
  BENCH_LIB="$HOME/R/pursue-bench-lib"; mkdir -p "$BENCH_LIB"
  BASE_LIB="$(R --no-echo -e 'cat(Sys.getenv("R_LIBS_USER"))' 2>/dev/null | tr ':' '\n' | grep -v "pursue-bench-lib" | paste -sd: -)"
  export R_LIBS_USER="$BENCH_LIB${BASE_LIB:+:$BASE_LIB}"
  echo ">> R library path: $R_LIBS_USER (export R_LIBS_USER before running jobs)"
  Rscript "$HERE/install_packages.R"
fi
echo ">> now: bash hpc/download_templates.sh && Rscript hpc/smoke_test.R"
