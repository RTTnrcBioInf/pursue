#!/usr/bin/env bash
# Create the benchmark environment on the HPC, then install every R package.
#
#   bash hpc/setup_env.sh                  # conda-family env named pursue-bench (recommended)
#   bash hpc/setup_env.sh conda --fresh    # delete and rebuild the env first
#   bash hpc/setup_env.sh module R/4.5.0   # a cluster module instead, if one exists
#
# Solver preference is micromamba > mamba > conda. That order is deliberate: plain conda
# crashed on 2026-09-11 linking a bioconda package with "ValueError: unsupported format
# character 'T'", a conda bug in its Python prefix-replacement path that mamba/micromamba do
# not share. If no solver is present, micromamba is bootstrapped into ~/.local/bin — a single
# static binary needing no admin rights.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; ROOT="$(cd "$HERE/.." && pwd)"
MODE="${1:-conda}"; FRESH=""
for a in "$@"; do [ "$a" = "--fresh" ] && FRESH=1; done
ENVNAME="pursue-bench"

bootstrap_micromamba() {
  echo ">> no conda/mamba/micromamba found — installing micromamba to ~/.local/bin (no admin needed)"
  mkdir -p "$HOME/.local/bin"
  curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xj -C "$HOME/.local" bin/micromamba
  export PATH="$HOME/.local/bin:$PATH"
  echo ">> add this to ~/.bashrc:  export PATH=\"\$HOME/.local/bin:\$PATH\""
}

if [ "$MODE" = "conda" ]; then
  if   command -v micromamba >/dev/null 2>&1; then SOLVER=micromamba
  elif command -v mamba      >/dev/null 2>&1; then SOLVER=mamba
  elif command -v conda      >/dev/null 2>&1; then SOLVER=conda
  else bootstrap_micromamba; SOLVER=micromamba; fi
  echo ">> solver: $SOLVER"

  if [ -n "$FRESH" ]; then
    echo ">> removing any existing $ENVNAME"
    $SOLVER env remove -y -n "$ENVNAME" 2>/dev/null || $SOLVER remove -y -n "$ENVNAME" --all 2>/dev/null || true
  fi

  # "prefix already exists" is not an error: update the existing env instead of failing.
  if [ "$SOLVER" = "micromamba" ]; then
    micromamba create -y -n "$ENVNAME" -f "$HERE/environment.yml" \
      || micromamba install -y -n "$ENVNAME" -f "$HERE/environment.yml"
    eval "$(micromamba shell hook -s bash)"; micromamba activate "$ENVNAME"
  else
    $SOLVER env create -n "$ENVNAME" -f "$HERE/environment.yml" \
      || $SOLVER env update -n "$ENVNAME" -f "$HERE/environment.yml"
    eval "$($SOLVER shell.bash hook)"; conda activate "$ENVNAME"
  fi
  echo ">> activate later with:  $SOLVER activate $ENVNAME"
else
  MOD="${2:-R}"; module load "$MOD" || true
  # Keep the user's existing library on the path: R_LIBS_USER is a colon-separated list, and
  # replacing it outright hides packages already installed on this cluster.
  BENCH_LIB="$HOME/R/pursue-bench-lib"; mkdir -p "$BENCH_LIB"
  BASE_LIB="$(R --no-echo -e 'cat(Sys.getenv("R_LIBS_USER"))' 2>/dev/null | tr ':' '\n' | grep -v "pursue-bench-lib" | paste -sd: -)"
  export R_LIBS_USER="$BENCH_LIB${BASE_LIB:+:$BASE_LIB}"
  echo ">> R library path: $R_LIBS_USER (export R_LIBS_USER before running jobs)"
fi

RV="$(Rscript -e 'cat(R.version.string)')"
echo ">> $RV"
case "$RV" in
  *"4.0"*|*"4.1"*|*"4.2"*|*"4.3"*|*"4.4"*)
    echo ">> WARNING: R < 4.5 — LOCOM2 (needs Deriv, R >= 4.5), MaAsLin 3 and ADAPT (Bioc >= 3.20)"
    echo ">>          will not install. Those are three of the comparators the benchmark needs." ;;
esac
Rscript "$HERE/install_packages.R"
echo ">> now: bash hpc/download_templates.sh && Rscript benchmarks/expected/make_expected.R && Rscript hpc/smoke_test.R"
