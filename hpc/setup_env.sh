#!/usr/bin/env bash
# Create the benchmark environment on the HPC, then install every R package.
#
#   bash hpc/setup_env.sh                  # conda-family env named pursue-bench (recommended)
#   --keep   reuse an existing pursue-bench instead of rebuilding it
#   (default: the env is always rebuilt from scratch)
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
KEEP=""
for a in "$@"; do [ "$a" = "--fresh" ] && FRESH=1; [ "$a" = "--keep" ] && KEEP=1; done
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

  # The env is ALWAYS built from scratch, never updated in place.
  # Updating an existing pursue-bench is what broke on 2026-09-11: the old env still held
  # bioconda bioconductor-* builds from the R 4.3 attempt, bumping r-base to 4.5 forced their
  # r45 rebuilds, and their post-link scripts ran `R CMD INSTALL` against a half-updated R
  # ("ERROR: loading failed for 'R', 'R.c~'"). conda's own logger then crashed formatting the
  # error ("unsupported format character 'T'"), hiding the real cause. A clean env built from
  # this conda-forge-only file never runs a bioconda post-link script at all.
  #
  # `$SOLVER env remove` is not trusted here: it reported success while leaving the prefix in
  # place, after which `env create` failed with "prefix already exists". Delete the directory.
  PREFIX="$($SOLVER env list 2>/dev/null | awk -v n="$ENVNAME" '$1 == n {print $NF}' | head -1)"
  [ -z "${PREFIX:-}" ] && for base in "$HOME/.conda/envs" "$HOME/micromamba/envs" "${MAMBA_ROOT_PREFIX:-}/envs" "${CONDA_PREFIX:-}/envs"; do
    [ -n "$base" ] && [ -d "$base/$ENVNAME" ] && PREFIX="$base/$ENVNAME" && break
  done
  if [ -n "${PREFIX:-}" ] && [ -d "$PREFIX" ] && [ -n "$KEEP" ]; then
    echo ">> --keep: reusing the existing environment at $PREFIX (skipping rebuild)"
    PREFIX=""
  fi
  if [ -n "${PREFIX:-}" ] && [ -d "$PREFIX" ]; then
    echo ">> removing existing environment at $PREFIX"
    $SOLVER env remove -y -n "$ENVNAME" >/dev/null 2>&1 || true
    rm -rf "$PREFIX"
    [ -d "$PREFIX" ] && { echo ">> ERROR: could not delete $PREFIX -- remove it by hand and re-run"; exit 1; }
  fi

  if [ "$SOLVER" = "micromamba" ]; then
    micromamba create -y -n "$ENVNAME" -f "$HERE/environment.yml"
  else
    $SOLVER env create -n "$ENVNAME" -f "$HERE/environment.yml"
  fi

  # Activation is NOT "$SOLVER shell.bash hook": that spelling is conda's alone. mamba rejects
  # it ("invalid choice: 'shell.bash'"), the eval then yields nothing, and the bare activate
  # fails with "Run 'conda init' before 'conda activate'". mamba v1/v2 create ordinary conda
  # envs, so activate them through conda's hook; micromamba and standalone mamba v2 use
  # `shell hook -s bash` instead.
  if [ "$SOLVER" = "micromamba" ]; then
    eval "$(micromamba shell hook -s bash)"; micromamba activate "$ENVNAME"
  elif command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"; conda activate "$ENVNAME"
  else
    eval "$(mamba shell hook -s bash)"; mamba activate "$ENVNAME"
  fi
  ACT=$([ "$SOLVER" = "micromamba" ] && echo micromamba || echo conda)
  echo ">> activate later with:  $ACT activate $ENVNAME"
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
