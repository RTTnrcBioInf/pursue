# Sourced by every sbatch script. Adapt to the cluster: either activate the conda env or
# load the R module and point R_LIBS_USER at the library created by setup_env.sh.
if command -v conda >/dev/null 2>&1 && conda env list 2>/dev/null | grep -q pursue-bench; then
  eval "$(conda shell.bash hook)"; conda activate pursue-bench
else
  module load R 2>/dev/null || true
  # Keep the user's existing library on the path: R_LIBS_USER is a colon-separated list, and
  # replacing it outright hides packages already installed on this cluster (that is what hid
  # MicrobiomeStat / corncob / LDM / LOCOM in the 2026-09-11 smoke test). New installs go to
  # the first entry; everything already installed stays visible.
  BENCH_LIB="$HOME/R/pursue-bench-lib"; mkdir -p "$BENCH_LIB"
  BASE_LIB="$(R --no-echo -e 'cat(Sys.getenv("R_LIBS_USER"))' 2>/dev/null | tr ':' '\n' | grep -v "pursue-bench-lib" | paste -sd: -)"
  export R_LIBS_USER="$BENCH_LIB${BASE_LIB:+:$BASE_LIB}"
fi
export PURSUE_BENCH_ROOT="$PWD/benchmarks"
export PURSUE_DATA_ROOT="$PWD/benchmarks/data"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p logs results cache
