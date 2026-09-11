# Sourced by every sbatch script. Adapt to the cluster: either activate the conda env or
# load the R module and point R_LIBS_USER at the library created by setup_env.sh.
if command -v conda >/dev/null 2>&1 && conda env list 2>/dev/null | grep -q pursue-bench; then
  eval "$(conda shell.bash hook)"; conda activate pursue-bench
else
  module load R 2>/dev/null || true
  export R_LIBS_USER="${R_LIBS_USER:-$HOME/R/pursue-bench-lib}"
fi
export PURSUE_BENCH_ROOT="$PWD/benchmarks"
export PURSUE_DATA_ROOT="$PWD/benchmarks/data"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p logs results cache
