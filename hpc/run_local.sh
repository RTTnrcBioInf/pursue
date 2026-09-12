#!/usr/bin/env bash
# Run a benchmark task list on ONE machine, N cells at a time. Replaces the SLURM array jobs
# where no scheduler exists; the task lists and the engine are identical either way.
#
#   bash hpc/run_local.sh hpc/tasks_pilot/axisA.txt            # N = all cores
#   bash hpc/run_local.sh hpc/tasks/axisA.txt 32               # N = 32 concurrent cells
#   DRY=1 bash hpc/run_local.sh hpc/tasks/axisA.txt 32         # show what would run
#
# Each cell is an independent Rscript, so concurrency is just how many run at once. Resumable:
# a cell whose manifest already exists is skipped, so re-running after an interruption (or
# after adding replicates) only does the missing work.
set -euo pipefail
TASKS="${1:?usage: bash hpc/run_local.sh <tasklist.txt> [n_concurrent]}"
NJOBS="${2:-$(nproc)}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
[ -f "$TASKS" ] || { echo "no such task list: $TASKS"; exit 1; }

# ONE thread per cell. Without this each of N concurrent R processes spawns its own BLAS
# threads and the machine thrashes -- the classic way to make 32 jobs slower than 8.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PURSUE_BENCH_ROOT="$ROOT/benchmarks" PURSUE_DATA_ROOT="$ROOT/benchmarks/data"
mkdir -p logs results cache

RV="$(Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || echo 0)"
case "$RV" in 4.5*|4.6*|5.*) ;; *) echo ">> WARNING: R $RV -- did you forget 'conda activate pursue-bench'?";; esac
Rscript -e 'if (!requireNamespace("PURSUE", quietly=TRUE)) { cat(">> ERROR: PURSUE is not installed in this library\n"); quit(status=1) }'

# Resume: derive the cell id the engine will write and skip cells already finished.
runner() {
  line="$1"
  ax=$(sed -n 's/.*--axis \([^ ]*\).*/\1/p' <<<"$line")
  sim=$(sed -n 's/.*--simulator \([^ ]*\).*/\1/p' <<<"$line")
  tpl=$(sed -n 's/.*--template \([^ ]*\).*/\1/p' <<<"$line")
  reg=$(sed -n 's/.*--regime \([^ ]*\).*/\1/p' <<<"$line")
  rep=$(sed -n 's/.*--replicate \([0-9]*\).*/\1/p' <<<"$line")
  tag="${TAG:-}"
  id=$(printf "%s__%s__%s__%s__r%03d%s" "$ax" "$sim" "$tpl" "$reg" "$rep" "${tag:+__$tag}")
  out="results/axis$ax"
  if [ -f "$out/$id.manifest.json" ]; then echo "skip  $id"; return 0; fi
  if [ -n "${DRY:-}" ]; then echo "would run  $id"; return 0; fi
  start=$SECONDS
  if Rscript benchmarks/R/engine/run_cell.R $line --out "$out" --cache cache \
       --master-seed "${PURSUE_MASTER_SEED:-1}" --timeout "${CELL_TIMEOUT:-3600}" \
       ${tag:+--tag "$tag"} > "logs/$id.log" 2>&1; then
    echo "ok    $id  $((SECONDS-start))s"
  else
    echo "FAIL  $id  $((SECONDS-start))s  (see logs/$id.log)"
  fi
}
export -f runner

total=$(grep -cve '^\s*$' "$TASKS")
echo ">> $TASKS: $total cells, $NJOBS at a time, R $RV"
echo ">> started $(date '+%F %T')"
grep -ve '^\s*$' "$TASKS" | xargs -d '\n' -P "$NJOBS" -I{} bash -c 'runner "$@"' _ {}
echo ">> finished $(date '+%F %T')"
done_n=$(ls results/axis*/ 2>/dev/null | grep -c 'manifest.json' || true)
echo ">> $done_n cells have results. Re-run this command to retry anything that failed."
