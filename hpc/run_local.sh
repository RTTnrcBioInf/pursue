#!/usr/bin/env bash
# Run a benchmark task list on ONE machine, N cells at a time. Replaces the SLURM array jobs
# where no scheduler exists; the task lists and the engine are identical either way.
#
#   bash hpc/run_local.sh hpc/tasks_pilot/axisA.txt            # N = all cores
#
# Cells are NOT equal in memory: a 2044-feature template costs many times a 300-feature one, and
# ZicoSeq is the peak. Split a task list by template and give the heavy half a smaller N:
#   grep -E 'risk_stool|hmp_stool' hpc/tasks/axisA.txt > hpc/tasks/axisA_big.txt
#   grep -vE 'risk_stool|hmp_stool' hpc/tasks/axisA.txt > hpc/tasks/axisA_small.txt
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
# GNU time gives peak RSS per cell. Without it the ok/FAIL lines carry no memory information,
# which is what left the 2026-09-12 silent kills undiagnosable.
TIME_BIN=""; [ -x /usr/bin/time ] && /usr/bin/time -f '%M' true 2>/dev/null && TIME_BIN=/usr/bin/time
[ -n "$TIME_BIN" ] || echo ">> note: /usr/bin/time not usable; per-cell memory will not be reported"
export TIME_BIN

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
  # A cell whose regime the simulator cannot produce writes <id>.skipped.json and no manifest.
  # Checking only for the manifest meant every resume re-ran all ~3125 of them at ~80s of R
  # startup each (~70 CPU-h per resume) and they could never become "done".
  if [ -f "$out/$id.manifest.json" ] || [ -f "$out/$id.skipped.json" ]; then echo "skip  $id"; return 0; fi
  if [ -n "${DRY:-}" ]; then echo "would run  $id"; return 0; fi
  start=$SECONDS
  mf="logs/$id.mem"
  # `rc=$?` on its own line would abort under `set -e`; capture it through the if instead.
  if [ -n "$TIME_BIN" ]; then
    if "$TIME_BIN" -f '%M' -o "$mf" Rscript benchmarks/R/engine/run_cell.R $line --out "$out" --cache cache \
         --master-seed "${PURSUE_MASTER_SEED:-1}" --timeout "${CELL_TIMEOUT:-3600}" \
         ${tag:+--tag "$tag"} > "logs/$id.log" 2>&1; then rc=0; else rc=$?; fi
  else
    if Rscript benchmarks/R/engine/run_cell.R $line --out "$out" --cache cache \
         --master-seed "${PURSUE_MASTER_SEED:-1}" --timeout "${CELL_TIMEOUT:-3600}" \
         ${tag:+--tag "$tag"} > "logs/$id.log" 2>&1; then rc=0; else rc=$?; fi
  fi
  # peak RSS in MB; a killed process still leaves the line GNU time already wrote
  mem=$( [ -f "$mf" ] && awk 'END{if($1+0>0) printf "%.0fMB", $1/1024; else print "?"}' "$mf" || echo "?" )
  if [ "$rc" -eq 0 ]; then
    echo "ok    $id  $((SECONDS-start))s  $mem"
  else
    # exit 137 = SIGKILL (the OOM killer); 139 = segfault. Both die with no R error in the log,
    # which is exactly what the 44 sd2/mid failures on risk_stool looked like.
    case "$rc" in 137) why=" KILLED(oom?)";; 139) why=" SEGFAULT";; *) why="";; esac
    echo "FAIL  $id  $((SECONDS-start))s  $mem  rc=$rc$why  (see logs/$id.log)"
  fi
}
export -f runner

total=$(grep -cve '^\s*$' "$TASKS")
echo ">> $TASKS: $total cells, $NJOBS at a time, R $RV"
echo ">> started $(date '+%F %T')"
grep -ve '^\s*$' "$TASKS" | xargs -d '\n' -P "$NJOBS" -I{} bash -c 'runner "$@"' _ {}
echo ">> finished $(date '+%F %T')"
done_n=$(ls results/axis*/ 2>/dev/null | grep -c 'manifest.json' || true)
skip_n=$(ls results/axis*/ 2>/dev/null | grep -c 'skipped.json' || true)
echo ">> $done_n cells have results; $skip_n regimes that simulator cannot produce."
echo ">> Re-run this command to retry anything that failed."
