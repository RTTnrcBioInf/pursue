#!/usr/bin/env bash
# run_p03.sh -- the outer loop for the PURSUE 0.3 candidates: full benchmark, axes A-C, evaluation
# pool, all five simulators, candidates ONLY, as a tagged solo re-run beside the stored comparator
# results (hpc/README.md section 7b). Then aggregates and pushes results/summary/.
#
#   nohup bash hpc/run_p03.sh > logs/run_p03.log 2>&1 &
#
# Resumable: re-running skips every cell whose tagged manifest exists. Memory-aware via
# run_adaptive.sh (MEM_PCT / MAXJOBS env as there). Candidates: benchmarks/R/methods/pursue03.R.
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs results
M="pursue03_erdl,pursue03_erdc,pursue03_erlc"; TAG="${TAG:-p03}"
export PURSUE_BENCH_ROOT="$ROOT/benchmarks" PURSUE_DATA_ROOT="$ROOT/benchmarks/data"

echo ">> smoke: one house cell, the three candidates"
SM=$(mktemp -d)
Rscript benchmarks/R/engine/run_cell.R --axis A --simulator house --template hmp_stool --regime R00 --replicate 1 \
  --methods "$M" --out "$SM" --tag smoke 2>&1 | tee "$SM/smoke.log"
if [ "$(grep -E '^  pursue03_[a-z]+ ' "$SM/smoke.log" | grep -vcE 'error|not_installed')" -ne 3 ]; then
  echo ">> smoke FAILED -- not launching the grid"; exit 1; fi
rm -rf "$SM"

Rscript hpc/make_tasklist.R --pool evaluation --simulators house,msq,mid,sd2,sps --methods "$M" --out hpc/tasks_p03
for ax in A B C; do echo ">> axis $ax"; TAG="$TAG" bash hpc/run_adaptive.sh "hpc/tasks_p03/axis$ax.txt"; done

echo ">> aggregate"
Rscript benchmarks/R/analysis/aggregate.R results
cp logs/run_p03.log results/summary/run_p03.log 2>/dev/null || true
git add results/summary && git commit -m "full benchmark, axes A-C: PURSUE 0.3 candidates (tag $TAG)" \
  && git pull --rebase --autostash && git push
echo ">> done"
