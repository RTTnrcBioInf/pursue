#!/usr/bin/env bash
# run_p15.sh -- full benchmark (seed 1, evaluation pool) for eperm_c10 where it can differ from eperm_c5q (p14).
# c10 = c5q with a plain-correlation permutation statistic for continuous exposure without covariates; every other
# design runs identical code (checked: identical p and q on binary cells). So: axis A R22 (continuous exposure) on all
# five simulators, axis B B_cont, plus R00 / B_ref replicate 1 as identity checks against p14's c5q. The c10 scoreboard is
# then p14's c5q cells with these replaced. Axis C (binary null permutations) cannot differ and is skipped.
# hpc/guard.sh: shared memory cap (MEMTOTAL), OOM-killer first choice, single-threaded BLAS; the scheduler stays within
# 90% of MEMTOTAL and counts cores without the OMP variables. Resumable (finished cells are skipped); pushes results/summary.
#   git pull && mkdir -p logs && (MEMTOTAL=130G nohup nice -n 10 bash hpc/run_p15.sh > logs/run_p15.log 2>&1 &)
source "$(dirname "${BASH_SOURCE[0]}")/guard.sh" "$@"
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs results
export PURSUE_BENCH_ROOT="$ROOT/benchmarks" PURSUE_DATA_ROOT="$ROOT/benchmarks/data" PURSUE_FOCUS="pursue03_epermc5q"
M=pursue03_epermc10; TAG=p15; T="hpc/tasks_$TAG"
echo ">> smoke: one house R22 cell ($(date '+%F %T'))"
SM=$(mktemp -d)
Rscript benchmarks/R/engine/run_cell.R --axis A --simulator house --template hmp_stool --regime R22 --replicate 1 --methods "$M" --out "$SM" --tag smoke 2>&1 | tee "$SM/smoke.log"
if ! grep -E "^  $M " "$SM/smoke.log" | grep -vqE 'error|not_installed'; then echo ">> smoke FAILED -- not launching"; exit 1; fi
rm -rf "$SM"
Rscript hpc/make_tasklist.R --pool evaluation --simulators house,msq,mid,sd2,sps --methods "$M" --regimes R22,R00 --out "$T"
grep -E -- "--regime R22 |--regime R00 --replicate 1 " "$T/axisA.txt" > "$T/axisA_run.txt"
grep -E -- "--regime B_cont |--regime B_ref --replicate 1 " "$T/axisB.txt" > "$T/axisB_run.txt"
wc -l "$T/axisA_run.txt" "$T/axisB_run.txt"
for ax in A B; do echo ">> axis $ax ($(date '+%F %T'))"; TAG="$TAG" bash hpc/run_adaptive.sh "$T/axis${ax}_run.txt"; done
echo ">> aggregate ($(date '+%F %T'))"
Rscript benchmarks/R/analysis/aggregate.R results
cp "logs/run_$TAG.log" "results/summary/run_$TAG.log" 2>/dev/null || true
git add results/summary && git commit -m "full benchmark: eperm_c10 on continuous exposure (axis A R22, axis B B_cont; tag $TAG)" \
  && git pull --rebase --autostash && git push
echo ">> done $(date '+%F %T')"
