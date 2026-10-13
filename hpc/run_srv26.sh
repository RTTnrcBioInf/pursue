#!/usr/bin/env bash
# run_srv26.sh -- it54 null check at scale. srv25: c10 (correlation statistic, continuous exposure) gains B_cont 13.0 ->
# 19.6, msq R22 8.65 -> 10.0 TP per cell, house R22 level; nulls P(any) c5q 0.042 vs c10 0.058 over 240 cells, house
# R22.null 5 -> 8 of 80 cells (bulk FPR equal, 0.052 / 0.051). The three continuous nulls at 160 reps per template; the
# first 40 are srv25n's cells (same seeds; their checkpoints are reused), so 120 are new.
# hpc/guard.sh: shared memory cap (MEMTOTAL), OOM-killer first choice, single-threaded BLAS. Resumable; pushes.
#   git pull && mkdir -p logs && (MEMTOTAL=130G CORES=24 nohup nice -n 10 bash hpc/run_srv26.sh > logs/run_srv26.log 2>&1 &)
source "$(dirname "${BASH_SOURCE[0]}")/guard.sh" "$@"
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
CORES="${CORES:-24}"
LN=srv26n; RN=benchmarks/dev/results/$LN
echo ">> null run $LN ($(date '+%F %T'))"
mkdir -p "$RN/cells" && cp -n benchmarks/dev/results/srv25n/cells/*.csv "$RN/cells/" 2>/dev/null || true
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq --settings none --extra house:R22.null,msq:R22.null,implant:B_cont.null \
  --reps 160 --cores "$CORES" --label "$LN" --candidates eperm_c5q,eperm_c10 \
  && Rscript benchmarks/dev/tools/null_summary.R "$RN/cells.csv" > "$RN/null_summary.txt" 2>&1
cp logs/run_srv26.log "$RN/run.log" 2>/dev/null || true
git add "$RN" && git commit -m "dev suite $LN: continuous-exposure nulls at 160 reps for eperm_c5q / c10" && git pull --rebase --autostash && git push
echo ">> done $(date '+%F %T')"
