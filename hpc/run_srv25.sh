#!/usr/bin/env bash
# run_srv25.sh -- it54: plain-correlation permutation statistic for continuous exposure (no covariates).
#   srv25:  eperm_c5q (reference) and eperm_c10 on the three continuous settings -- house R22, msq R22, implant B_cont --
#           10 reps (nothing else differs between the two candidates);
#   srv25n: their nulls (house R22, msq R22, B_cont; da_frac 0) at 40 reps.
# bank3 locally: twinsuk B_cont 21.2 -> 32.8 TP per cell (ZicoSeq 37.8, LDM 33.8), hmp_tongue 3.6 -> 4.0 (5.8 / 5.4).
# hpc/guard.sh: shared memory cap for all our runs (MEMTOTAL, default 100G; MEMMAX per run), OOM-killer first choice,
# single-threaded BLAS. Cores from CORES. Resumable; pushes.
#   git pull && mkdir -p logs && (MEMTOTAL=130G CORES=16 nohup nice -n 10 bash hpc/run_srv25.sh > logs/run_srv25.log 2>&1 &)
source "$(dirname "${BASH_SOURCE[0]}")/guard.sh" "$@"
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
CORES="${CORES:-16}"
C=eperm_c5q,eperm_c10
L=srv25; R=benchmarks/dev/results/$L
echo ">> dev suite $L ($(date '+%F %T'))"
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq --settings house:R22 --extra msq:R22,implant:B_cont \
  --reps 10 --cores "$CORES" --label "$L" --candidates "$C" || exit 1
Rscript benchmarks/dev/compare.R --labels "$L" > "$R/compare.txt" 2>&1
LN=srv25n; RN=benchmarks/dev/results/$LN
echo ">> null run $LN ($(date '+%F %T'))"
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq --settings none --extra house:R22.null,msq:R22.null,implant:B_cont.null \
  --reps 40 --cores "$CORES" --label "$LN" --candidates "$C" \
  && Rscript benchmarks/dev/tools/null_summary.R "$RN/cells.csv" > "$RN/null_summary.txt" 2>&1
cp logs/run_srv25.log "$R/run.log" 2>/dev/null || true
git add "$R" "$RN" && git commit -m "dev suite $L / $LN: eperm_c10 (correlation statistic for continuous exposure)" && git pull --rebase --autostash && git push
echo ">> done $(date '+%F %T')"
