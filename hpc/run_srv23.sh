#!/usr/bin/env bash
# run_srv23.sh -- it52: the unadjusted frequency kernel (lin0) and a sqrt / log omnibus as selection options (msq gap).
#   srv23:  eperm_c5q (reference, same cells as srv21 -- also a reproducibility check), eperm_c8 (+ lin0), eperm_c8s
#           (+ lin0 and max(|t_sqrt|, |t_log|)) on the whole dev suite + sd2 / sps + B_conf04 / B_conf07 / B_cont + msq R22;
#   srv23n: nulls at 40 reps (srv20n's six + house R22, B_cont: lin0 also enters the continuous path) and bloom x4 at 20.
# Single-threaded BLAS is exported (srv22: without it msq / mid draws change and runs stop being paired); the dev suite
# takes its core count from CORES, so nothing here counts cores after the export. Resumable; pushes its results.
#   git pull && mkdir -p logs && (nohup bash hpc/run_srv23.sh > logs/run_srv23.log 2>&1 &)
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
CORES="${CORES:-56}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
C=eperm_c5q,eperm_c8,eperm_c8s
L=srv23; R=benchmarks/dev/results/$L
echo ">> dev suite $L ($(date '+%F %T'))"
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq,mid,sd2,sps --extra implant:B_conf04,implant:B_conf07,implant:B_cont,msq:R22 \
  --reps 10 --cores "$CORES" --label "$L" --candidates "$C" || exit 1
Rscript benchmarks/dev/compare.R --labels "srv21,$L" > "$R/compare.txt" 2>&1
cp logs/run_srv23.log "$R/run.log" 2>/dev/null || true
git add "$R" && git commit -m "dev suite $L: eperm_c8 / c8s (unadjusted frequency kernel, sqrt-log omnibus)" && git pull --rebase --autostash && git push
LN=srv23n; RN=benchmarks/dev/results/$LN
echo ">> null run $LN ($(date '+%F %T'))"
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq,mid,sd2,sps --settings sd2:R06,sps:R06,house:R06,implant:B_null,msq:R06,mid:R06 \
  --extra house:R22.null,implant:B_cont.null --reps 40 --cores "$CORES" --label "$LN" --candidates "$C" \
  && Rscript benchmarks/dev/tools/null_summary.R "$RN/cells.csv" > "$RN/null_summary.txt" 2>&1
Rscript benchmarks/dev/devsuite.R --suite full --sims house --settings bloom:x4 --reps 20 --cores "$CORES" --label "${LN}b" --candidates "$C"
cp logs/run_srv23.log "$RN/run.log" 2>/dev/null || true
git add "$RN" "benchmarks/dev/results/${LN}b" && git commit -m "dev suite $LN: nulls and bloom for eperm_c8 / c8s" && git pull --rebase --autostash && git push
echo ">> done $(date '+%F %T')"
