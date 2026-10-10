#!/usr/bin/env bash
# run_srv20.sh -- Phase 2 after the conf1 confirmation (eperm_cs1 held up: 0 failing cells at a fresh seed).
#   probe1: what makes LDM and ZicoSeq stronger on msq / sd2 -- their installed sources, and LDM (each scale's p / F / q
#           kept apart) and ZicoSeq re-run on the banked sd2 / msq cells (dev/tools/probe_comparators.R);
#   srv20:  eperm_c4 (eperm_c3 + an abundance-given-presence kernel as a selection option -- PURSUE 0.2's abundance arm
#           has 3-5x our sd2 power at both seeds; locally on banked sd2 cells c3 -> c4: R08 2.0 -> 10.0, R13 6.0 -> 17.0
#           TP per cell, msq unchanged) and eperm_c4q (+ adaptive BH), with eperm_c3 as the reference, on the whole dev
#           suite + sd2 / sps + implant B_conf04 / B_conf07 / B_cont + msq R22, 10 reps;
#   srv20n: nulls at 40 reps (sd2 R06, sps R06, house R06, implant B_null, msq R06, mid R06) and the x4 bloom at 20.
# The confirmation is finished, so the whole machine is free: CORES defaults to 48. Resumable; pushes its results.
#   git pull && mkdir -p logs && (nohup bash hpc/run_srv20.sh > logs/run_srv20.log 2>&1 &)
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
CORES="${CORES:-48}"; C=eperm_c3,eperm_c4,eperm_c4q
O=benchmarks/dev/results/probe1
echo ">> probe ($(date '+%F %T'))"
CORES=16 Rscript benchmarks/dev/tools/probe_comparators.R benchmarks/dev/bank/bank1 "$O" && ls -la "$O"
git add "$O" && git commit -m "probe1: LDM / ZicoSeq sources and per-scale outputs on bank1" && git pull --rebase --autostash && git push
L=srv20; R=benchmarks/dev/results/$L
echo ">> dev suite $L ($(date '+%F %T'))"
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq,mid,sd2,sps --extra implant:B_conf04,implant:B_conf07,implant:B_cont,msq:R22 \
  --reps 10 --cores "$CORES" --label "$L" --candidates "$C" || exit 1
Rscript benchmarks/dev/compare.R --labels "srv19,$L" > "$R/compare.txt" 2>&1
cp logs/run_srv20.log "$R/run.log" 2>/dev/null || true
git add "$R" && git commit -m "dev suite $L: eperm_c4 / c4q (abundance-given-presence kernel)" && git pull --rebase --autostash && git push
LN=srv20n; RN=benchmarks/dev/results/$LN
echo ">> null run $LN ($(date '+%F %T'))"
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq,mid,sd2,sps --settings sd2:R06,sps:R06,house:R06,implant:B_null,msq:R06,mid:R06 \
  --reps 40 --cores "$CORES" --label "$LN" --candidates "$C" \
  && Rscript benchmarks/dev/tools/null_summary.R "$RN/cells.csv" > "$RN/null_summary.txt" 2>&1
Rscript benchmarks/dev/devsuite.R --suite full --sims house --settings bloom:x4 --reps 20 --cores "$CORES" --label "${LN}b" --candidates "$C"
cp logs/run_srv20.log "$RN/run.log" 2>/dev/null || true
git add "$RN" "benchmarks/dev/results/${LN}b" && git commit -m "dev suite $LN: nulls and bloom for eperm_c4 / c4q" && git pull --rebase --autostash && git push
echo ">> done $(date '+%F %T')"
