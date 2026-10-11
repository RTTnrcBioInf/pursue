#!/usr/bin/env bash
# run_srv21.sh -- it48 after srv20 (eperm_c4: sd2 R08 1.15 -> 4.20, R13 3.40 -> 7.65 TP per cell; null P(any) c3 0.052 ->
# c4 0.061, every extra null cell exactly two false calls = a spurious pair selecting its own option).
#   srv21:  eperm_c5   = c4 + an option needs >= 2 leave-one-out discoveries before it is used (no self-supporting pairs);
#           eperm_c5q  = c5 + Storey's adaptive BH, applied only when plain BH rejects (global-null P(any) = BH's);
#           eperm_c5w  = c5q without winsorising (banked msq R00 / R11 tongue +0.6 / +0.4 TP per cell);
#           eperm_c5ws = c5w without size factors (banked msq R00 tongue p_full 9.4 -> 10.2, linear kernel 5.8 -> 8.6);
#           whole dev suite + sd2 / sps + implant B_conf04 / B_conf07 / B_cont + msq R22, 10 reps (same cells as srv20);
#   srv21n: nulls at 40 reps (srv20n's six + mid R17 / R19, implant B_conf07 / B_cont) and the x4 bloom at 20.
# Resumable; pushes its results.
#   git pull && mkdir -p logs && (nohup bash hpc/run_srv21.sh > logs/run_srv21.log 2>&1 &)
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
CORES="${CORES:-48}"; C=eperm_c5,eperm_c5q,eperm_c5w,eperm_c5ws
L=srv21; R=benchmarks/dev/results/$L
echo ">> dev suite $L ($(date '+%F %T'))"
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq,mid,sd2,sps --extra implant:B_conf04,implant:B_conf07,implant:B_cont,msq:R22 \
  --reps 10 --cores "$CORES" --label "$L" --candidates "$C" || exit 1
Rscript benchmarks/dev/compare.R --labels "srv20,$L" > "$R/compare.txt" 2>&1
cp logs/run_srv21.log "$R/run.log" 2>/dev/null || true
git add "$R" && git commit -m "dev suite $L: eperm_c5 / c5q / c5w / c5ws (selection floor, gated adaptive BH, no winsorising, no size factors)" && git pull --rebase --autostash && git push
LN=srv21n; RN=benchmarks/dev/results/$LN
echo ">> null run $LN ($(date '+%F %T'))"
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq,mid,sd2,sps --settings sd2:R06,sps:R06,house:R06,implant:B_null,msq:R06,mid:R06 \
  --extra mid:R17.null,mid:R19.null,implant:B_conf07.null,implant:B_cont.null \
  --reps 40 --cores "$CORES" --label "$LN" --candidates "$C" \
  && Rscript benchmarks/dev/tools/null_summary.R "$RN/cells.csv" > "$RN/null_summary.txt" 2>&1
Rscript benchmarks/dev/devsuite.R --suite full --sims house --settings bloom:x4 --reps 20 --cores "$CORES" --label "${LN}b" --candidates "$C"
cp logs/run_srv21.log "$RN/run.log" 2>/dev/null || true
git add "$RN" "benchmarks/dev/results/${LN}b" && git commit -m "dev suite $LN: nulls and bloom for eperm_c5 family" && git pull --rebase --autostash && git push
echo ">> done $(date '+%F %T')"
