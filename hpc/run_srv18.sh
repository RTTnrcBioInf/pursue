#!/usr/bin/env bash
# run_srv18.sh -- Phase 2, it42 + it44 (benchmarks/dev/candidates/70_epg.R, 71_eps.R): the exact permutation path for
# every design (eperm_g) and detection alone as a selection option (eperm_gd), against eperm_cs1.
#   srv18:  the whole dev suite (house, implant, msq, mid; core + ext) plus the implant specs the suite lacked
#           (B_conf04, B_conf07, B_cont), 10 reps -- power and the usual gate;
#   srv18n: nulls of the designs that newly take the permutation path (depth-confounded R17/R19, confounder R21,
#           continuous R22, repeated R23, on house / msq / mid, and implant B_conf07 / B_cont) next to the
#           reference nulls (house R06, implant B_null, msq R06, mid R06), 20 reps -- P(any) and null FPR.
# Development only: the benchmark does not source these candidate files, so the conf1 confirmation is untouched.
# 16 cores by default (conf1 A-C uses 40, D/E 8). Resumable per cell; pushes its results.
#   git pull && mkdir -p logs && (nohup bash hpc/run_srv18.sh > logs/run_srv18.log 2>&1 &)
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
CORES="${CORES:-16}"; C=eperm_cs1,eperm_g,eperm_gd
L=srv18; R=benchmarks/dev/results/$L
echo ">> dev suite $L ($(date '+%F %T'))"
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq,mid --extra implant:B_conf04,implant:B_conf07,implant:B_cont \
  --reps 10 --cores "$CORES" --label "$L" --candidates "$C" || exit 1
Rscript benchmarks/dev/compare.R --labels "$L" > "$R/compare.txt" 2>&1
cp logs/run_srv18.log "$R/run.log" 2>/dev/null || true
git add "$R" && git commit -m "dev suite $L: it42/it44 -- permutation path for every design, detection as an option" \
  && git pull --rebase --autostash && git push
LN=srv18n; RN=benchmarks/dev/results/$LN
echo ">> null run $LN ($(date '+%F %T'))"
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq,mid --settings house:R06,implant:B_null,msq:R06,mid:R06 \
  --extra house:R17.null,house:R19.null,house:R21.null,house:R22.null,house:R23.null,msq:R19.null,mid:R17.null,mid:R19.null,implant:B_conf07.null,implant:B_cont.null \
  --reps 20 --cores "$CORES" --label "$LN" --candidates "$C" \
  && Rscript benchmarks/dev/tools/null_summary.R "$RN/cells.csv" > "$RN/null_summary.txt" 2>&1
cp logs/run_srv18.log "$RN/run.log" 2>/dev/null || true
git add "$RN" && git commit -m "dev suite $LN: nulls for the it42 permutation paths" && git pull --rebase --autostash && git push
echo ">> done $(date '+%F %T')"
