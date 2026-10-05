#!/usr/bin/env bash
# run_srv17.sh -- parallel development track: is p13's sd2 / sps gap visible on the TUNING templates?
# eperm_cs1 and the reference comparators (ZicoSeq, LDM, ADAPT, LinDA) on the new sd2 / sps dev settings
# (R00, R06, R08, R13 x hmp_tongue, twinsuk_stool x 5 reps), then msq R00/R11 for the same references at 10
# reps. 16 cores by default, so it can share the machine with hpc/run_confirm.sh (40). Pushes its results.
#   git pull && mkdir -p logs && (nohup bash hpc/run_srv17.sh > logs/run_srv17.log 2>&1 &)
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
L=srv17; R=benchmarks/dev/results/$L
C=eperm_cs1,ref_zicoseq,ref_ldm,ref_adapt,ref_linda
Rscript benchmarks/dev/devsuite.R --suite full --sims sd2,sps --reps 5 --cores "${CORES:-16}" --label "$L" --candidates "$C"
Rscript benchmarks/dev/devsuite.R --suite full --sims msq --settings msq:R00,msq:R11 --reps 10 --cores "${CORES:-16}" --label "${L}m" --candidates "$C"
{ Rscript benchmarks/dev/compare.R --labels "$L" --brief; Rscript benchmarks/dev/compare.R --labels "${L}m" --brief; } > "$R/compare.txt" 2>&1
cp logs/run_srv17.log "$R/run.log" 2>/dev/null || true
git add "$R" "benchmarks/dev/results/${L}m" && git commit -m "dev suite $L: sd2 / sps / msq gap on the tuning templates (eperm_cs1 vs references)" \
  && git pull --rebase --autostash && git push
