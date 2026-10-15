#!/usr/bin/env bash
# run_srv31.sh -- c13's within-subject path with more than two visits (2026-10-15). house's generator with R23's
# subject structure (per-taxon subject intercept, ICC 0.3) but the exposure assigned within subject: HW5 (5 visits,
# 1-4 of them exposed per subject) and HWv (2-8 visits, 20% of subjects never exposed); signal (10 reps) and nulls
# (40 reps) on both tuning templates. A local smoke run found the first version of the path miscalibrated with > 2
# visits (its pair graph depended on the labels); fixed, two-visit results unchanged. c12, c13 and the comparators,
# with and without subject structure.
# hpc/guard.sh: shared memory cap (MEMTOTAL), OOM-killer first choice, single-threaded BLAS. Resumable; pushes.
#   git pull && mkdir -p logs && (MEMTOTAL=130G CORES=12 nohup nice -n 10 bash hpc/run_srv31.sh > logs/run_srv31.log 2>&1 &)
source "$(dirname "${BASH_SOURCE[0]}")/guard.sh" "$@"
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
CORES="${CORES:-12}"
L=srv31; R=benchmarks/dev/results/$L
C=eperm_c12,eperm_c13,ref_linda,ref_linda_re,ref_ldm,ref_ldm_cl,ref_zicoseq,ref_zicoseq_st,ref_adapt,ref_logistic_presence,ref_logistic_cr
echo ">> rm suite $L ($(date '+%F %T'))"
Rscript benchmarks/dev/rm_suite.R --candidates "$C" --settings HW5,HWv,HW5.null,HWv.null --reps 10 --null-reps 40 --cores "$CORES" --label "$L" || exit 1
cp logs/run_srv31.log "$R/run.log" 2>/dev/null || true
git add "$R" && git commit -q -m "rm suite srv31: within-subject exposure with more than two visits (house), c13 beside c12 and comparators" \
  && git pull -q --rebase --autostash && git push -q
echo ">> done $(date '+%F %T')"
