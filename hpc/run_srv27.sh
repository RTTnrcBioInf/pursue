#!/usr/bin/env bash
# run_srv27.sh -- spread-only implants on the tuning templates (benchmarks/dev/tools/spread.R, spread_suite.R).
# The real-data transfer check (2026-10-13) found one realistic stressor the scoreboard does not see: groups that differ
# in spread (R25, house only; c10 FDR 0.19, best of the field). Here real samples get a spread difference with no
# location change -- 30 / 70 and 50 / 50 groups, 20% or all taxa, presence kept or redrawn, and a continuous exposure
# whose spread trends with z -- so every call is false. c10 beside ZicoSeq, LDM, LinDA, ADAPT, logistic presence,
# Wilcoxon and a t-test on log proportions; then c5q (HC0 t) on the continuous setting, paired with c10 (plain
# correlation, it54), which was never checked against exposure-dependent spread. 20 reps per template.
# hpc/guard.sh: shared memory cap (MEMTOTAL), OOM-killer first choice, single-threaded BLAS. Resumable; pushes.
#   git pull && mkdir -p logs && (MEMTOTAL=130G CORES=40 nohup nice -n 10 bash hpc/run_srv27.sh > logs/run_srv27.log 2>&1 &)
source "$(dirname "${BASH_SOURCE[0]}")/guard.sh" "$@"
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
CORES="${CORES:-40}"
TR="$(printf 'Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>\n')"
C=eperm_c10,ref_zicoseq,ref_ldm,ref_linda,ref_adapt,ref_logistic_presence,ref_wilcoxon_tss,ref_lm_logtss
L=srv27; R=benchmarks/dev/results/$L
echo ">> spread suite $L ($(date '+%F %T'))"
Rscript benchmarks/dev/spread_suite.R --candidates "$C" --reps 20 --cores "$CORES" --label "$L" || exit 1
L2=srv27c; R2=benchmarks/dev/results/$L2
echo ">> spread suite $L2: c5q on H_cont ($(date '+%F %T'))"
Rscript benchmarks/dev/spread_suite.R --candidates eperm_c5q --settings H_cont --reps 20 --cores "$CORES" --label "$L2" || exit 1
cp logs/run_srv27.log "$R/run.log" 2>/dev/null || true
git add "$R" "$R2" && git commit -m "dev suite $L: spread-only implants (c10 beside seven references; c5q on H_cont)" -m "$TR" \
  && git pull --rebase --autostash && git push
echo ">> done $(date '+%F %T')"
