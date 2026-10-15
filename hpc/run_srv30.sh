#!/usr/bin/env bash
# run_srv30.sh -- repeated measures on real data, first look (2026-10-15). hmp_tongue (tuning template) gets HMP subject
# and visit labels (dev/tools/hmp_subjects.R -> benchmarks/devdata/hmp16s_subjects.csv, from the HMP16SData sites the
# longitudinal fetch saved). Then dev/rm_suite.R: 40 subjects x 2 visits per cell, exposure fixed per subject (P_subj)
# or varying within subject (P_pair), with implanted signal (20 reps) and as exact nulls (40 reps). c12 beside each
# comparator twice: as the benchmark runs it (no subject structure) and with the subject handling its package supports
# (LinDA random intercept, LDM cluster.id, ZicoSeq strata, logistic presence with a subject-clustered sandwich).
# hpc/guard.sh: shared memory cap (MEMTOTAL), OOM-killer first choice, single-threaded BLAS. Resumable; pushes.
#   git pull && mkdir -p logs && (MEMTOTAL=130G CORES=8 nohup nice -n 10 bash hpc/run_srv30.sh > logs/run_srv30.log 2>&1 &)
source "$(dirname "${BASH_SOURCE[0]}")/guard.sh" "$@"
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
CORES="${CORES:-8}"
L=srv30; R=benchmarks/dev/results/$L; mkdir -p "$R"
echo ">> subject map ($(date '+%F %T'))"
Rscript benchmarks/dev/tools/hmp_subjects.R 2>&1 | tee "$R/subject_map.txt" || exit 1
C=eperm_c12,ref_linda,ref_linda_re,ref_ldm,ref_ldm_cl,ref_zicoseq,ref_zicoseq_st,ref_adapt,ref_logistic_presence,ref_logistic_cr
echo ">> rm suite $L ($(date '+%F %T'))"
Rscript benchmarks/dev/rm_suite.R --candidates "$C" --reps 20 --null-reps 40 --cores "$CORES" --label "$L" || exit 1
cp logs/run_srv30.log "$R/run.log" 2>/dev/null || true
git add benchmarks/devdata/hmp16s_subjects.csv "$R" && git commit -q -m "rm suite srv30: repeated measures on hmp_tongue with HMP subject labels" \
  && git pull -q --rebase --autostash && git push -q
echo ">> done $(date '+%F %T')"
