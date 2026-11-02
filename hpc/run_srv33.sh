#!/usr/bin/env bash
# run_srv33.sh -- c14 (c12 + the corrected within-subject path, benchmarks/dev/candidates/84_c14.R) checked where srv32 left
# questions, on development data only (rd-charter 2026-10-30). Two suites, each pushed when done:
#   srv33a -- two samples per subject, more reps where c13 was worst: within-subject nulls (120 reps) and signal (20) on
#             thin-hmp_tongue (one member thinned to 10% depth), rm_throat, rm_baboon and a scrambled-site stress
#             (mix-hmp_tongue-scr-rm_buccal: tongue paired with buccal whose taxon labels are permuted -- the server version
#             of the local check that reproduced the ear's size of failure). Reps 1-40 repeat srv32's cells exactly.
#   srv33b -- more than two samples per subject, exactly half exposed (so the exact within-subject path is taken):
#             baboons with 4 samples per animal (P_pair4: 80 null, 20 signal reps) and house's balanced crossover
#             (HW4h: 4 visits, 2 exposed, subject intercept ICC 0.3; tongue and twinsuk templates).
# Candidates: c12, c13, c14, c13g0 (no centre: exact reference), LinDA (random intercept), LDM (cluster). Dev seed 9001.
# hpc/guard.sh: shared memory cap (MEMTOTAL), OOM-killer first choice, single-threaded BLAS. Resumable; pushes.
#   git pull && mkdir -p logs && (MEMTOTAL=90G CORES=12 nohup nice -n 10 bash hpc/run_srv33.sh > logs/run_srv33.log 2>&1 &)
source "$(dirname "${BASH_SOURCE[0]}")/guard.sh" "$@"
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
CORES="${CORES:-12}"
C=eperm_c12,eperm_c13,eperm_c14,eperm_c13g0,ref_linda_re,ref_ldm_cl
suite() {  # suite <label> <message> <rm_suite arguments...>
  local L="$1" msg="$2"; shift 2; local R=benchmarks/dev/results/$L
  echo ">> rm suite $L ($(date '+%F %T'))"
  Rscript benchmarks/dev/rm_suite.R --candidates "$C" --cores "$CORES" --label "$L" "$@" || { echo ">> $L failed"; return 1; }
  cp logs/run_srv33.log "$R/run.log" 2>/dev/null || true
  git add "$R" && git commit -q -m "$msg" && git pull -q --rebase --autostash && git push -q
}
suite srv33a "rm suite srv33a: c14 on two-sample pairs, more reps where c13 was worst" \
  --settings P_pair.null,P_pair --template thin-hmp_tongue,rm_throat,rm_baboon,mix-hmp_tongue-scr-rm_buccal \
  --pick-random rm_baboon --reps 20 --null-reps 120
suite srv33b "rm suite srv33b: c14 with four samples per subject (baboons, house crossover)" \
  --settings P_pair4.null,P_pair4,HW4h.null,HW4h --template rm_baboon --pick-random rm_baboon --reps 20 --null-reps 80
echo ">> done $(date '+%F %T')"
