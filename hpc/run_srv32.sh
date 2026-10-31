#!/usr/bin/env bash
# run_srv32.sh -- why c13's within-subject path is not exact (conf2: hmp_skin_ear within-subject nulls, P(any) 0.45), on
# development data only (rd-charter 2026-10-30: hmp_tongue, rm_buccal, rm_throat, rm_baboon). Within-subject exposure,
# exact nulls (40 reps) and implanted signal (10 reps), on each site and on two stress templates fixed before any exposure
# is drawn: mix-A-B (one sample from each of two sites per subject: pairs that differ strongly) and thin-hmp_tongue (one
# member of each pair thinned to 10% depth). A local smoke run on a fake second site reproduced the failure (c13: 17 false
# calls per null cell on mixed pairs, 2.5 on thinned) and removed it with the observed statistic computed by the
# permutations' own re-centring (c13lin: 0). Candidates (benchmarks/dev/candidates/83_pairdiag.R): c13, c13lin, c13lf
# (+ label-free kap), c13fx (+ rho 1), c13ns (no kernel choice), c13g0 (no centre: exact reference), c12 (clustered path),
# and the subject-aware comparators.
# hpc/guard.sh: shared memory cap (MEMTOTAL), OOM-killer first choice, single-threaded BLAS. Resumable; pushes.
#   git pull && mkdir -p logs && (MEMTOTAL=90G CORES=12 nohup nice -n 10 bash hpc/run_srv32.sh > logs/run_srv32.log 2>&1 &)
source "$(dirname "${BASH_SOURCE[0]}")/guard.sh" "$@"
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
CORES="${CORES:-12}"
L=srv32; R=benchmarks/dev/results/$L
C=eperm_c12,eperm_c13,eperm_c13lin,eperm_c13lf,eperm_c13fx,eperm_c13ns,eperm_c13g0,ref_linda_re,ref_zicoseq_st,ref_ldm_cl
TP=hmp_tongue,rm_buccal,rm_throat,rm_baboon,mix-hmp_tongue-rm_buccal,mix-hmp_tongue-rm_throat,mix-rm_buccal-rm_throat,thin-hmp_tongue
echo ">> rm suite $L ($(date '+%F %T'))"
Rscript benchmarks/dev/rm_suite.R --candidates "$C" --settings P_pair.null,P_pair --template "$TP" --pick-random rm_baboon \
  --reps 10 --null-reps 40 --cores "$CORES" --label "$L" || exit 1
cp logs/run_srv32.log "$R/run.log" 2>/dev/null || true
git add "$R" && git commit -q -m "rm suite srv32: c13's within-subject path, exactness diagnostics on development data" \
  && git pull -q --rebase --autostash && git push -q
echo ">> done $(date '+%F %T')"
