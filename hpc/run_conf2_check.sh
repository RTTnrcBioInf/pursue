#!/usr/bin/env bash
# run_conf2_check.sh -- pre-launch check of conf2's data (protocol 13, 2026-10-21): every new template loads through the
# engine (benchmarks/dev/tools/conf2_check.R: dimensions, depths, raw counts, contrast sizes, GEMS1 shared taxa) and the
# repeated-measures templates have the subjects their cells need (rm_suite.R --dry). Format and design only: no method
# runs, no taxon is compared. Report pushed as benchmarks/dev/results/conf2_check/report.txt. A few minutes, < 10 GB.
#   git pull && mkdir -p logs && (nohup bash hpc/run_conf2_check.sh > logs/run_conf2_check.log 2>&1 &)
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"
export PURSUE_BENCH_ROOT="$ROOT/benchmarks" PURSUE_DATA_ROOT="${PURSUE_DATA_ROOT:-$ROOT/benchmarks/data}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
R=benchmarks/dev/results/conf2_check; mkdir -p "$R" logs
RMT=hmp_stool,hmp_gingiva,hmp_skin_ear,rm_buccal,rm_throat,rm_baboon     # protocol 13: repeated measures, held-out data
{ Rscript benchmarks/dev/tools/conf2_check.R; a=$?
  echo; echo "== repeated-measures templates (rm_suite.R --dry; cells need 40 subjects with two samples)"
  Rscript benchmarks/dev/rm_suite.R --dry --settings P_subj,P_pair --template "$RMT" --pick-random rm_baboon --candidates eperm_c13; b=$?
  echo; echo "commit $(git rev-parse --short HEAD); data check exit $a, rm check exit $b"; } 2>&1 | tee "$R/report.txt"
git add "$R" && git commit -q -m "conf2 pre-launch data check (format and design only)" && git pull -q --rebase --autostash && git push -q
echo ">> done $(date '+%F %T')"
