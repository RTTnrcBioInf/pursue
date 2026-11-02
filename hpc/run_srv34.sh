#!/usr/bin/env bash
# run_srv34.sh -- sign diagnostics for c14's reported effect (benchmarks/dev/tools/sign_diag.R; no method change). For every
# taxon c14 calls that truly changes: true direction, sign of the reported estimate, the kernel that made the call and the
# sign of its statistic; LinDA on the same cells. Tuning templates, dev seed: house R00 / R11 / R13 / R16, implants B_ref,
# B_prev, B_abund, B_allup, and two new implants with presence and abundance moving in opposite directions (B_opp, B_mix4).
# 10 settings x 2 templates x 20 reps = 400 cells (~35 s each for c14). Pushes summary.txt and taxa.csv (cells/ dropped).
# hpc/guard.sh: shared memory cap (MEMTOTAL), OOM-killer first choice, single-threaded BLAS. Resumable until it pushes.
#   git pull && mkdir -p logs && (MEMTOTAL=90G CORES=12 nohup nice -n 10 bash hpc/run_srv34.sh > logs/run_srv34.log 2>&1 &)
source "$(dirname "${BASH_SOURCE[0]}")/guard.sh" "$@"
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
CORES="${CORES:-12}"; L=srv34; R=benchmarks/dev/results/$L
echo ">> sign diagnostics $L ($(date '+%F %T'))"
Rscript benchmarks/dev/tools/sign_diag.R --reps 20 --cores "$CORES" --label "$L" || exit 1
rm -rf "$R/cells"; cp logs/run_srv34.log "$R/run.log" 2>/dev/null || true
git add "$R" && git commit -q -m "srv34: sign diagnostics for c14's reported effect (development data)" \
  && git pull -q --rebase --autostash && git push -q
echo ">> done $(date '+%F %T')"
