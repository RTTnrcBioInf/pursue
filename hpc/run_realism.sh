#!/usr/bin/env bash
# run_realism.sh -- the protocol's realism gate (section 4.4) for all five simulators on the five evaluation templates and
# the two tuning templates, with descriptive statistics on why a simulator is separable (benchmarks/dev/tools/realism_all.R).
# It was planned as a SLURM array (hpc/slurm/realism.sbatch) that hpc/run_local.sh cannot run, so it was never computed
# here. sparseDOSSA2 reuses the cached fits in cache/. Light (35 tasks); safe beside another run if MEMTOTAL matches it.
# Resumable; pushes the summary only (benchmarks/dev/results/realism/).
#   git pull && mkdir -p logs && (MEMTOTAL=130G CORES=6 nohup nice -n 10 bash hpc/run_realism.sh > logs/run_realism.log 2>&1 &)
source "$(dirname "${BASH_SOURCE[0]}")/guard.sh" "$@"
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs
CORES="${CORES:-6}"
R=benchmarks/dev/results/realism
echo ">> realism gate ($(date '+%F %T'))"
Rscript benchmarks/dev/tools/realism_all.R --cores "$CORES" || exit 1
cp logs/run_realism.log "$R/run.log" 2>/dev/null || true
git add "$R" && git commit -q -m "realism gate: five simulators x seven templates, with descriptive statistics" \
  && git pull -q --rebase --autostash && git push -q
echo ">> done $(date '+%F %T')"
