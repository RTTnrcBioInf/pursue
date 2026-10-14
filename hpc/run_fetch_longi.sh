#!/usr/bin/env bash
# run_fetch_longi.sh -- download public 16S datasets with repeated samples per subject (HMP16SData V35 sites, HMP2Data
# cohorts, the wild-baboon time series) into benchmarks/data/longitudinal/ (not committed) and push only the summary
# of what arrived (benchmarks/dev/results/longi_fetch/summary.txt). Installs the three Bioconductor data packages into
# the current R library if they are missing. Network and disk only; no analysis.
#   git pull && mkdir -p logs && (MEMTOTAL=40G nohup nice -n 10 bash hpc/run_fetch_longi.sh > logs/run_fetch_longi.log 2>&1 &)
source "$(dirname "${BASH_SOURCE[0]}")/guard.sh" "$@"
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs
R=benchmarks/dev/results/longi_fetch
Rscript benchmarks/dev/tools/fetch_longitudinal.R
du -sh benchmarks/data/longitudinal 2>/dev/null | tee -a "$R/summary.txt"
cp logs/run_fetch_longi.log "$R/run.log" 2>/dev/null || true
git add "$R" && git commit -q -m "longitudinal 16S datasets: what was fetched (summary only; data stay on the server)" \
  && git pull -q --rebase --autostash && git push -q
echo ">> done $(date '+%F %T')"
