#!/usr/bin/env bash
# run_conf2_qmpfix.sh -- conf2 deviation (protocol 13, 2026-10-28). Axis D's measured-load scoring (run_axisDE.R --qmp)
# stopped with an error for the 7 methods that report no effect estimate in some arm (pursue, ancombc2, maaslin3, corncob,
# ldm, locom2, zicoseq): the Spearman correlation of an all-NA vector ("no complete element pairs"). It now returns NA
# there; nothing else changed. This re-runs only those 7 Vandeputte jobs (the other 11 methods' outputs exist and are
# kept), same master seed, with every other file as tagged conf2-code; re-merges the D/E summary and pushes, with its log.
#   git pull && mkdir -p logs && (nohup bash hpc/run_conf2_qmpfix.sh > logs/run_conf2_qmpfix.log 2>&1 &)
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs
FROZEN="benchmarks/R benchmarks/dev/candidates benchmarks/dev/tools benchmarks/dev/rm_suite.R benchmarks/templates.tsv benchmarks/expected hpc"
ch=$(git diff --name-only conf2-code -- $FROZEN | grep -v '^hpc/run_conf2_qmpfix.sh$')
echo ">> code differing from conf2-code: ${ch:-none}"
[ "$ch" = "benchmarks/R/engine/run_axisDE.R" ] || { echo ">> expected run_axisDE.R only -- stopping"; exit 1; }
git diff conf2-code -- benchmarks/R/engine/run_axisDE.R
SET=conf2 ONLY='^vandeputte$' SEED=20261021 LEAD=pursue03_epermc13 P=7 OUT=results_conf2/axisDE \
  SUMMARY=results_conf2/summary/axisDE DE_LOG=logs/run_conf2_DE.log bash hpc/run_DE.sh
cp logs/run_conf2_qmpfix.log results_conf2/summary/axisDE/run_qmpfix.log 2>/dev/null || true
git add results_conf2/summary/axisDE && git commit -q -m "conf2 deviation: Vandeputte measured-load scoring for the 7 methods without estimates (protocol 13, 2026-10-28)" \
  && git pull -q --rebase --autostash && git push -q
echo ">> done $(date '+%F %T')"
