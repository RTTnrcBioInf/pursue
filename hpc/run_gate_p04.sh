#!/usr/bin/env bash
# run_gate_p04.sh -- it12's pairwise U-statistic candidates: the msq/mid half of the dev suite first
# (label srv4, pushed), and ONLY if erdl_u passes the calibration gate there, the full benchmark
# (axes A-C, tag p04, via run_p03.sh, which pushes results/summary/ when done).
#   git pull && mkdir -p logs && (nohup bash hpc/run_gate_p04.sh > logs/run_gate_p04.log 2>&1 &)
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
L=srv4; R=benchmarks/dev/results/$L
Rscript benchmarks/dev/devsuite.R --sims msq,mid --reps 5 --cores 16 --label "$L" --candidates erdl_u,erd_u2,erl_u || exit 1
Rscript benchmarks/dev/compare.R --labels "$L" --brief | tee "$R/gate.txt"
cp logs/run_gate_p04.log "$R/run.log" 2>/dev/null || true
git add "$R" && git commit -m "dev suite $L: pairwise U-statistics on msq/mid" && git pull --rebase --autostash && git push
if grep -qE '^ *erdl_u +PASS' "$R/gate.txt"; then
  echo ">> erdl_u passed on msq/mid -- launching the full benchmark (tag p04)"
  M=pursue03_erdlu,pursue03_erdu,pursue03_erlu TAG=p04 bash hpc/run_p03.sh > logs/run_p04.log 2>&1
  rc=$?; if [ $rc -eq 0 ]; then echo ">> full benchmark finished and pushed: see logs/run_p04.log"
  else echo ">> full benchmark FAILED (exit $rc): see logs/run_p04.log"; fi
else
  echo ">> erdl_u did not pass on msq/mid -- full benchmark NOT launched"
fi
