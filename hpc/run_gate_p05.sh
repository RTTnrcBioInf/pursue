#!/usr/bin/env bash
# run_gate_p05.sh -- it13: pairwise U-statistics with both pair members thinned (rho 0.7 and 0.5).
# 1. the whole dev suite on the server -- house, implant, stress suite, msq, mid (label srv5, pushed);
# 2. picks the larger rho whose erdl_u variant passes the calibration gate on ALL of it;
# 3. only then runs the full benchmark (axes A-C, tag p05) for that rho, which pushes results/summary/.
#   git pull && mkdir -p logs && (nohup bash hpc/run_gate_p05.sh > logs/run_gate_p05.log 2>&1 &)
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
L=srv5; R=benchmarks/dev/results/$L
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq,mid --reps 5 --cores 24 --label "$L" \
  --candidates erdl_u_r07,erd_u_r07,erdl_u_r05,erd_u_r05 || exit 1
Rscript benchmarks/dev/compare.R --labels "$L" --brief | tee "$R/gate.txt"
cp logs/run_gate_p05.log "$R/run.log" 2>/dev/null || true
git add "$R" && git commit -m "dev suite $L: pairwise U-statistics, rho 0.7 / 0.5, full suite on the server" && git pull --rebase --autostash && git push
TAGR=""
if grep -qE '^ *erdl_u_r07 +PASS' "$R/gate.txt"; then TAGR=r07; elif grep -qE '^ *erdl_u_r05 +PASS' "$R/gate.txt"; then TAGR=r05; fi
if [ -n "$TAGR" ]; then
  echo ">> erdl_u_$TAGR passed everywhere -- launching the full benchmark (tag p05)"
  M="pursue03_erdlu_$TAGR,pursue03_erdu_$TAGR" TAG=p05 bash hpc/run_p03.sh > logs/run_p05.log 2>&1
  echo ">> full benchmark finished: see logs/run_p05.log"
else
  echo ">> neither rho passed -- full benchmark NOT launched"
fi
