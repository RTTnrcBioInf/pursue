#!/usr/bin/env bash
# run_gate_p11.sh -- it38: exact permutation p-values in exchangeable designs (eperm, eperm_r).
# 1. the whole dev suite on the server -- house, implant, stress suite, msq, mid (label srv14, pushed);
# 2. ranks them with the srv13 results (efull_b3w and the reference comparators were scored there);
# 3. if a candidate passes the calibration gate everywhere, runs the full benchmark (axes A-C, tag p11)
#    for it -- eperm_r first, it is the superset -- which pushes results/summary/.
#   git pull && mkdir -p logs && (nohup bash hpc/run_gate_p11.sh > logs/run_gate_p11.log 2>&1 &)
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
L=srv14; R=benchmarks/dev/results/$L
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq,mid --reps 5 --cores "${CORES:-24}" --label "$L" \
  --candidates eperm,eperm_r || exit 1
Rscript benchmarks/dev/compare.R --labels "$L" --brief | tee "$R/gate.txt"
PREV=$(ls -d benchmarks/dev/results/srv1[0-3] 2>/dev/null | xargs -n1 basename | paste -sd, -)
[ -n "$PREV" ] && Rscript benchmarks/dev/compare.R --labels "$PREV,$L" --brief > "$R/compare_with_srv10-13.txt" 2>&1
cp logs/run_gate_p11.log "$R/run.log" 2>/dev/null || true
git add "$R" && git commit -m "dev suite $L: exact permutation p-values (eperm, eperm_r), full suite on the server" && git pull --rebase --autostash && git push
M=""                                                       # one candidate (a solo re-run takes ~10 h): eperm_r if it passes, else eperm
if grep -qE '^ *eperm_r +PASS' "$R/gate.txt"; then M="pursue03_epermr"; elif grep -qE '^ *eperm +PASS' "$R/gate.txt"; then M="pursue03_eperm"; fi
if [ -n "$M" ]; then
  echo ">> passed the gate: $M -- launching the full benchmark (tag p11)"
  PURSUE_FOCUS="pursue03_efullb3w,$M" M="$M" TAG=p11 bash hpc/run_p03.sh > logs/run_p11.log 2>&1
  rc=$?; if [ $rc -eq 0 ]; then echo ">> full benchmark finished and pushed: see logs/run_p11.log"
  else echo ">> full benchmark FAILED (exit $rc): see logs/run_p11.log"; fi
else
  echo ">> neither candidate passed the gate -- full benchmark NOT launched"
fi
