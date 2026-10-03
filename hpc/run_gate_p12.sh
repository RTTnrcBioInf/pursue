#!/usr/bin/env bash
# run_gate_p12.sh -- it39/it40: the compositional centre re-estimated in every permutation (eperm_c), and the scale
# chosen by BH discoveries among the other taxa (eperm_cs).
# 1. the whole dev suite on the server -- house, implant, stress suite, msq, mid (label srv15, pushed);
# 2. ranks them with srv10-14;
# 3. runs the full benchmark (axes A-C, tag p12) for the passing candidate with the most TP per cell.
#   git pull && mkdir -p logs && (nohup bash hpc/run_gate_p12.sh > logs/run_gate_p12.log 2>&1 &)
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
L=srv15; R=benchmarks/dev/results/$L
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq,mid --reps 5 --cores "${CORES:-24}" --label "$L" \
  --candidates eperm_c,eperm_cs || exit 1
Rscript benchmarks/dev/compare.R --labels "$L" --brief | tee "$R/gate.txt"
PREV=$(ls -d benchmarks/dev/results/srv1[0-4] 2>/dev/null | xargs -n1 basename | paste -sd, -)
[ -n "$PREV" ] && Rscript benchmarks/dev/compare.R --labels "$PREV,$L" --brief > "$R/compare_with_srv10-14.txt" 2>&1
cp logs/run_gate_p12.log "$R/run.log" 2>/dev/null || true
git add "$R" && git commit -m "dev suite $L: centre re-estimated per permutation (eperm_c, eperm_cs), full suite on the server" && git pull --rebase --autostash && git push
# null calibration at scale: the dev suite's 10 cells per null setting cannot separate P(any rejection) 0.05 from
# 0.10 (p11: eperm_r 0.100 on axis B B_null, 0.092 on msq R06); 80 cells per setting here, all four PURSUE builds
LN=srv15n; RN=benchmarks/dev/results/$LN
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq,mid --settings house:R06,implant:B_null,msq:R06,mid:R06 \
  --reps 40 --cores "${CORES:-24}" --label "$LN" --candidates efull_b3w,eperm_r,eperm_c,eperm_cs \
  && Rscript -e 'd <- read.csv("benchmarks/dev/results/srv15n/cells.csv"); a <- aggregate(cbind(any = d$fp > 0, fp = d$fp) ~ setting + candidate, d, mean); print(a, row.names = FALSE); print(aggregate(cbind(any = d$fp > 0) ~ candidate, d, mean), row.names = FALSE)' > "$RN/null_any.txt" 2>&1
git add "$RN" && git commit -m "dev suite $LN: null calibration at scale (80 cells per null setting)" && git pull --rebase --autostash && git push
BEST=$(awk '$2 == "PASS" && ($1 == "eperm_c" || $1 == "eperm_cs") {print $8, $1}' "$R/gate.txt" | sort -nr | head -1 | awk '{print $2}')
case "$BEST" in eperm_c) M=pursue03_epermc ;; eperm_cs) M=pursue03_epermcs ;; *) M="" ;; esac
if [ -n "$M" ]; then
  echo ">> $BEST passed with the most TP -- launching the full benchmark (tag p12): $M"
  PURSUE_FOCUS="pursue03_efullb3w,pursue03_epermr,$M" M="$M" TAG=p12 bash hpc/run_p03.sh > logs/run_p12.log 2>&1
  rc=$?; if [ $rc -eq 0 ]; then echo ">> full benchmark finished and pushed: see logs/run_p12.log"
  else echo ">> full benchmark FAILED (exit $rc): see logs/run_p12.log"; fi
else
  echo ">> neither candidate passed the gate -- full benchmark NOT launched"
fi
