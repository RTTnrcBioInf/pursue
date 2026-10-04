#!/usr/bin/env bash
# run_gate_p13.sh -- it41: selection margin 1 (eperm_cs1) and the unthinned lin / sqp scales (eperm_x1, eperm_x2).
# 1. the whole dev suite on the server (label srv16, pushed), ranked with srv10-15;
# 2. null calibration at scale (srv16n: house R06, implant B_null, msq R06, mid R06 x 2 templates x 40 reps, pushed);
# 3. the full benchmark (axes A-C, tag p13) for the candidate with the most TP per cell among those that pass the
#    dev-suite gate AND keep P(any rejection) <= 0.075 over the 320 null cells.
#   git pull && mkdir -p logs && (nohup bash hpc/run_gate_p13.sh > logs/run_gate_p13.log 2>&1 &)
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
C=eperm_cs1,eperm_x1,eperm_x2
L=srv16; R=benchmarks/dev/results/$L
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq,mid --reps 5 --cores "${CORES:-24}" --label "$L" --candidates "$C" || exit 1
Rscript benchmarks/dev/compare.R --labels "$L" --brief | tee "$R/gate.txt"
PREV=$(ls -d benchmarks/dev/results/srv1[0-5] 2>/dev/null | xargs -n1 basename | paste -sd, -)
[ -n "$PREV" ] && Rscript benchmarks/dev/compare.R --labels "$PREV,$L" --brief > "$R/compare_with_srv10-15.txt" 2>&1
cp logs/run_gate_p13.log "$R/run.log" 2>/dev/null || true
git add "$R" && git commit -m "dev suite $L: margin 1 and unthinned scales (eperm_cs1, eperm_x1, eperm_x2)" && git pull --rebase --autostash && git push
LN=srv16n; RN=benchmarks/dev/results/$LN
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq,mid --settings house:R06,implant:B_null,msq:R06,mid:R06 \
  --reps 40 --cores "${CORES:-24}" --label "$LN" --candidates "$C" \
  && Rscript -e 'd <- read.csv("benchmarks/dev/results/srv16n/cells.csv"); print(aggregate(cbind(any = d$fp > 0, fp = d$fp) ~ setting + candidate, d, mean), row.names = FALSE); a <- aggregate(cbind(any = d$fp > 0) ~ candidate, d, mean); print(a, row.names = FALSE); writeLines(a$candidate[a$any <= 0.075], "benchmarks/dev/results/srv16n/null_ok.txt")' > "$RN/null_any.txt" 2>&1
git add "$RN" && git commit -m "dev suite $LN: null calibration at scale for it41" && git pull --rebase --autostash && git push
BEST=$(awk '$2 == "PASS" {print $8, $1}' "$R/gate.txt" | grep -wFf "$RN/null_ok.txt" | sort -nr | head -1 | awk '{print $2}')
case "$BEST" in eperm_cs1) M=pursue03_epermcs1 ;; eperm_x1) M=pursue03_epermx1 ;; eperm_x2) M=pursue03_epermx2 ;; *) M="" ;; esac
if [ -n "$M" ]; then
  echo ">> $BEST passed the gate and the null check with the most TP -- launching the full benchmark (tag p13): $M"
  PURSUE_FOCUS="pursue03_efullb3w,pursue03_epermcs,$M" M="$M" TAG=p13 bash hpc/run_p03.sh > logs/run_p13.log 2>&1
  rc=$?; if [ $rc -eq 0 ]; then echo ">> full benchmark finished and pushed: see logs/run_p13.log"
  else echo ">> full benchmark FAILED (exit $rc): see logs/run_p13.log"; fi
else
  echo ">> no candidate passed both the gate and the null check -- full benchmark NOT launched"
fi
