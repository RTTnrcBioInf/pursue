#!/usr/bin/env bash
# run_srv22.sh -- it49, beside p14 (8 cores; p14 holds 56): designed depth confounding on the efull_nodet path.
#   srv22:  eperm_c6a (= c5q, log alone under designed depth confounding) and eperm_c6b (+ pairs thinned to 0.35) on
#           house / msq / mid R17 and R19, 10 reps -- the same cells as srv21, so c5q there is the reference;
#   srv22n: mid R17 / R19, msq R19, house R19 nulls at 40 reps, with eperm_c5q alongside.
# Everything else in c6a / c6b is c5q, so only these settings can differ. Resumable; pushes its results.
#   git pull && mkdir -p logs && (CORES=8 nohup bash hpc/run_srv22.sh > logs/run_srv22.log 2>&1 &)
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
CORES="${CORES:-8}"
L=srv22; R=benchmarks/dev/results/$L
echo ">> dev suite $L ($(date '+%F %T'))"
Rscript benchmarks/dev/devsuite.R --suite full --sims house,msq,mid --settings house:R17,house:R19,msq:R17,msq:R19,mid:R17,mid:R19 \
  --reps 10 --cores "$CORES" --label "$L" --candidates eperm_c6a,eperm_c6b || exit 1
Rscript benchmarks/dev/compare.R --labels "srv21,$L" > "$R/compare.txt" 2>&1
LN=srv22n; RN=benchmarks/dev/results/$LN
echo ">> null run $LN ($(date '+%F %T'))"
Rscript benchmarks/dev/devsuite.R --suite full --sims house,msq,mid --settings none --extra mid:R17.null,mid:R19.null,msq:R19.null,house:R19.null \
  --reps 40 --cores "$CORES" --label "$LN" --candidates eperm_c5q,eperm_c6a,eperm_c6b \
  && Rscript benchmarks/dev/tools/null_summary.R "$RN/cells.csv" > "$RN/null_summary.txt" 2>&1
cp logs/run_srv22.log "$R/run.log" 2>/dev/null || true
git add "$R" "$RN" && git commit -m "dev suite $L / $LN: eperm_c6a / c6b (log alone, harder thinning under designed depth confounding)" && git pull --rebase --autostash && git push
echo ">> done $(date '+%F %T')"
