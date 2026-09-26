#!/usr/bin/env bash
# run_gate.sh -- the promotion path in one command (replaces the one-off run_gate_p0x.sh scripts):
#   1. the whole dev suite on the server for CANDS -- house, implant, stress suite, msq, mid (label L),
#      committed and pushed;
#   2. if GATE passes the calibration gate there, the full benchmark (axes A-C) for the benchmark
#      methods M under tag TAG, via run_p03.sh, which pushes results/summary/ when done.
# usage:
#   L=srv6 CANDS=erdl_u_ad,erd_u_ad GATE=erdl_u_ad M=pursue03_erdlu_ad,pursue03_erdu_ad TAG=p06 \
#     bash hpc/run_gate.sh > logs/run_gate_p06.log 2>&1
set -o pipefail
: "${L:?}" "${CANDS:?}" "${GATE:?}" "${M:?}" "${TAG:?}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
R=benchmarks/dev/results/$L; GLOG="logs/run_gate_$TAG.log"
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq,mid --reps 5 --cores 24 --label "$L" --candidates "$CANDS" || exit 1
Rscript benchmarks/dev/compare.R --labels "$L" --brief | tee "$R/gate.txt"
cp "$GLOG" "$R/run.log" 2>/dev/null || true
git add "$R" && git commit -m "dev suite $L ($CANDS), full suite on the server" && git pull --rebase --autostash && git push
if grep -qE "^ *$GATE +PASS" "$R/gate.txt"; then
  echo ">> $GATE passed the whole dev suite -- launching the full benchmark (tag $TAG)"
  M="$M" TAG="$TAG" bash hpc/run_p03.sh > "logs/run_$TAG.log" 2>&1
  rc=$?; if [ $rc -eq 0 ]; then echo ">> full benchmark finished and pushed: see logs/run_$TAG.log"
  else echo ">> full benchmark FAILED (exit $rc): see logs/run_$TAG.log"; fi
else
  echo ">> $GATE did not pass -- full benchmark NOT launched"
fi
