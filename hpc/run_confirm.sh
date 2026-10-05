#!/usr/bin/env bash
# run_confirm.sh -- the confirmation run of a PURSUE release candidate (rd-charter: one run, never used for
# choosing; protocol sections 8-9 and 11).
#   1. axes D/E -- real data, never run before: every comparator plus the candidate (hpc/run_DE.sh);
#   2. axes A-C under a FRESH master seed: evaluation pool, five simulators, the candidate, PURSUE 0.2 and
#      every comparator that was competitive on seed 1 (calibrated power within reach, or a mandatory
#      baseline). Written to results_conf1/ with tag conf1, so it neither skips nor mixes with the seed-1
#      cells; aggregated into results_conf1/summary/ and pushed.
# The heavy, dominated comparators (LOCOM, LOCOM2, ANCOM-BC2, ALDEx2, corncob, MaAsLin 3, fastEmu: ~70% of
# the CPU of a full re-run, none within 10% of the candidate's calibrated power on seed 1) are left out of
# A-C by default; M=... overrides the list. MAXJOBS defaults to 40 so dev-suite runs can share the machine.
#
#   git pull && mkdir -p logs && (nohup bash hpc/run_confirm.sh > logs/run_confirm.log 2>&1 &)
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs
export PURSUE_BENCH_ROOT="$ROOT/benchmarks" PURSUE_DATA_ROOT="$ROOT/benchmarks/data"
SEED="${SEED:-20261005}"; LEAD="${LEAD:-pursue03_epermcs1}"; TAG=conf1; RR=results_conf1
M="${M:-$LEAD,pursue,zicoseq,adapt,ldm,linda,wilcoxon_tss,logistic_presence,limma_logtss,lm_logtss,fastancom}"
echo ">> confirmation of $LEAD -- master seed $SEED -- started $(date '+%F %T')"

if [ -z "${SKIP_DE:-}" ]; then
  echo ">> 1. axes D/E (all comparators + $LEAD; log: logs/run_DE.log)"
  SEED="$SEED" LEAD="$LEAD" P="${P_DE:-24}" bash hpc/run_DE.sh > logs/run_DE.log 2>&1
  tail -2 logs/run_DE.log
fi

echo ">> 2. axes A-C, master seed $SEED, methods: $M"
SM=$(mktemp -d)
Rscript benchmarks/R/engine/run_cell.R --axis A --simulator house --template hmp_stool --regime R00 --replicate 1 \
  --methods "$M" --master-seed "$SEED" --out "$SM" --tag smoke 2>&1 | tee "$SM/smoke.log"
bad=0; for m in $(echo "$M" | tr ',' ' '); do l=$(grep -E "^  $m " "$SM/smoke.log")
  if [ -z "$l" ] || echo "$l" | grep -qE 'error|not_installed'; then echo ">> smoke: $m did not run (${l:-no line})"; bad=1; fi; done
if [ $bad -ne 0 ]; then echo ">> smoke FAILED -- not launching A-C"; exit 1; fi
rm -rf "$SM"
Rscript hpc/make_tasklist.R --pool evaluation --simulators house,msq,mid,sd2,sps --methods "$M" --out hpc/tasks_conf1
for ax in A B C; do
  echo ">> axis $ax"
  PURSUE_MASTER_SEED="$SEED" TAG="$TAG" RESULTS_ROOT="$RR" MAXJOBS="${MAXJOBS:-40}" bash hpc/run_adaptive.sh "hpc/tasks_conf1/axis$ax.txt"
done
echo ">> aggregate"
PURSUE_FOCUS="$LEAD,pursue" Rscript benchmarks/R/analysis/aggregate.R "$RR"
echo "master_seed=$SEED lead=$LEAD methods=$M finished=$(date '+%F %T')" > "$RR/summary/CONFIRMATION.txt"
cp logs/run_confirm.log "$RR/summary/run_confirm.log" 2>/dev/null || true
git add .gitignore "$RR/summary" && git commit -m "confirmation run (axes A-C, fresh master seed $SEED): $LEAD" \
  && git pull --rebase --autostash && git push
echo ">> done $(date '+%F %T')"
