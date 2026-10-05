#!/usr/bin/env bash
# run_DE.sh -- the confirmation run on the real-data axes (protocol sections 8-9): D, biological and
# experimental truth (gingival aerobes, BV taxa, Stammler spike-ins); E, replicability over random
# 50/50 splits (crc_genus, risk_ileum). Every comparator plus the PURSUE 0.3 lead.
#
#   P=24 nohup bash hpc/run_DE.sh > logs/run_DE.log 2>&1 &
#   SEED: master seed for the random splits and spike-in groupings (default 1; the confirmation passes its own)
#
# One process per (task, method), P at a time, one BLAS thread each. Splits and spike-in groupings
# are drawn before any method runs, so every method sees the same ones. Resumable: a (task, method)
# whose CSV exists is skipped. A smoke pass (wilcoxon_tss on every task) runs first and stops the
# run if any task fails. Merged tables and a summary go to results/summary/axisDE/ and are pushed.
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs/DE results/axisDE
export PURSUE_BENCH_ROOT="$ROOT/benchmarks" PURSUE_DATA_ROOT="${PURSUE_DATA_ROOT:-$ROOT/benchmarks/data}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
P="${P:-24}"; ONLY="${ONLY:-.}"; LEAD="${LEAD:-pursue03_epermcs1}"; SEED="${SEED:-1}"
M="${M:-wilcoxon_tss,lm_logtss,limma_logtss,logistic_presence,pursue,linda,ancombc2,maaslin3,aldex2,corncob,ldm,locom,locom2,zicoseq,fastancom,adapt,fastemu,$LEAD}"
R=benchmarks/R/engine/run_axisDE.R; O=results/axisDE
SPK="$(tr -d '[:space:]' < benchmarks/expected/stammler_spikein_ids.txt)"
# name | output-file stem (method appended) | arguments
TASKS=(
  "gingival|D__mbd_gingival_v35|--axis D --template mbd_gingival_v35 --group body_subsite --case ^supra --expected expected/gingival_aerobes.tsv|biotruth"
  "bv|D__mbd_ravel_bv|--axis D --template mbd_ravel_bv --group study_condition --case vagin|bv --expected expected/bv_taxa.tsv|biotruth"
  "spikein|D__mbd_stammler_spikein|--axis D --template mbd_stammler_spikein --spikein $SPK|spikein"
  "crc|E__crc_genus|--axis E --template crc_genus --group diagnosis --levels control,CRC --splits 5|replicability"
  "risk|E__risk_ileum|--axis E --template risk_ileum --group diagnosis --levels no,CD --splits 5|replicability"
)
SEL=(); for t in "${TASKS[@]}"; do [[ "${t%%|*}" =~ $ONLY ]] && SEL+=("$t"); done; TASKS=("${SEL[@]}")
jobline() {  # $1 task spec, $2 method -> one shell command (skips when its CSV exists)
  # fields: name|stem|args|kind -- args may contain '|' (the --case regex), so split from both ends
  name="${1%%|*}"; rest="${1#*|}"; stem="${rest%%|*}"; rest="${rest#*|}"; kind="${rest##*|}"; args="${rest%|*}"
  local out="$O/${stem}__$2__${kind}.csv"
  local a; a="$(sed -E "s/--case ([^ ]+)/--case '\1'/" <<<"$args")"
  echo "[ -f $out ] || Rscript $R $a --methods $2 --tag $2 --out $O --master-seed $SEED > logs/DE/${name}__$2.log 2>&1; [ -f $out ] && echo \"ok   $name $2\" || echo \"FAIL $name $2\""
}

echo ">> smoke: wilcoxon_tss on every task"; fail=0
for t in "${TASKS[@]}"; do bash -c "$(jobline "$t" wilcoxon_tss)" | tee -a logs/DE/smoke.txt | grep -q '^ok' || fail=1; done
if [ $fail -ne 0 ]; then echo ">> smoke FAILED (see logs/DE/*__wilcoxon_tss.log) -- not launching"; exit 1; fi
grep -h ">> contrast" logs/DE/*__wilcoxon_tss.log

echo ">> $(echo "$M" | tr ',' '\n' | grep -c .) methods x ${#TASKS[@]} tasks, $P at a time"
JOBS=$(mktemp)
for m in $(echo "$M" | tr ',' ' '); do for t in "${TASKS[@]}"; do jobline "$t" "$m"; done; done > "$JOBS"
xargs -d '\n' -P "$P" -I CMD bash -c CMD < "$JOBS" | tee logs/DE/status.txt
rm -f "$JOBS"

echo ">> merge"
Rscript hpc/merge_DE.R "$O" results/summary/axisDE | tee logs/DE/summary.txt
cp logs/run_DE.log results/summary/axisDE/run_DE.log 2>/dev/null || true
git add results/summary/axisDE && git commit -m "axes D/E confirmation run: comparators + $LEAD" \
  && git pull --rebase --autostash && git push
echo ">> done: $(grep -c '^ok' logs/DE/status.txt) ok, $(grep -c '^FAIL' logs/DE/status.txt) failed"
