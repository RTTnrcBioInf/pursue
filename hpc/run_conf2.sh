#!/usr/bin/env bash
# run_conf2.sh -- the second confirmation run (protocol section 13, pre-registered 2026-10-21; rd-charter: one run, never
# used for choosing). Candidate pursue03_epermc13 (eperm_c12 plus the exact within-subject path), master seed 20261021.
#   0. freeze -- tags the code conf2-code (pushed) on first launch; before every stage, checks that the frozen paths
#                (benchmarks/R, the dev candidates and tools, rm_suite.R, templates.tsv, expected/, hpc/) are as tagged,
#                and stops if not;
#   1. DE     -- hpc/run_DE.sh SET=conf2: conf1's real-data tasks plus Vandeputte 2017 (truth from measured load), six
#                replicability datasets, GEMS1's four sites and their six cross-site pairs; every comparator and the lead;
#   2. RM     -- repeated measures on held-out real data (hmp_stool, hmp_gingiva, hmp_skin_ear, rm_buccal, rm_throat,
#                rm_baboon): subject-level and within-subject exposure, signal (20 reps) and exact nulls (40 reps);
#   3. AC     -- axes A-C, the five simulators on the evaluation pool under the fresh seed, as conf1 (results_conf2/, tag conf2);
#   4. results_conf2/summary/CONFIRMATION.txt (seed, lead, tagged commit, freeze status) and push.
# Resumable: finished jobs and cells are kept, so relaunching continues. STAGES="RM AC" runs a subset. Each stage pushes.
# hpc/guard.sh: shared memory cap (MEMTOTAL), OOM-killer first choice, single-threaded BLAS.
#   git pull && bash hpc/run_conf2_check.sh                 # first: the data check (minutes)
#   git pull && mkdir -p logs && (MEMTOTAL=200G CORES=40 nohup nice -n 10 bash hpc/run_conf2.sh > logs/run_conf2.log 2>&1 &)
source "$(dirname "${BASH_SOURCE[0]}")/guard.sh" "$@"
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs
export PURSUE_BENCH_ROOT="$ROOT/benchmarks" PURSUE_DATA_ROOT="$ROOT/benchmarks/data"
SEED=20261021; LEAD=pursue03_epermc13; TAG=conf2; RR=results_conf2; CODE=conf2-code
CORES="${CORES:-40}"; STAGES="${STAGES:-DE RM AC}"
M="${M:-$LEAD,pursue,zicoseq,adapt,ldm,linda,wilcoxon_tss,logistic_presence,limma_logtss,lm_logtss,fastancom}"
RC=eperm_c13,ref_linda,ref_linda_re,ref_ldm,ref_ldm_cl,ref_zicoseq,ref_zicoseq_st,ref_adapt,ref_logistic_presence,ref_logistic_cr
RMT=hmp_stool,hmp_gingiva,hmp_skin_ear,rm_buccal,rm_throat,rm_baboon
FROZEN="benchmarks/R benchmarks/dev/candidates benchmarks/dev/tools benchmarks/dev/rm_suite.R benchmarks/templates.tsv benchmarks/expected hpc"
has() { [[ " $STAGES " == *" $1 "* ]]; }
push() { git add "$@" && git commit -q -m "$MSG" && git pull -q --rebase --autostash && git push -q; }
echo ">> conf2: $LEAD, master seed $SEED, stages $STAGES -- started $(date '+%F %T')"

# 0. freeze
git rev-parse -q --verify "refs/tags/$CODE" >/dev/null || git fetch -q origin "refs/tags/$CODE:refs/tags/$CODE" 2>/dev/null
if git rev-parse -q --verify "refs/tags/$CODE" >/dev/null; then
  echo ">> code: $CODE at $(git rev-parse --short "$CODE^{commit}") (resuming)"
else
  if ! git diff --quiet HEAD -- $FROZEN; then echo ">> uncommitted changes in frozen paths -- commit or discard them first:"; git status --short -- $FROZEN; exit 1; fi
  git tag -a "$CODE" -m "conf2 (protocol 13): code frozen for the confirmation run; lead $LEAD, master seed $SEED" \
    && git push -q origin "$CODE" && echo ">> code: tagged $CODE at $(git rev-parse --short HEAD) and pushed" \
    || { echo ">> could not create or push the tag $CODE"; exit 1; }
fi
frozen() { git diff --quiet "$CODE" -- $FROZEN && return 0
  echo ">> frozen code changed since $CODE:"; git diff --stat "$CODE" -- $FROZEN
  echo ">> stopping: the confirmation runs on the tagged code only (restore it with: git checkout $CODE -- <paths>)"; exit 1; }

# 1. axes D/E
if has DE; then frozen
  echo ">> 1. axes D/E -- every comparator + $LEAD (log: logs/run_conf2_DE.log) $(date '+%F %T')"
  SET=conf2 SEED="$SEED" LEAD="$LEAD" P="${P_DE:-24}" OUT="$RR/axisDE" SUMMARY="$RR/summary/axisDE" DE_LOG=logs/run_conf2_DE.log \
    bash hpc/run_DE.sh > logs/run_conf2_DE.log 2>&1 || { echo ">> D/E stopped (smoke failed?): see logs/run_conf2_DE.log"; tail -5 logs/run_conf2_DE.log; exit 1; }
  tail -1 logs/run_conf2_DE.log
fi

# 2. repeated measures, held-out real data
if has RM; then frozen
  echo ">> 2. repeated measures on $RMT (log: logs/run_conf2_RM.log) $(date '+%F %T')"
  Rscript benchmarks/dev/rm_suite.R --candidates "$RC" --settings P_subj,P_pair,P_subj.null,P_pair.null --template "$RMT" \
    --pick-random rm_baboon --reps 20 --null-reps 40 --cores "${CORES_RM:-$CORES}" --seed "$SEED" --label conf2_rm --out "$RR/rm" \
    > logs/run_conf2_RM.log 2>&1 || { echo ">> RM failed: see logs/run_conf2_RM.log"; tail -5 logs/run_conf2_RM.log; exit 1; }
  mkdir -p "$RR/summary/rm"; cp "$RR"/rm/summary* "$RR/rm/cells.csv" "$RR/summary/rm/"; cp logs/run_conf2_RM.log "$RR/summary/rm/run.log"
  MSG="conf2 repeated measures (held-out real data, seed $SEED): $LEAD (eperm_c13) and comparators" push "$RR/summary/rm"
fi

# 3. axes A-C under the fresh master seed
if has AC; then frozen
  echo ">> 3. axes A-C, master seed $SEED, methods: $M $(date '+%F %T')"
  SM=$(mktemp -d)
  Rscript benchmarks/R/engine/run_cell.R --axis A --simulator house --template hmp_stool --regime R00 --replicate 1 \
    --methods "$M" --master-seed "$SEED" --out "$SM" --tag smoke 2>&1 | tee "$SM/smoke.log"
  bad=0; for m in $(echo "$M" | tr ',' ' '); do l=$(grep -E "^  $m " "$SM/smoke.log")
    if [ -z "$l" ] || echo "$l" | grep -qE 'error|not_installed'; then echo ">> smoke: $m did not run (${l:-no line})"; bad=1; fi; done
  if [ $bad -ne 0 ]; then echo ">> smoke FAILED -- not launching A-C"; exit 1; fi
  rm -rf "$SM"
  Rscript hpc/make_tasklist.R --pool evaluation --simulators house,msq,mid,sd2,sps --methods "$M" --out hpc/tasks_conf2
  for ax in A B C; do frozen; echo ">> axis $ax $(date '+%F %T')"
    PURSUE_MASTER_SEED="$SEED" TAG="$TAG" RESULTS_ROOT="$RR" MAXJOBS="${MAXJOBS:-$CORES}" bash hpc/run_adaptive.sh "hpc/tasks_conf2/axis$ax.txt"
  done
  echo ">> aggregate"; PURSUE_FOCUS="$LEAD,pursue" Rscript benchmarks/R/analysis/aggregate.R "$RR"
fi

# 4. record and push
mkdir -p "$RR/summary"
fz=$(git diff --quiet "$CODE" -- $FROZEN && echo yes || echo NO)
echo "master_seed=$SEED lead=$LEAD code=$CODE@$(git rev-parse --short "$CODE^{commit}") frozen=$fz stages=$STAGES methods_AC=$M rm_candidates=$RC rm_templates=$RMT finished=$(date '+%F %T')" \
  >> "$RR/summary/CONFIRMATION.txt"
cp logs/run_conf2.log "$RR/summary/run_conf2.log" 2>/dev/null || true
MSG="confirmation run conf2 (fresh master seed $SEED, stages $STAGES): $LEAD" push "$RR/summary"
echo ">> done $(date '+%F %T')"
