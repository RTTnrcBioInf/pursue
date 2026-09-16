#!/usr/bin/env bash
# run_adaptive.sh -- run a benchmark task list with MEMORY-AWARE concurrency.
#
# A fixed -P is wrong in both directions here: axis-A cells range from ~0.1 GB (hmp_vagina) to
# 11.3 GB (sd2/risk_stool), so -P 44 kills the big ones and -P 10 wastes 40 cores on the small
# ones. This schedules against a memory budget instead, and estimates each cell's cost from the
# logs/<cell>.mem files previous runs already wrote -- so it gets more accurate as it goes and
# needs no table to be maintained by hand.
#
#   bash hpc/run_adaptive.sh hpc/tasks/axisC.txt                  # 75% of RAM, all cores
#   MEM_PCT=60 MAXJOBS=32 bash hpc/run_adaptive.sh hpc/tasks/axisC.txt
#   DRY=1 bash hpc/run_adaptive.sh hpc/tasks/axisC.txt            # show the schedule only
#
# Env:
#   MEM_PCT        share of total RAM this run may hold            (default 75)
#   MAXJOBS        hard cap on concurrent cells                    (default nproc)
#   FLOOR_MB       MemAvailable never allowed below this           (default 8192)
#   DEFAULT_EST_MB estimate for a cell class never seen before     (default 4096)
#   SAFETY_PCT     percent of the learned estimate to reserve       (default 125)
#   TAG            re-run tag, as in run_local.sh
#   CELL_CMD       override the per-cell command (testing only)
# NOT `set -u`: in bash before 4.4, ${#arr[@]} on an empty associative array is an unbound-variable
# error, and both the estimate table and the running-PID table are legitimately empty at the start
# and end of the run. Every variable below carries its own ${X:-default} instead.
set -o pipefail
TASKS="${1:?usage: bash hpc/run_adaptive.sh <tasklist.txt>}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"
[ -f "$TASKS" ] || { echo "no such task list: $TASKS"; exit 1; }

MEM_PCT="${MEM_PCT:-75}"; MAXJOBS="${MAXJOBS:-$(nproc)}"; FLOOR_MB="${FLOOR_MB:-8192}"
DEFAULT_EST_MB="${DEFAULT_EST_MB:-4096}"; SAFETY_PCT="${SAFETY_PCT:-125}"; TAG="${TAG:-}"
MEM_TOTAL_MB=$(awk '/^MemTotal:/{printf "%d", $2/1024}' /proc/meminfo)
BUDGET_MB=$(( MEM_TOTAL_MB * MEM_PCT / 100 ))

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PURSUE_BENCH_ROOT="$ROOT/benchmarks" PURSUE_DATA_ROOT="$ROOT/benchmarks/data"
mkdir -p logs results cache

TIME_BIN=""; [ -x /usr/bin/time ] && /usr/bin/time -f '%M' true 2>/dev/null && TIME_BIN=/usr/bin/time
RV="$(Rscript -e 'cat(as.character(getRversion()))' 2>/dev/null || echo 0)"
case "$RV" in 4.5*|4.6*|5.*) ;; *) echo ">> WARNING: R $RV -- did you forget 'conda activate pursue-bench'?";; esac

# ---- learn per-cell-class memory from every .mem this machine has ever written ---------------
# Key is axis__simulator__template__regime: the replicate does not change the table's size, the
# regime and template do. Take the MAX seen, because the budget has to survive the worst cell,
# not the average one.
declare -A EST
if ls logs/*.mem >/dev/null 2>&1; then
  while IFS=$'\t' read -r k mb; do EST["$k"]=$mb; done < <(
    for m in logs/*.mem; do
      id=$(basename "$m" .mem); kb=$(awk 'END{print $1+0}' "$m" 2>/dev/null)
      [ "${kb:-0}" -gt 0 ] && echo "$id" | awk -F'__' -v kb="$kb" 'NF>=4{printf "%s__%s__%s__%s\t%d\n",$1,$2,$3,$4,kb/1024}'
    done | awk -F'\t' '{if ($2+0 > m[$1]) m[$1]=$2+0} END{for (k in m) printf "%s\t%d\n", k, m[k]}')
fi
echo ">> $TASKS"
echo ">> RAM ${MEM_TOTAL_MB}MB total, budget ${BUDGET_MB}MB (${MEM_PCT}%), floor ${FLOOR_MB}MB, max ${MAXJOBS} jobs"
echo ">> learned memory for ${#EST[@]} cell classes from logs/*.mem (reserve ${SAFETY_PCT}% of observed, default ${DEFAULT_EST_MB}MB)"

# These run once per task line during the sort and again on every scan, so on a 13 750-line list
# a sed-per-field costs tens of thousands of processes. Pure bash, no subshells.
# NOTE: these set globals (_KEY, _EST, _AVAIL, _ID) rather than echoing. Command substitution
# forks a subshell even around a pure-bash function, and the scan below can call them a few
# hundred times per launch -- on a 13 750-cell list that was ~100 s of pure fork overhead.
key_of()   { _KEY="${1%__r[0-9][0-9][0-9]*}"; }
est_of()   { key_of "$1"; local v="${EST[$_KEY]:-}"
             if [ -n "$v" ]; then _EST=$(( v * SAFETY_PCT / 100 + 1 )); else _EST=$DEFAULT_EST_MB; fi; }
avail_mb() { local k v; while read -r k v _; do [ "$k" = "MemAvailable:" ] && { _AVAIL=$(( v / 1024 )); return; }; done < /proc/meminfo; _AVAIL=0; }

_parse() {  # sets _ax _sim _tpl _reg _rep from a task line
  local -a w=($1); local k
  _ax=; _sim=; _tpl=; _reg=; _rep=0
  for ((k=0; k<${#w[@]}-1; k++)); do
    case "${w[$k]}" in
      --axis)      _ax=${w[$((k+1))]};;  --simulator) _sim=${w[$((k+1))]};;
      --template)  _tpl=${w[$((k+1))]};; --regime)    _reg=${w[$((k+1))]};;
      --replicate) _rep=${w[$((k+1))]};;
    esac
  done
}
id_of()   { _parse "$1"; printf -v _ID "%s__%s__%s__%s__r%03d%s" "$_ax" "$_sim" "$_tpl" "$_reg" "$_rep" "${TAG:+__$TAG}"; }
axis_of() { _parse "$1"; printf "%s" "$_ax"; }

run_one() {   # $1 = task line, $2 = cell id, $3 = estimate MB
  local line="$1" id="$2" est="$3" out mf rc mem start
  out="results/axis$(axis_of "$line")"; mf="logs/$id.mem"; start=$SECONDS
  if [ -n "${CELL_CMD:-}" ]; then
    if $CELL_CMD "$line" "$id" > "logs/$id.log" 2>&1; then rc=0; else rc=$?; fi
  elif [ -n "$TIME_BIN" ]; then
    if "$TIME_BIN" -f '%M' -o "$mf" Rscript benchmarks/R/engine/run_cell.R $line --out "$out" --cache cache \
         --master-seed "${PURSUE_MASTER_SEED:-1}" --timeout "${CELL_TIMEOUT:-3600}" \
         ${TAG:+--tag "$TAG"} > "logs/$id.log" 2>&1; then rc=0; else rc=$?; fi
  else
    if Rscript benchmarks/R/engine/run_cell.R $line --out "$out" --cache cache \
         --master-seed "${PURSUE_MASTER_SEED:-1}" --timeout "${CELL_TIMEOUT:-3600}" \
         ${TAG:+--tag "$TAG"} > "logs/$id.log" 2>&1; then rc=0; else rc=$?; fi
  fi
  mem=$( [ -f "$mf" ] && awk 'END{if($1+0>0) printf "%d", $1/1024; else print 0}' "$mf" || echo 0 )
  if [ "$rc" -eq 0 ]; then echo "ok    $id  $((SECONDS-start))s  est ${est}MB used ${mem}MB"
  else case "$rc" in 137) w=" KILLED(oom?)";; 139) w=" SEGFAULT";; *) w="";; esac
       echo "FAIL  $id  $((SECONDS-start))s  est ${est}MB used ${mem}MB  rc=$rc$w  (see logs/$id.log)"; fi
}

mapfile -t RAW < <(grep -ve '^\s*$' "$TASKS")
# First-fit DECREASING: place the heaviest cells while the machine is empty and let the light ones
# fill the gaps. Ordering is free -- cells are independent -- and it is the difference between
# packing 107 GB well and running one 30 GB cell at a time.
if [ "${SORT:-desc}" = "desc" ]; then
  mapfile -t LINES < <(
    for l in "${RAW[@]}"; do id_of "$l"; est_of "$_ID"; printf '%s\t%s\n' "$_EST" "$l"; done \
      | sort -k1,1nr -s | cut -f2-)
else LINES=("${RAW[@]}"); fi

total=${#LINES[@]}; head=0; skip_n=0; RESERVED=0; peak_jobs=0
declare -A PID_EST PID_ID PID_KEY
declare -a TAKEN; for ((k=0;k<total;k++)); do TAKEN[$k]=0; done
echo ">> $total cells, started $(date '+%F %T')"

reap() {   # free the reservation of anything that has exited, and learn its real cost
  local pid
  for pid in "${!PID_ID[@]}"; do
    if ! kill -0 "$pid" 2>/dev/null; then
      wait "$pid" 2>/dev/null
      local id="${PID_ID[$pid]}" k="${PID_KEY[$pid]}" mb=0
      [ -f "logs/$id.mem" ] && mb=$(awk 'END{printf "%d", $1/1024}' "logs/$id.mem" 2>/dev/null)
      if [ "${mb:-0}" -gt "${EST[$k]:-0}" ]; then EST["$k"]=$mb; fi
      RESERVED=$(( RESERVED - PID_EST[$pid] )); [ "$RESERVED" -lt 0 ] && RESERVED=0
      unset 'PID_EST[$pid]' 'PID_ID[$pid]' 'PID_KEY[$pid]'
    fi
  done
}

# PENDING is maintained rather than recounted: recounting from head on every loop turn made the
# scheduler O(n^2) and cost 99 s of pure bookkeeping on a 13 750-cell list.
PENDING=$total
while [ "$PENDING" -gt 0 ] || [ "${#PID_ID[@]}" -gt 0 ]; do
  reap
  launched_any=0
  while :; do
    running=${#PID_ID[@]}
    [ "$running" -ge "$MAXJOBS" ] && break
    avail_mb; picked=-1
    # Scan a window rather than stopping at the first cell that does not fit: a single heavy cell
    # at the head must not stall the light ones behind it (the 2026-09-16 bug -- peak concurrency
    # of 1 on a mixed list). The window keeps the scan cheap on a 13 000-line task list.
    for ((j=head; j<total && j<head+${WINDOW:-400}; j++)); do
      [ "${TAKEN[$j]}" -eq 1 ] && continue
      line="${LINES[$j]}"; id_of "$line"; id=$_ID; ax=$_ax
      if [ -f "results/axis$ax/$id.manifest.json" ] || [ -f "results/axis$ax/$id.skipped.json" ]; then
        echo "skip  $id"; TAKEN[$j]=1; PENDING=$((PENDING-1)); skip_n=$((skip_n+1)); continue
      fi
      est_of "$id"; est=$_EST
      # Reserve at most the whole budget: a cell estimated above it still has to run, alone.
      [ "$est" -gt "$BUDGET_MB" ] && est=$BUDGET_MB
      if [ "$running" -eq 0 ] || { [ $(( RESERVED + est )) -le "$BUDGET_MB" ] && [ $(( _AVAIL - est )) -ge "$FLOOR_MB" ]; }; then
        picked=$j; break
      fi
    done
    [ "$picked" -lt 0 ] && break
    line="${LINES[$picked]}"; id_of "$line"; id=$_ID; est_of "$id"; est=$_EST
    [ "$est" -gt "$BUDGET_MB" ] && est=$BUDGET_MB
    TAKEN[$picked]=1; PENDING=$((PENDING-1))
    while [ "$head" -lt "$total" ] && [ "${TAKEN[$head]}" -eq 1 ]; do head=$((head+1)); done
    if [ -n "${DRY:-}" ]; then
      echo "would run  $id  est ${est}MB  (reserved $(( RESERVED + est ))/${BUDGET_MB}MB, avail ${_AVAIL}MB)"
      continue
    fi
    run_one "$line" "$id" "$est" &
    pid=$!; key_of "$id"; PID_EST[$pid]=$est; PID_ID[$pid]=$id; PID_KEY[$pid]=$_KEY
    RESERVED=$(( RESERVED + est )); launched_any=1
    [ "${#PID_ID[@]}" -gt "$peak_jobs" ] && peak_jobs=${#PID_ID[@]}
  done
  [ -n "${DRY:-}" ] && [ "$PENDING" -eq 0 ] && break
  if [ "$launched_any" -eq 0 ]; then sleep "${POLL_S:-3}"; fi
done
wait
echo ">> finished $(date '+%F %T'); peak concurrency $peak_jobs, $skip_n skipped"
d=$(ls results/axis*/ 2>/dev/null | grep -c 'manifest.json' || true)
s=$(ls results/axis*/ 2>/dev/null | grep -c 'skipped.json' || true)
echo ">> $d cells have results; $s not producible by that simulator."
echo ">> Re-run this command to retry anything that failed."
