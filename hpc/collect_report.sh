#!/usr/bin/env bash
# collect_report.sh -- gather everything worth looking at after a run into report/, small enough
# to travel through git. Run from the repo root; safe to re-run (report/ is rebuilt each time).
#
#   bash hpc/collect_report.sh            # then commit and push report/
#
# Deliberately does NOT copy the per-cell feature tables: those are the bulk of results/ and are
# only needed once a specific question points at specific cells.
set -uo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
R=report; rm -rf "$R"; mkdir -p "$R/summary"
say() { printf '  %-34s %s\n' "$1" "$2"; }
echo ">> collecting into $R/"

# --- 1. run logs (the ok/FAIL/skip lines, with peak RSS where instrumented) -----------------
for f in logs/axisA.run logs/axisA.resume logs/axisA_final.run logs/axisA_big.run \
         logs/axisB.run logs/axisC.run; do
  [ -f "$f" ] && { gzip -c "$f" > "$R/$(basename "$f").gz"; say "$(basename "$f")" "$(wc -l < "$f") lines"; }
done

# --- 2. aggregate output --------------------------------------------------------------------
for f in results/summary/*.csv results/summary/*.txt results/summary/*.log; do
  [ -f "$f" ] || continue
  b=$(basename "$f"); sz=$(stat -c%s "$f" 2>/dev/null || echo 0)
  # metrics_long is the one big file; ship it compressed, and only if it stays sane
  if [ "$b" = "metrics_long.csv" ]; then
    gzip -c "$f" > "$R/summary/$b.gz"
    gz=$(stat -c%s "$R/summary/$b.gz" 2>/dev/null || echo 0)
    if [ "$gz" -gt 18000000 ]; then rm -f "$R/summary/$b.gz"; say "$b" "SKIPPED (${gz}B gzipped, too big)"
    else say "$b" "$((sz/1048576))MB -> $((gz/1048576))MB gzipped"; fi
  else
    cp "$f" "$R/summary/"; say "$b" "$((sz/1024))KB"
  fi
done

# --- 3. peak memory per cell ----------------------------------------------------------------
# run_local.sh writes logs/<cell>.mem (GNU time %M, in KB). Cell id encodes axis/sim/template/
# regime, so this is the memory-by-design table we never had when cells were being killed.
if ls logs/*.mem >/dev/null 2>&1; then
  { echo "axis,simulator,template,regime_id,replicate,peak_rss_mb"
    for m in logs/*.mem; do
      id=$(basename "$m" .mem); kb=$(awk 'END{print $1+0}' "$m")
      echo "$id" | awk -F'__' -v kb="$kb" 'NF>=5{sub("^r","",$5); printf "%s,%s,%s,%s,%s,%.0f\n",$1,$2,$3,$4,$5,kb/1024}'
    done
  } > "$R/memory_by_cell.csv"
  say "memory_by_cell.csv" "$(( $(wc -l < "$R/memory_by_cell.csv") - 1 )) cells"
  { echo "== peak RSS MB by template =="
    awk -F, 'NR>1{v[$3]=v[$3]" "$6} END{for(t in v){n=split(v[t],a," "); asort_n=0
      for(i=1;i<=n;i++){s+=a[i]; if(a[i]+0>mx)mx=a[i]+0}
      printf "%-22s max %6.0f  mean %6.0f  n %5d\n", t, mx, s/n, n; s=0; mx=0}}' "$R/memory_by_cell.csv"
    echo; echo "== 15 heaviest cells =="
    tail -n +2 "$R/memory_by_cell.csv" | sort -t, -k6 -rn | head -15
  } > "$R/memory_summary.txt"
else
  echo "no logs/*.mem -- run_local.sh was not the instrumented version" > "$R/memory_summary.txt"
fi
say "memory_summary.txt" "written"

# --- 4. what each method actually returned, across every cell log ---------------------------
# Each cell log has one "  <method>  <seconds>s  <status>" line per method. Cheaper and more
# complete than re-reading the feature tables.
if ls logs/A__*.log >/dev/null 2>&1; then
  grep -h -E '^  [a-z0-9_]+ +[0-9.]+s ' logs/*.log 2>/dev/null \
    | awk '{m=$1; t=$2; $1=""; $2=""; sub(/^ +/,""); s=$0; key=m"\t"s; c[key]++; sec[m]+=t; n[m]++}
           END{for(k in c) printf "%7d\t%s\n", c[k], k}' \
    | sort -rn > "$R/method_status_counts.tsv"
  say "method_status_counts.tsv" "$(wc -l < "$R/method_status_counts.tsv") distinct method/status pairs"
  grep -h -E '^  [a-z0-9_]+ +[0-9.]+s ' logs/*.log 2>/dev/null \
    | awk '{sec[$1]+=$2; n[$1]++} END{printf "%-20s %8s %10s %12s\n","method","cells","mean_s","total_CPUh"
            for(m in sec) printf "%-20s %8d %10.1f %12.1f\n", m, n[m], sec[m]/n[m], sec[m]/3600}' \
    | (read -r h; echo "$h"; sort -k4 -rn) > "$R/method_runtime.txt"
  say "method_runtime.txt" "per-method cost"
fi

# --- 5. inventory ---------------------------------------------------------------------------
{ echo "== results inventory =="
  for d in results/axis*; do
    [ -d "$d" ] || continue
    printf "%-16s %6d manifests  %6d not-producible  %6d metrics\n" "$(basename "$d")" \
      "$(ls "$d" | grep -c 'manifest\.json$')" "$(ls "$d" | grep -c 'skipped\.json$')" \
      "$(ls "$d" | grep -c 'metrics\.csv$')"
  done
  echo; echo "== disk =="; du -sh results logs cache 2>/dev/null
  echo; echo "== git =="; git rev-parse HEAD; git status --short | head -20
  echo; echo "== R =="; Rscript -e 'cat(R.version.string, "\n")' 2>/dev/null
} > "$R/inventory.txt"
say "inventory.txt" "written"

# --- 6. the fastEmu reference-set question ---------------------------------------------------
[ -f hpc/fastemu_refset.R ] && { Rscript hpc/fastemu_refset.R > "$R/fastemu_refset.txt" 2>&1
                                 say "fastemu_refset.txt" "written"; }

echo
echo ">> total: $(du -sh "$R" | cut -f1)"
echo ">> push it:"
echo "     git add -f report/ && git commit -m 'run report' && git push"
