#!/usr/bin/env bash
# run_srv28.sh -- it56: eperm_c12 = c10 with the pos kernel centred on its own median taxon (dev/candidates/81_own.R).
# srv27: c10 failed H_all (spread-only, every taxon) on twinsuk_stool because the selection moved every taxon to pos,
# whose bulk is shifted under the log kernel's centre. Locally c12 makes 0 calls on the two worst stool cells (c10: 163,
# 68) and keeps pos's sd2 gain (10 banked cells: TP / FP 131 / 5 vs c10 129 / 3). c10 is bit-identical (option off).
#   srv28s: spread suite, c12 alone (srv27 has c10 and the comparators on the same seeds); 20 reps per template.
#   srv28:  whole dev suite + sd2 / sps + B_conf04 / B_conf07 / B_cont + msq R22, c10 beside c12; 10 reps.
#   srv28n: two-group nulls at 40 reps (house / msq / mid / sd2 / sps R06, B_null) and bloom x4 at 20 (srv28nb).
# hpc/guard.sh: shared memory cap (MEMTOTAL), OOM-killer first choice, single-threaded BLAS. Resumable; pushes after
# each stage.
#   git pull && mkdir -p logs && (MEMTOTAL=130G CORES=24 nohup nice -n 10 bash hpc/run_srv28.sh > logs/run_srv28.log 2>&1 &)
source "$(dirname "${BASH_SOURCE[0]}")/guard.sh" "$@"
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
CORES="${CORES:-24}"
push() { git add "$@" && git commit -q -m "$MSG" && git pull -q --rebase --autostash && git push -q; }
LS=srv28s; RS=benchmarks/dev/results/$LS
echo ">> spread suite $LS ($(date '+%F %T'))"
Rscript benchmarks/dev/spread_suite.R --candidates eperm_c12 --reps 20 --cores "$CORES" --label "$LS" || exit 1
MSG="dev suite $LS: spread-only implants for eperm_c12"; push "$RS"
L=srv28; R=benchmarks/dev/results/$L
echo ">> dev suite $L ($(date '+%F %T'))"
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq,mid,sd2,sps --extra implant:B_conf04,implant:B_conf07,implant:B_cont,msq:R22 \
  --reps 10 --cores "$CORES" --label "$L" --candidates eperm_c10,eperm_c12 || exit 1
Rscript benchmarks/dev/compare.R --labels "$L" > "$R/compare.txt" 2>&1
MSG="dev suite $L: eperm_c12 beside eperm_c10"; push "$R"
LN=srv28n; RN=benchmarks/dev/results/$LN
echo ">> null run $LN ($(date '+%F %T'))"
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq,mid,sd2,sps --settings sd2:R06,sps:R06,house:R06,implant:B_null,msq:R06,mid:R06 \
  --reps 40 --cores "$CORES" --label "$LN" --candidates eperm_c10,eperm_c12 \
  && Rscript benchmarks/dev/tools/null_summary.R "$RN/cells.csv" > "$RN/null_summary.txt" 2>&1
Rscript benchmarks/dev/devsuite.R --suite full --sims house --settings bloom:x4 --reps 20 --cores "$CORES" --label "${LN}b" --candidates eperm_c10,eperm_c12
cp logs/run_srv28.log "$RN/run.log" 2>/dev/null || true
MSG="dev suite $LN: two-group nulls and bloom for eperm_c10 / c12"; push "$RN" "benchmarks/dev/results/${LN}b"
echo ">> done $(date '+%F %T')"
