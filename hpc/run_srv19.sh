#!/usr/bin/env bash
# run_srv19.sh -- Phase 2 after srv18 (benchmarks/dev/candidates/70_epg.R):
#   srv19:  eperm_c2 (eperm_cs1 + the continuous-exposure permutation path + detection as a selection option) and
#           eperm_c2q (the same with Storey's adaptive BH as its own FDR procedure) and eperm_c3 (eperm_c2 without the
#           detection kernel under designed depth confounding) on the whole dev suite + implant B_conf04 / B_conf07 /
#           B_cont + msq R22, 10 reps (c2q and c3 reuse eperm_c2's fit, so they cost almost nothing extra);
#   srv19n: nulls at 40 reps -- the continuous path (house R22, implant B_cont, msq R22), the depth-confounded nulls
#           (house / msq / mid R19, mid R17) for eperm_c3, and the reference nulls (house R06, implant B_null, msq R06,
#           mid R06), where the adaptive pi0 must stay at 1;
#   bank1:  cells banked with p-values (dev/tools/bank_cells.R) for local diagnosis of the gaps that need the server's
#           simulators: sd2 R00/R08/R13 (hmp_tongue), msq R00/R11 and mid R19 (both templates), 5 reps;
#           eperm_cs1 (per-option p), LDM, ZicoSeq.
# Development only (the benchmark does not source these files). 16 cores by default; resumable; pushes its results.
#   git pull && mkdir -p logs && (nohup bash hpc/run_srv19.sh > logs/run_srv19.log 2>&1 &)
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results benchmarks/dev/bank
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
CORES="${CORES:-16}"; C=eperm_c2,eperm_c2q,eperm_c3
B=benchmarks/dev/bank/bank1
echo ">> bank $B ($(date '+%F %T'))"
Rscript benchmarks/dev/tools/bank_cells.R --settings sd2:R00,sd2:R08,sd2:R13 --templates hmp_tongue --reps 5 \
  --candidates eperm_cs1,ref_ldm,ref_zicoseq --out "$B" --cores "$CORES"
Rscript benchmarks/dev/tools/bank_cells.R --settings msq:R00,msq:R11,mid:R19 --templates hmp_tongue,twinsuk_stool --reps 5 \
  --candidates eperm_cs1,ref_ldm,ref_zicoseq --out "$B" --cores "$CORES"
git add "$B" && git commit -m "dev bank1: sd2 / msq / mid R19 cells with eperm_cs1, LDM, ZicoSeq p-values" && git pull --rebase --autostash && git push
L=srv19; R=benchmarks/dev/results/$L
echo ">> dev suite $L ($(date '+%F %T'))"
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq,mid --extra implant:B_conf04,implant:B_conf07,implant:B_cont,msq:R22 \
  --reps 10 --cores "$CORES" --label "$L" --candidates "$C" || exit 1
Rscript benchmarks/dev/compare.R --labels "srv18,$L" > "$R/compare.txt" 2>&1
cp logs/run_srv19.log "$R/run.log" 2>/dev/null || true
git add "$R" && git commit -m "dev suite $L: eperm_c2, eperm_c2q (adaptive BH), eperm_c3 (no detection under depth confounding)" && git pull --rebase --autostash && git push
LN=srv19n; RN=benchmarks/dev/results/$LN
echo ">> null run $LN ($(date '+%F %T'))"
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq,mid --settings house:R06,implant:B_null,msq:R06,mid:R06 \
  --extra house:R22.null,implant:B_cont.null,msq:R22.null,house:R19.null,msq:R19.null,mid:R19.null,mid:R17.null --reps 40 --cores "$CORES" --label "$LN" --candidates "$C" \
  && Rscript benchmarks/dev/tools/null_summary.R "$RN/cells.csv" > "$RN/null_summary.txt" 2>&1
cp logs/run_srv19.log "$RN/run.log" 2>/dev/null || true
git add "$RN" && git commit -m "dev suite $LN: nulls for eperm_c2 / eperm_c2q / eperm_c3" && git pull --rebase --autostash && git push
echo ">> done $(date '+%F %T')"
