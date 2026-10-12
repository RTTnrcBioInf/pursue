#!/usr/bin/env bash
# run_srv24.sh -- comparator targets on the tuning templates where comparators still lead (p14): msq (ZicoSeq / LDM)
# and the implants (logistic presence; ZicoSeq / LDM on B_cont). eperm_c5q beside the five references, so every gap is
# paired, cell by cell. Then bank3: implant B_cont / B_conf07 / B_prev cells with per-feature p of c5q, LDM,
# ZicoSeq and logistic presence, for local diagnosis (formulas as text; small).
#   srv24:  msq R00 / R11 / R17 / R22, implant B_ref / B_prev / B_abund / B_allup / B_conf04 / B_conf07 / B_cont; 10 reps.
# Single-threaded BLAS exported; cores from CORES. Resumable; pushes its results.
#   git pull && mkdir -p logs && (CORES=24 nohup nice -n 10 bash hpc/run_srv24.sh > logs/run_srv24.log 2>&1 &)
set -o pipefail
# Shared server: if memory runs out, the kernel's OOM killer should take these jobs first, not other users' (2026-10-12).
# A non-root process may raise its own oom_score_adj; every child (Rscript, its forked workers) inherits it.
echo 1000 > /proc/self/oom_score_adj 2>/dev/null || true
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
CORES="${CORES:-48}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
C=eperm_c5q,ref_zicoseq,ref_ldm,ref_adapt,ref_linda,ref_logistic_presence
L=srv24; R=benchmarks/dev/results/$L
echo ">> dev suite $L ($(date '+%F %T'))"
Rscript benchmarks/dev/devsuite.R --suite full --sims implant,msq --settings msq:R00,msq:R11,msq:R17,implant:B_ref,implant:B_prev,implant:B_abund,implant:B_allup \
  --extra msq:R22,implant:B_conf04,implant:B_conf07,implant:B_cont --reps 10 --cores "$CORES" --label "$L" --candidates "$C" || exit 1
Rscript benchmarks/dev/compare.R --labels "$L" > "$R/compare.txt" 2>&1
cp logs/run_srv24.log "$R/run.log" 2>/dev/null || true
git add "$R" && git commit -m "dev suite $L: comparator targets on msq and implants (c5q beside ZicoSeq, LDM, ADAPT, LinDA, logistic presence)" && git pull --rebase --autostash && git push
O=benchmarks/dev/bank/bank3
echo ">> bank3 ($(date '+%F %T'))"
Rscript benchmarks/dev/tools/bank_cells.R --settings implant:B_cont,implant:B_conf07,implant:B_prev --reps 5 \
  --candidates eperm_c5q,ref_ldm,ref_zicoseq,ref_logistic_presence --out "$O" --cores "$CORES" && du -sh "$O" \
  && git add "$O" && git commit -m "bank3: B_cont / B_conf07 / B_prev cells with c5q, LDM, ZicoSeq, logistic presence p-values" && git pull --rebase --autostash && git push
echo ">> done $(date '+%F %T')"
