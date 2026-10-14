#!/usr/bin/env bash
# run_srv29.sh -- covariate adjustment, first look (tuning templates). On the evaluation templates c12 trails ADAPT by
# ~15% in covariate-confounded regimes (R20 / R21: house, msq) with FDR 0.01-0.03, while on the implants (B_conf04 / 07)
# it is level with the field. Locally, house R20 / R21 at 3 reps: null-taxon FPR 0.056-0.061, so not obviously
# conservative per taxon. This run gives the comparator targets on the tuning templates, the covariate nulls, and banked
# cells with per-feature p for local diagnosis.
#   srv29:  house R20 / R21, msq R20 / R21, B_conf04 / B_conf07; c12 beside ADAPT, LinDA, ZicoSeq, LDM, logistic presence; 10 reps.
#   srv29n: the same six settings as nulls (20 reps), same methods.
#   bank4:  house R20 / R21 and msq R20 cells (5 reps per template) with per-feature p of c12, ADAPT, LinDA, LDM, ZicoSeq.
# hpc/guard.sh: shared memory cap (MEMTOTAL), OOM-killer first choice, single-threaded BLAS. Resumable; pushes after each stage.
#   git pull && mkdir -p logs && (MEMTOTAL=130G CORES=24 nohup nice -n 10 bash hpc/run_srv29.sh > logs/run_srv29.log 2>&1 &)
source "$(dirname "${BASH_SOURCE[0]}")/guard.sh" "$@"
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
CORES="${CORES:-24}"
push() { git add "$@" && git commit -q -m "$MSG" && git pull -q --rebase --autostash && git push -q; }
C=eperm_c12,ref_adapt,ref_linda,ref_zicoseq,ref_ldm,ref_logistic_presence
L=srv29; R=benchmarks/dev/results/$L
echo ">> dev suite $L ($(date '+%F %T'))"
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq --settings house:R21 \
  --extra house:R20,msq:R20,msq:R21,implant:B_conf04,implant:B_conf07 --reps 10 --cores "$CORES" --label "$L" --candidates "$C" || exit 1
Rscript benchmarks/dev/compare.R --labels "$L" > "$R/compare.txt" 2>&1
MSG="dev suite $L: covariate settings, c12 beside five references"; push "$R"
LN=srv29n; RN=benchmarks/dev/results/$LN
echo ">> null run $LN ($(date '+%F %T'))"
Rscript benchmarks/dev/devsuite.R --suite full --sims house,implant,msq --settings none \
  --extra house:R20.null,house:R21.null,msq:R20.null,msq:R21.null,implant:B_conf04.null,implant:B_conf07.null \
  --reps 20 --cores "$CORES" --label "$LN" --candidates "$C" \
  && Rscript benchmarks/dev/tools/null_summary.R "$RN/cells.csv" > "$RN/null_summary.txt" 2>&1
MSG="dev suite $LN: covariate nulls, c12 beside five references"; push "$RN"
O=benchmarks/dev/bank/bank4
echo ">> bank4 ($(date '+%F %T'))"
Rscript benchmarks/dev/tools/bank_cells.R --settings house:R20,house:R21,msq:R20 --reps 5 \
  --candidates eperm_c12,ref_adapt,ref_linda,ref_ldm,ref_zicoseq --out "$O" --cores "$CORES" && du -sh "$O"
cp logs/run_srv29.log "$RN/run.log" 2>/dev/null || true
MSG="bank4: covariate cells with per-feature p of c12, ADAPT, LinDA, LDM, ZicoSeq"; push "$O" "$RN"
echo ">> done $(date '+%F %T')"
