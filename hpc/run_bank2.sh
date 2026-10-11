#!/usr/bin/env bash
# run_bank2.sh -- bank mid cells under designed depth confounding for local work on the efull path (the one place the
# dev suite's calibration is still weak: mid R19.null P(any) 0.16, mid R17 FDR 0.124 for every c3-based build):
# mid R17.null, mid R19.null, mid R17; both tuning templates; 10 reps; per-feature p of eperm_cs1 (per option) and
# eperm_c5. Small (formulas stored as text); pushes.
#   git pull && mkdir -p logs && (nohup bash hpc/run_bank2.sh > logs/run_bank2.log 2>&1 &)
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
O=benchmarks/dev/bank/bank2
Rscript benchmarks/dev/tools/bank_cells.R --settings mid:R17.null,mid:R19.null,mid:R17 --reps 10 \
  --candidates eperm_cs1,eperm_c5 --out "$O" --cores "${CORES:-24}" || exit 1
du -sh "$O"
git add "$O" && git commit -m "bank2: mid R17 / R19 depth-confounded cells (nulls + R17) with eperm_cs1 / c5 p-values" && git pull --rebase --autostash && git push
