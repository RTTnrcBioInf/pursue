#!/usr/bin/env bash
# run_p14.sh -- the outer loop for it48's lead: full benchmark (axes A-C, evaluation pool, seed 1, all five simulators)
# for pursue03_epermc5 (eperm_c5: pos kernel + selection floor, BH) and pursue03_epermc5q (+ Storey's adaptive BH when
# plain BH rejects), as a tagged solo re-run beside the stored comparator results; aggregated with eperm_cs1 (p13) in
# focus. Dev suite (srv20/21, 44 settings, excl. R25): eperm_c3 9.57, c5 9.76, c5q 10.11 TP per cell; null P(any) on
# srv20n's 440 cells c3 0.052, c5 / c5q 0.045. Resumable (run_adaptive.sh skips finished cells); pushes results/summary.
#   git pull && mkdir -p logs && (nohup bash hpc/run_p14.sh > logs/run_p14.log 2>&1 &)
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs
# No OMP_NUM_THREADS export here: GNU nproc honours it, so run_adaptive.sh computed MAXJOBS = 1 (first p14 launch,
# 2026-10-11). run_adaptive.sh sets the thread variables itself, after it has counted the cores.
export PURSUE_FOCUS="pursue03_epermcs1,pursue03_epermc5,pursue03_epermc5q"
M=pursue03_epermc5,pursue03_epermc5q TAG=p14 bash hpc/run_p03.sh
