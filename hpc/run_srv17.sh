#!/usr/bin/env bash
# run_srv17.sh -- parallel development track: is p13's sd2 / sps gap visible on the TUNING templates?
# eperm_cs1 and the reference comparators (ZicoSeq, LDM, ADAPT, LinDA) on the new sd2 / sps dev settings
# (R00, R06, R08, R13 x hmp_tongue, twinsuk_stool x 5 reps), then msq R00/R11 for the same references at 10
# reps. 16 cores by default, so it can share the machine with hpc/run_confirm.sh (40). Pushes its results.
#   git pull && mkdir -p logs && (nohup bash hpc/run_srv17.sh > logs/run_srv17.log 2>&1 &)
# Resumable: cells already checkpointed under benchmarks/dev/results/srv17/cells are not re-run.
#
# 2026-10-06: every worker single-threaded (the first launch let each forked worker start a multi-threaded
# BLAS/OpenMP pool: load 620 on 64 cores), and the sparseDOSSA2 fit of each tuning template made ONCE, before
# the suite, in its own process (1-2.5 h each on the evaluation templates); otherwise every worker that draws
# an sd2 cell misses the cache at the same moment and refits the same template in parallel.
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results cache
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
L=srv17; R=benchmarks/dev/results/$L
C=eperm_cs1,ref_zicoseq,ref_ldm,ref_adapt,ref_linda

echo ">> warm the sd2 fits of the tuning templates ($(date '+%F %T'))"
for t in hmp_tongue twinsuk_stool; do
  if [ -f "cache/sd2_$t.rds" ]; then echo "   cache/sd2_$t.rds exists"; continue; fi
  # two fits, 8 threads each: 16 cores, the same share as the suite
  ( OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8 Rscript -e '
    root <- normalizePath("benchmarks"); Sys.setenv(PURSUE_BENCH_ROOT = root)
    for (f in c("R/engine/templates.R", "R/engine/regimes.R", "R/engine/metrics.R", "R/simulators/sim_house.R",
                "R/simulators/implant.R", "R/simulators/dispatch.R", "R/methods/elementary.R")) source(file.path(root, f))
    t <- commandArgs(TRUE)[1]; tpl <- readRDS(file.path(root, "devdata", paste0(t, ".rds")))
    reg <- read.delim(file.path(root, "regimes.tsv"), stringsAsFactors = FALSE, comment.char = "#")
    t0 <- Sys.time(); x <- simulate_cell_data("sd2", tpl, reg[reg$regime_id == "R00", ], 1L, normalizePath("cache"))
    cat(t, "sd2 fit cached in", round(as.numeric(Sys.time() - t0, units = "mins")), "min; cell has", nrow(x$counts), "features\n")' "$t" \
    > "logs/warm_sd2_$t.log" 2>&1; echo "   $t: $(tail -1 logs/warm_sd2_$t.log)" ) &
done
wait
ls -la cache/sd2_hmp_tongue.rds cache/sd2_twinsuk_stool.rds

echo ">> dev suite $L ($(date '+%F %T'))"
Rscript benchmarks/dev/devsuite.R --suite full --sims sd2,sps --reps 5 --cores "${CORES:-16}" --label "$L" --candidates "$C"
Rscript benchmarks/dev/devsuite.R --suite full --sims msq --settings msq:R00,msq:R11 --reps 10 --cores "${CORES:-16}" --label "${L}m" --candidates "$C"
{ Rscript benchmarks/dev/compare.R --labels "$L" --brief; Rscript benchmarks/dev/compare.R --labels "${L}m" --brief; } > "$R/compare.txt" 2>&1
cp logs/run_srv17.log "$R/run.log" 2>/dev/null || true
git add "$R" "benchmarks/dev/results/${L}m" && git commit -m "dev suite $L: sd2 / sps / msq gap on the tuning templates (eperm_cs1 vs references)" \
  && git pull --rebase --autostash && git push
echo ">> done $(date '+%F %T')"
