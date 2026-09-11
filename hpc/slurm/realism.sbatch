#!/usr/bin/env bash
#SBATCH --job-name=pursue-realism
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=logs/realism_%A_%a.out
#SBATCH --error=logs/realism_%A_%a.err
# usage: sbatch --array=1-$(wc -l < hpc/tasks/realism.txt) hpc/slurm/realism.sbatch
# Also fits and caches the simulator templates (sparseDOSSA2 fits take hours), so run this
# BEFORE the axis A array.
set -euo pipefail
source hpc/slurm/env.sh
read SIM TPL < <(sed -n "${SLURM_ARRAY_TASK_ID}p" hpc/tasks/realism.txt)
mkdir -p results/realism cache
Rscript -e "
Sys.setenv(PURSUE_BENCH_ROOT='benchmarks'); setwd('benchmarks')
for (f in c('R/engine/templates.R','R/engine/regimes.R','R/engine/realism.R','R/simulators/sim_house.R','R/simulators/dispatch.R')) source(f)
tpl <- load_template('$TPL'); g <- realism_gate(tpl, '$SIM', n_sim = 200, seed = 1, cache_dir = '../cache')
write.csv(g, '../results/realism/${SIM}__${TPL}.csv', row.names = FALSE); print(g)"
