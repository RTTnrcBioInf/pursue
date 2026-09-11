#!/usr/bin/env bash
#SBATCH --job-name=pursue-B
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=logs/axisB_%A_%a.out
#SBATCH --error=logs/axisB_%A_%a.err
# usage: sbatch --array=1-$(wc -l < hpc/tasks/axisB.txt)%200 hpc/slurm/axisB.sbatch
set -euo pipefail
source hpc/slurm/env.sh
ARGS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" hpc/tasks/axisB.txt)
echo "[$(date)] task $SLURM_ARRAY_TASK_ID: $ARGS"
Rscript benchmarks/R/engine/run_cell.R $ARGS --out results/axisB --cache cache --n-cores "$SLURM_CPUS_PER_TASK" --master-seed "${PURSUE_MASTER_SEED:-1}"
