#!/usr/bin/env bash
#SBATCH --job-name=pursue-A
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/axisA_%A_%a.out
#SBATCH --error=logs/axisA_%A_%a.err
# usage: sbatch --array=1-$(wc -l < hpc/tasks/axisA.txt)%200 hpc/slurm/axisA.sbatch
set -euo pipefail
source hpc/slurm/env.sh
ARGS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" hpc/tasks/axisA.txt)
echo "[$(date)] task $SLURM_ARRAY_TASK_ID: $ARGS"
Rscript benchmarks/R/engine/run_cell.R $ARGS --out results/axisA --cache cache --n-cores "$SLURM_CPUS_PER_TASK" --master-seed "${PURSUE_MASTER_SEED:-1}"
