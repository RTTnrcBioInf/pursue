#!/usr/bin/env bash
# Submit the full benchmark in dependency order. Run from the repository root after
# setup_env.sh, download_templates.sh and smoke_test.R have succeeded.
set -euo pipefail
Rscript hpc/make_tasklist.R --pool "${POOL:-evaluation}" --simulators "${SIMS:-house,msq,mid,sd2,sps}" --methods "${METHODS:-all}"
nA=$(wc -l < hpc/tasks/axisA.txt); nB=$(wc -l < hpc/tasks/axisB.txt); nC=$(wc -l < hpc/tasks/axisC.txt); nR=$(wc -l < hpc/tasks/realism.txt)
jR=$(sbatch --parsable --array=1-$nR hpc/slurm/realism.sbatch)
jA=$(sbatch --parsable --dependency=afterok:$jR --array=1-$nA%${MAXC:-200} hpc/slurm/axisA.sbatch)
jB=$(sbatch --parsable --array=1-$nB%${MAXC:-200} hpc/slurm/axisB.sbatch)
jC=$(sbatch --parsable --array=1-$nC%${MAXC:-200} hpc/slurm/axisC.sbatch)
jD=$(sbatch --parsable hpc/slurm/axisDE.sbatch)
echo "submitted: realism=$jR axisA=$jA axisB=$jB axisC=$jC axisDE=$jD"
echo "when done:  Rscript benchmarks/R/analysis/aggregate.R results"
