#!/usr/bin/env bash
#SBATCH --job-name=pursue-DE
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/axisDE_%j.out
#SBATCH --error=logs/axisDE_%j.err
set -euo pipefail
source hpc/slurm/env.sh
R=benchmarks/R/engine/run_axisDE.R; O=results/axisDE; mkdir -p $O
Rscript $R --axis D --template mbd_gingival_v35 --group body_subsite --expected expected/gingival_aerobes.tsv --out $O --n-cores 4
Rscript $R --axis D --template mbd_ravel_bv --group study_condition --expected expected/bv_taxa.tsv --out $O --n-cores 4
Rscript $R --axis D --template mbd_stammler_spikein --spikein "$(cat benchmarks/expected/stammler_spikein_ids.txt)" --n-groupings 50 --out $O
Rscript $R --axis E --template crc_genus --group diagnosis --levels control,CRC --splits 5 --out $O --n-cores 4
Rscript $R --axis E --template risk_ileum --group diagnosis --levels no,CD --splits 5 --out $O --n-cores 4
