#!/usr/bin/env bash
# run_srv36.sh -- it59 on the dev suite (benchmarks/dev/candidates/85_c15.R; tuning templates, dev seed 9001):
#   eperm_c15c  c14 with the compositional centre as a |kap|-weighted median, re-solved on reference taxa (|centred t| <= 2)
#               -- the one-directional gap (mid R10, B_allup). Same rule in the observed statistic and every permutation.
#   (it60, the plain-proportion kernel, failed the bloom and R25 locally and is not run.)
# 1. srv36: srv35's cells (whole dev suite + the larger-n msq settings, 10 reps) -- compare.R against c14, LDM, MaAsLin 3, LOCOM2;
# 2. srv36n: srv35n's global nulls (40 reps);
# 3. srv36s: the spread suite (H_unb / H_bal / H_all / H_cont, 20 reps; srv28s has c12 = c14 on the same cells) as a
#    regression check: the centre moves every taxon's statistic.
# Needs srv35 pushed first (its cells.csv are the reference). ~40 CPU-hours: ~3.5 h at 12 cores. hpc/guard.sh as usual.
#   git pull && mkdir -p logs && (MEMTOTAL=90G CORES=12 nohup nice -n 10 bash hpc/run_srv36.sh > logs/run_srv36.log 2>&1 &)
source "$(dirname "${BASH_SOURCE[0]}")/guard.sh" "$@"
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
CORES="${CORES:-12}"; R=benchmarks/dev/results; C="${C:-eperm_c15c}"
SUITE="--suite full --sims house,implant,msq,mid,sd2,sps"
EXTRA="msq:R02,msq:R03,msq:R09,msq:R10,msq:R13,msq:R22,mid:R02,mid:R10,mid:R13,house:R10,implant:B_cont,implant:B_conf04,implant:B_conf07"
NULLS="house:R06,msq:R06,mid:R06,sd2:R06,sps:R06,implant:B_null"
push() { git add "$@" && git commit -q -m "$MSG" && git pull -q --rebase --autostash && git push -q; }
[ -f "$R/srv35/cells.csv" ] || { echo ">> srv35 results not found -- run hpc/run_srv35.sh first"; exit 1; }
echo ">> srv36 ($(date '+%F %T')): $C, $CORES cores"
L=srv36; Rscript benchmarks/dev/devsuite.R $SUITE --extra "$EXTRA" --reps 10 --cores "$CORES" --label "$L" --candidates "$C" || exit 1
Rscript benchmarks/dev/compare.R --labels "srv35,srv35m,$L" > "$R/$L/compare.txt" 2>&1
rm -rf "$R/$L/cells"; MSG="srv36: it59 (c15c, weighted compositional centre) on the dev suite" push "$R/$L"
L=srv36n; Rscript benchmarks/dev/devsuite.R $SUITE --settings "$NULLS" --reps 40 --cores "$CORES" --label "$L" --candidates "$C" || exit 1
Rscript -e 'a <- read.csv("benchmarks/dev/results/srv35n/cells.csv"); b <- read.csv("benchmarks/dev/results/srv36n/cells.csv"); write.csv(rbind(a, b), "benchmarks/dev/results/srv36n/cells_all.csv", row.names = FALSE)'
Rscript benchmarks/dev/tools/null_summary.R "$R/$L/cells_all.csv" > "$R/$L/null_summary.txt" 2>&1; rm -f "$R/$L/cells_all.csv"
rm -rf "$R/$L/cells"; MSG="srv36n: global nulls for c15c (40 reps)" push "$R/$L"
L=srv36s; Rscript benchmarks/dev/spread_suite.R --candidates "$C" --settings H_unb,H_bal,H_all,H_cont --reps 20 --cores "$CORES" --label "$L" || exit 1
rm -rf "$R/$L/cells"; cp logs/run_srv36.log "$R/$L/run.log" 2>/dev/null || true
MSG="srv36s: spread suite for c15c" push "$R/$L"
echo ">> done $(date '+%F %T')"
