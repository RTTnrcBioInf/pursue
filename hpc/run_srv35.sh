#!/usr/bin/env bash
# run_srv35.sh -- after conf2 (2026-11-04): (1) are the comparators conf1 / conf2 left out of A-C, at their CURRENT versions,
# ahead of c14 anywhere on the tuning templates? (2) does conf2's largest remaining gap -- SimulateMSeq against LDM, widest at
# 100 / 200 per group -- exist on the tuning templates? (3) bank cells for local work on that gap and on the compositional
# centre (one-directional designs). Development data only (tuning templates, dev seed 9001); nothing held out is touched.
#   0. install LOCOM2 (CRAN, current: 1.1) and MaAsLin 3 (GitHub biobakery/maaslin3 at d5f36fb, version 1.5.8, the newest; the
#      benchmark ran 1.0 and 1.2.0) into a separate library, NEWLIB (default ~/R/pursue-srv35); the benchmark environment is
#      not changed. Only the MaAsLin 3 / LOCOM2 runs put NEWLIB first on the library path. Versions -> srv35m/versions.txt.
#   1. bank5: msq R00 / R02 / R03 / R10 / R13 / R06 (null), mid R10 / R11; 5 reps x both templates; c14, LDM and ZicoSeq p.
#   2. srv35: the whole dev suite (house, implants, stress, msq, mid, sd2, sps) + msq R02 / R03 / R09 / R10 / R13 / R22, mid R02 /
#      R10 / R13, house R10, implant B_cont / B_conf04 / B_conf07; 10 reps; c14 and LDM (current environment).
#   3. srv35m: the same cells, MaAsLin 3 (combined arm, as the benchmark scores it) and LOCOM2 (NEWLIB).
#   4. srv35n / srv35mn: global nulls (house / msq / mid / sd2 / sps R06, implant B_null), 40 reps, the same four methods.
# Cells are identical across the runs (same seeds, single-threaded BLAS), so compare.R merges them. ~90 CPU-hours: ~8 h at
# 12 cores -- an underestimate: SimulateMSeq refitted its template for every cell (~25 min per twinsuk_stool cell). Since
# 2026-11-04 that fit is cached (benchmarks/dev/tools/msq_cache.R; step 1b checks the cells are unchanged), so a relaunch
# resumes without it. hpc/guard.sh: shared memory cap (MEMTOTAL), OOM-killer first choice, single-threaded BLAS. Resumable: finished
# cells are kept, so relaunching continues. Each stage pushes.
#   git pull && mkdir -p logs && (MEMTOTAL=90G CORES=12 nohup nice -n 10 bash hpc/run_srv35.sh > logs/run_srv35.log 2>&1 &)
source "$(dirname "${BASH_SOURCE[0]}")/guard.sh" "$@"
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; mkdir -p logs benchmarks/dev/results benchmarks/dev/bank
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
CORES="${CORES:-12}"; NEWLIB="${NEWLIB:-$HOME/R/pursue-srv35}"; R=benchmarks/dev/results
M3=d5f36fbc0a6a77001f0736ae389618fba07b1e1f
SUITE="--suite full --sims house,implant,msq,mid,sd2,sps"
EXTRA="msq:R02,msq:R03,msq:R09,msq:R10,msq:R13,msq:R22,mid:R02,mid:R10,mid:R13,house:R10,implant:B_cont,implant:B_conf04,implant:B_conf07"
NULLS="house:R06,msq:R06,mid:R06,sd2:R06,sps:R06,implant:B_null"
push() { git add "$@" && git commit -q -m "$MSG" && git pull -q --rebase --autostash && git push -q; }
echo ">> srv35 ($(date '+%F %T')), $CORES cores, NEWLIB $NEWLIB"

# 0. current LOCOM2 and MaAsLin 3 in their own library
mkdir -p "$NEWLIB" "$R/srv35m"
NEWLIB="$NEWLIB" M3="$M3" Rscript -e '
  lib <- Sys.getenv("NEWLIB"); .libPaths(c(lib, .libPaths())); options(repos = c(CRAN = "https://cloud.r-project.org"))
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", lib = lib)
  if (!requireNamespace("LOCOM2", lib.loc = lib, quietly = TRUE)) install.packages("LOCOM2", lib = lib)
  if (!requireNamespace("maaslin3", lib.loc = lib, quietly = TRUE)) {
    ok <- tryCatch({ remotes::install_github(paste0("biobakery/maaslin3@", Sys.getenv("M3")), lib = lib, upgrade = "never", dependencies = NA); TRUE }, error = function(e) FALSE)
    if (!ok || !requireNamespace("maaslin3", lib.loc = lib, quietly = TRUE))   # a dependency below the version it needs: newer copies into lib only
      remotes::install_github(paste0("biobakery/maaslin3@", Sys.getenv("M3")), lib = lib, upgrade = "always", dependencies = NA) }
  v <- function(p, l = NULL) tryCatch(as.character(utils::packageVersion(p, lib.loc = l)), error = function(e) "not installed")
  writeLines(c(sprintf("new library %s: LOCOM2 %s, maaslin3 %s (GitHub %s)", lib, v("LOCOM2", lib), v("maaslin3", lib), substr(Sys.getenv("M3"), 1, 7)),
    sprintf("benchmark environment: LOCOM2 %s, maaslin3 %s, LDM %s, R %s", v("LOCOM2", .libPaths()[-1]), v("maaslin3", .libPaths()[-1]), v("LDM"), getRversion())),
    file.path("benchmarks/dev/results/srv35m", "versions.txt"))' > logs/run_srv35_install.log 2>&1
cat "$R/srv35m/versions.txt" 2>/dev/null || { echo ">> install step failed -- see logs/run_srv35_install.log"; exit 1; }
grep -q "not installed" <(head -1 "$R/srv35m/versions.txt") && { echo ">> install FAILED -- see logs/run_srv35_install.log"; exit 1; }
SM="srv35_smoke"
R_LIBS="$NEWLIB" Rscript benchmarks/dev/devsuite.R --settings house:R00 --templates hmp_tongue --reps 1 --cores 1 --label "$SM" \
  --candidates ref_maaslin3,ref_locom2 > logs/run_srv35_smoke.log 2>&1
Rscript -e 'd <- read.csv(file.path("benchmarks/dev/results", "srv35_smoke", "cells.csv")); print(d[, c("candidate", "tested", "rej", "error")]); if (any(!is.na(d$error)) || any(d$tested == 0)) quit(status = 1)' \
  || { echo ">> smoke FAILED for MaAsLin 3 / LOCOM2 -- see logs/run_srv35_smoke.log"; exit 1; }
rm -rf "$R/$SM"

# 1. bank5 -- cells for local work (SimulateMSeq vs LDM; the centre under one-directional change)
B=benchmarks/dev/bank/bank5
echo ">> bank $B ($(date '+%F %T'))"
Rscript benchmarks/dev/tools/bank_cells.R --settings msq:R00,msq:R02,msq:R03,msq:R10,msq:R13,msq:R06,mid:R10,mid:R11 --templates hmp_tongue,twinsuk_stool \
  --reps 5 --candidates eperm_c14,ref_ldm,ref_zicoseq --out "$B" --cores "$CORES" || exit 1
MSG="dev bank5: msq R00 / R02 / R03 / R10 / R13 / R06, mid R10 / R11 with c14, LDM, ZicoSeq p-values (tuning templates)" push "$B"

# 1b. the SimulateMSeq fit cache (added 2026-11-04, on relaunch): rebuild two of bank5's cells (made without the cache) with it
#     and stop unless they are identical -- counts, truth and metadata. The bank step above has already fitted the templates.
CK="$(mktemp -d)"
Rscript benchmarks/dev/tools/bank_cells.R --settings msq:R00 --templates hmp_tongue,twinsuk_stool --reps 1 --candidates lm_logtss \
  --out "$CK" --cores 2 > logs/run_srv35_msqcheck.log 2>&1
CK="$CK" B="$B" Rscript -e 'ck <- Sys.getenv("CK"); b <- Sys.getenv("B"); f <- list.files(ck, "[.]rds$")
  same <- vapply(f, function(x) { u <- readRDS(file.path(ck, x)); v <- readRDS(file.path(b, x))
    identical(u$counts, v$counts) && identical(u$truth, v$truth) && identical(u$meta, v$meta) }, logical(1))
  print(same); quit(status = if (length(f) == 2L && all(same)) 0 else 1)' \
  || { echo ">> msq cache check FAILED -- cells differ from bank5 (see logs/run_srv35_msqcheck.log); rerun with PURSUE_MSQ_CACHE=0"; exit 1; }
rm -rf "$CK"; echo ">> msq cache check passed: cells identical to bank5 ($(date '+%F %T'))"

# 2-3. the dev suite: c14 and LDM in the benchmark environment; MaAsLin 3 and LOCOM2 from NEWLIB
L=srv35; echo ">> dev suite $L ($(date '+%F %T'))"
Rscript benchmarks/dev/devsuite.R $SUITE --extra "$EXTRA" --reps 10 --cores "$CORES" --label "$L" --candidates eperm_c14,ref_ldm || exit 1
rm -rf "$R/$L/cells"; MSG="srv35: c14 and LDM on the whole dev suite + larger-n msq settings (tuning templates)" push "$R/$L"
L=srv35m; echo ">> dev suite $L, new library ($(date '+%F %T'))"
R_LIBS="$NEWLIB" Rscript benchmarks/dev/devsuite.R $SUITE --extra "$EXTRA" --reps 10 --cores "$CORES" --label "$L" --candidates ref_maaslin3,ref_locom2 || exit 1
Rscript benchmarks/dev/compare.R --labels "srv35,srv35m" > "$R/$L/compare.txt" 2>&1
rm -rf "$R/$L/cells"; MSG="srv35m: MaAsLin 3 and LOCOM2 (current versions) on the same dev-suite cells" push "$R/$L"

# 4. global nulls at 40 reps
L=srv35n; echo ">> nulls $L ($(date '+%F %T'))"
Rscript benchmarks/dev/devsuite.R $SUITE --settings "$NULLS" --reps 40 --cores "$CORES" --label "$L" --candidates eperm_c14,ref_ldm || exit 1
L=srv35mn; R_LIBS="$NEWLIB" Rscript benchmarks/dev/devsuite.R $SUITE --settings "$NULLS" --reps 40 --cores "$CORES" --label "$L" --candidates ref_maaslin3,ref_locom2 || exit 1
Rscript -e 'a <- read.csv("benchmarks/dev/results/srv35n/cells.csv"); b <- read.csv("benchmarks/dev/results/srv35mn/cells.csv"); write.csv(rbind(a, b), "benchmarks/dev/results/srv35mn/cells_all.csv", row.names = FALSE)'
Rscript benchmarks/dev/tools/null_summary.R "$R/srv35mn/cells_all.csv" > "$R/srv35mn/null_summary.txt" 2>&1
rm -f "$R/srv35mn/cells_all.csv"; rm -rf "$R/srv35n/cells" "$R/srv35mn/cells"; cp logs/run_srv35.log "$R/srv35mn/run.log" 2>/dev/null || true
MSG="srv35n / srv35mn: global nulls at 40 reps for c14, LDM, MaAsLin 3, LOCOM2 (tuning templates)" push "$R/srv35n" "$R/srv35mn"
echo ">> done $(date '+%F %T')"
