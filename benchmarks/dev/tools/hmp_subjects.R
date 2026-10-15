#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# hmp_subjects.R -- subject and visit labels for the HMP V35 samples (2026-10-15), so the MicrobeDS-based templates
# (hmp_tongue, the HMP evaluation templates) can be used for repeated-measures cells. Reads the HMP16SData sites that
# dev/tools/fetch_longitudinal.R saved (benchmarks/data/longitudinal/hmp16s_*.rds) and writes the sample map
# benchmarks/devdata/hmp16s_subjects.csv: SRS_SAMPLE_ID, RSID, VISITNO, SEX, RUN_CENTER, HMP_BODY_SUBSITE (public
# HMP metadata; no counts). Templates match it by the SRS id inside their X.SampleID. Prints how the tuning template
# hmp_tongue matches.
#   Rscript benchmarks/dev/tools/hmp_subjects.R
# -----------------------------------------------------------------------------
root <- normalizePath(Sys.getenv("PURSUE_BENCH_ROOT", unset = "benchmarks"), mustWork = TRUE)
f <- list.files(file.path(root, "data", "longitudinal"), pattern = "^hmp16s_.*[.]rds$", full.names = TRUE)
if (!length(f)) stop("no HMP16SData sites in benchmarks/data/longitudinal: run hpc/run_fetch_longi.sh first")
cols <- c("SRS_SAMPLE_ID", "RSID", "VISITNO", "SEX", "RUN_CENTER", "HMP_BODY_SUBSITE")
map <- do.call(rbind, lapply(f, function(x) { m <- readRDS(x)$meta; m <- m[, intersect(cols, names(m)), drop = FALSE]
  for (c in setdiff(cols, names(m))) m[[c]] <- NA; m[, cols] }))
map <- map[!is.na(map$SRS_SAMPLE_ID) & !duplicated(map$SRS_SAMPLE_ID), ]
out <- file.path(root, "devdata", "hmp16s_subjects.csv"); utils::write.csv(map, out, row.names = FALSE)
cat(sprintf("%d samples, %d subjects, %d sites -> %s\n", nrow(map), length(unique(map$RSID)), length(unique(map$HMP_BODY_SUBSITE)), out))
source(file.path(root, "dev", "tools", "rm.R"))
t <- rm_template(readRDS(file.path(root, "devdata", "hmp_tongue.rds")), map)
cat(sprintf("hmp_tongue: %d of 350 samples matched; %d after de-duplication; %d subjects, %d with >= 2 visits; visits %s\n",
            t$n_matched, t$n_samples, t$n_subjects, t$n_subj2, paste(names(table(t$meta$visit)), table(t$meta$visit), sep = "=", collapse = " ")))
