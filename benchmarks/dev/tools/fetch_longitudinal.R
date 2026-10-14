#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# fetch_longitudinal.R -- public 16S datasets with several samples per subject, for the repeated-measures work
# (2026-10-14): real-data nulls with subject-level labels, subject-level implants, and paired (within-subject) designs.
#
# Pulls from Bioconductor data packages (ExperimentHub; needs internet, so it runs on the server) and writes one
# template per dataset (or body site) in the templates' shape -- list(counts = taxa x samples integer matrix, meta =
# data.frame, id, subject_var, time_var) -- to benchmarks/data/longitudinal/<id>.rds. Taxa present in < 5% of samples
# are dropped to keep files small. Nothing is analysed. The summary (samples, subjects, samples per subject, metadata
# columns) goes to benchmarks/dev/results/longi_fetch/summary.txt so it can be read without the data.
#
# Data discipline: HMP V35 sites overlap the evaluation templates (hmp_stool, hmp_gingiva, hmp_skin_ear, hmp_vagina
# come from the same HMP V35 samples). Only the tongue (a tuning template) may be used for development; the summary
# marks the others. The HMP2 cohorts and the baboon series are new to the project.
# -----------------------------------------------------------------------------
options(timeout = 3600, warn = 1)
root <- normalizePath(Sys.getenv("PURSUE_BENCH_ROOT", unset = "benchmarks"), mustWork = TRUE)
out <- file.path(root, "data", "longitudinal"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
rdir <- file.path(root, "dev", "results", "longi_fetch"); dir.create(rdir, recursive = TRUE, showWarnings = FALSE)
log <- file.path(rdir, "summary.txt"); cat(format(Sys.time()), "fetch_longitudinal\n\n", file = log)
say <- function(...) { s <- paste0(...); cat(s, "\n"); cat(s, "\n", file = log, append = TRUE) }

pk <- c("HMP16SData", "HMP2Data", "microbiomeDataSets", "phyloseq", "SummarizedExperiment")
need <- pk[!vapply(pk, requireNamespace, logical(1), quietly = TRUE)]
if (length(need)) { if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
  BiocManager::install(need, ask = FALSE, update = FALSE) }
say("packages: ", paste(sprintf("%s %s", pk, vapply(pk, function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) "MISSING"), "")), collapse = ", "))

guess <- function(cols, pat) { h <- grep(pat, cols, ignore.case = TRUE, value = TRUE); if (length(h)) h[1] else NA_character_ }
save_tpl <- function(id, counts, meta, subject_var = NA, time_var = NA, note = "") {
  counts <- as.matrix(counts); storage.mode(counts) <- "integer"; meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  stopifnot(ncol(counts) == nrow(meta)); rownames(meta) <- colnames(counts)
  keep <- rowMeans(counts > 0) >= 0.05; counts <- counts[keep, , drop = FALSE]
  ok <- colSums(counts) > 0; counts <- counts[, ok, drop = FALSE]; meta <- meta[ok, , drop = FALSE]
  cols <- names(meta)
  if (is.na(subject_var)) subject_var <- guess(cols, "^(rsid|subject|subject_id|host_subject_id|participant|individual|baboon|animal|host_id|patient)")
  if (is.na(time_var)) time_var <- guess(cols, "visit|week|day|time|collection_date|age")
  saveRDS(list(counts = counts, meta = meta, id = id, subject_var = subject_var, time_var = time_var), file.path(out, paste0(id, ".rds")), compress = "xz")
  sp <- if (!is.na(subject_var)) table(meta[[subject_var]]) else NULL
  say(sprintf("== %s: %d taxa (>= 5%% prevalence) x %d samples; median depth %.0f %s", id, nrow(counts), ncol(counts), stats::median(colSums(counts)), note))
  say("   subject column: ", subject_var, "; time column: ", time_var)
  if (!is.null(sp)) say(sprintf("   subjects %d; samples per subject: median %.0f, max %d; subjects with >= 2 samples %d, >= 5 samples %d",
                                length(sp), stats::median(sp), max(sp), sum(sp >= 2), sum(sp >= 5)))
  say("   metadata columns: ", paste(cols, collapse = ", "))
  bin <- cols[vapply(meta, function(v) { u <- unique(v[!is.na(v)]); length(u) >= 2 && length(u) <= 6 }, logical(1))]
  for (b in head(bin, 12)) say(sprintf("     %-28s %s", b, paste(utils::head(sprintf("%s=%d", names(table(meta[[b]])), as.integer(table(meta[[b]]))), 6), collapse = "  ")))
}
from_se <- function(se) list(counts = SummarizedExperiment::assay(se), meta = as.data.frame(SummarizedExperiment::colData(se)))
from_ps <- function(ps) { otu <- as(phyloseq::otu_table(ps), "matrix"); if (!phyloseq::taxa_are_rows(ps)) otu <- t(otu)
  list(counts = otu, meta = data.frame(phyloseq::sample_data(ps), check.names = FALSE)) }
try_ds <- function(label, expr) tryCatch(expr, error = function(e) say("!! ", label, " failed: ", conditionMessage(e)))

# 1. HMP16SData V35, one template per body subsite with enough repeat visits
try_ds("HMP16SData V35", {
  x <- from_se(HMP16SData::V35()); site <- x$meta$HMP_BODY_SUBSITE
  for (s in names(sort(table(site), decreasing = TRUE))) {
    i <- which(site == s); m <- x$meta[i, , drop = FALSE]; rep2 <- sum(table(m$RSID) >= 2)
    if (length(i) < 150 || rep2 < 40) next
    id <- paste0("hmp16s_", tolower(gsub("[^A-Za-z]+", "_", s)))
    note <- if (grepl("tongue", s, ignore.case = TRUE)) "(tuning site: usable for development)" else "(overlaps an HMP evaluation site or untouched: not for development)"
    save_tpl(id, x$counts[, i, drop = FALSE], m, "RSID", "VISITNO", note) } })

# 2. HMP2 cohorts (MOMS-PI vaginal, IBDMDB stool, T2D prediabetes); function names vary between versions, so try each
for (fn in c("momspi16S", "IBD16S", "T2D16S")) try_ds(paste0("HMP2Data ", fn), {
  f <- getExportedValue("HMP2Data", fn); obj <- f()
  x <- if (inherits(obj, "phyloseq")) from_ps(obj) else if (methods::is(obj, "SummarizedExperiment")) from_se(obj) else stop("unexpected class ", class(obj)[1])
  save_tpl(paste0("hmp2_", tolower(fn)), x$counts, x$meta, note = "(HMP2; new to the project)") })

# 3. microbiomeDataSets: wild baboon gut time series (Grieneisen et al. 2021)
try_ds("microbiomeDataSets GrieneisenTSData", {
  obj <- microbiomeDataSets::GrieneisenTSData()
  se <- if (methods::is(obj, "MultiAssayExperiment")) obj[[1]] else obj
  an <- SummarizedExperiment::assayNames(se); se_counts <- SummarizedExperiment::assay(se, if ("counts" %in% an) "counts" else an[1])
  save_tpl("baboon_grieneisen", se_counts, as.data.frame(SummarizedExperiment::colData(se)), note = "(wild baboons; new to the project)") })

say("\nfiles: ", paste(list.files(out), collapse = ", "))
