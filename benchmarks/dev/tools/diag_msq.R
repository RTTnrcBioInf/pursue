#!/usr/bin/env Rscript
# diag_msq.R -- per-feature comparison of the PURSUE 0.3 candidate with ZicoSeq and LDM on dev-suite msq
# cells (same seeds as the dev suite), to see which true DA features the comparators find and we miss.
#   Rscript benchmarks/dev/tools/diag_msq.R --cores 10
# Writes benchmarks/dev/results/diag_msq/features.csv (one row per feature x cell).
suppressMessages(library(optparse))
opt <- parse_args(OptionParser(option_list = list(make_option("--cores", type = "integer", default = 10L),
  make_option("--settings", default = "msq:R00,msq:R17"), make_option("--reps", default = "1,2,3,4,5"), make_option("--templates", default = "hmp_tongue,twinsuk_stool"))))
Sys.setenv(PURSUE_BENCH_ROOT = normalizePath("benchmarks"))
suppressMessages(library(PURSUE)); source("benchmarks/dev/tools/cell.R")
grid <- expand.grid(sid = strsplit(opt$settings, ",")[[1]], tid = strsplit(opt$templates, ",")[[1]], rep = as.integer(strsplit(opt$reps, ",")[[1]]), stringsAsFactors = FALSE)
one <- function(k) {
  g <- grid[k, ]; cl <- suppressMessages(dev_cell(g$sid, g$tid, g$rep)); ct <- cl$counts; dep <- cl$meta$depth; g1 <- cl$meta$group == "case"
  v <- .efull_fit(ct, cl$meta, cl$formula, cl$tested_term, sf_balanced_only = TRUE, opts3 = TRUE, winsor = 0.03)
  zs <- tryCatch(.cands$ref_zicoseq$fn(ct, cl$meta, cl$formula, cl$tested_term), error = function(e) NULL)
  ld <- tryCatch(.cands$ref_ldm$fn(ct, cl$meta, cl$formula, cl$tested_term), error = function(e) NULL)
  rel <- sweep(ct, 2, dep, "/"); m <- function(r, col) if (is.null(r)) NA else r[[col]][match(rownames(ct), r$feature)]
  data.frame(setting = g$sid, template = g$tid, rep = g$rep, feature = rownames(ct), truth = cl$truth$truth_abs, truth_rel_lfc2 = cl$truth$truth_rel_lfc2,
    prev0 = rowMeans(ct[, !g1] > 0), prev1 = rowMeans(ct[, g1] > 0), mrel = rowMeans(rel), mcount = rowMeans(ct),
    sd_logrel_pos = apply(rel, 1, function(x) if (sum(x > 0) > 2) sd(log(x[x > 0])) else NA),
    top3_share = apply(rel, 1, function(x) { s <- sort(x, decreasing = TRUE); sum(s[1:3]) / max(sum(s), 1e-12) }),
    p_ours = v$p_full, p_det = v$p_det, p_log = v$p_log, p_sqrt = v$p_sqrt, p_zs = m(zs, "p"), q_zs = m(zs, "q"), p_ldm = m(ld, "p"), q_ldm = m(ld, "q"),
    gam = v$gam, depth_cv = sd(dep) / mean(dep), stringsAsFactors = FALSE)
}
res <- parallel::mclapply(seq_len(nrow(grid)), function(k) tryCatch(one(k), error = function(e) { message("cell ", k, ": ", conditionMessage(e)); NULL }), mc.cores = opt$cores)
out <- do.call(rbind, res); dir.create("benchmarks/dev/results/diag_msq", recursive = TRUE, showWarnings = FALSE)
write.csv(out, "benchmarks/dev/results/diag_msq/features.csv", row.names = FALSE)
q <- function(p) ave(p, out$setting, out$template, out$rep, FUN = function(z) p.adjust(z, "BH"))
cat(sprintf("TP: ours %d  zicoseq %d  ldm %d  (FP %d / %d / %d)\n", sum(q(out$p_ours) <= .05 & out$truth == 1, na.rm = TRUE), sum(out$q_zs <= .05 & out$truth == 1, na.rm = TRUE),
  sum(out$q_ldm <= .05 & out$truth == 1, na.rm = TRUE), sum(q(out$p_ours) <= .05 & out$truth == 0, na.rm = TRUE), sum(out$q_zs <= .05 & out$truth == 0, na.rm = TRUE), sum(out$q_ldm <= .05 & out$truth == 0, na.rm = TRUE)))
