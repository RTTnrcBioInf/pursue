#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# run_axisDE.R -- drivers for the real-data axes with no simulated truth.
#
# Axis D (biological / experimental truth):
#   --axis D --template mbd_gingival_v35 --group body_subsite --expected expected_gingival.tsv
#   runs every method on the real grouping and scores: enrichment of the expected taxa
#   among calls (odds ratio, Fisher p), fraction of calls in the expected direction.
#   --axis D --template mbd_stammler_spikein --spikein "S1,S2,S3" runs 50 random binary
#   groupings and counts spike-in taxa called at q <= 0.05 (any call is a false positive).
#
# Axis E (replicability, Pelto et al. 2025):
#   --axis E --template crc_genus --group diagnosis --splits 5
#   for each of `splits` random 50/50 splits, runs every method on both halves and
#   reports replication% (same sign, q <= 0.05 in both), conflict% (opposite sign, both
#   significant) and NHits per method and arm.
# -----------------------------------------------------------------------------
suppressMessages({ library(optparse); library(jsonlite) })
op <- OptionParser(option_list = list(
  make_option("--axis", type = "character"), make_option("--template", type = "character"),
  make_option("--group", type = "character", default = NULL), make_option("--levels", type = "character", default = NULL),
  make_option("--expected", type = "character", default = NULL), make_option("--spikein", type = "character", default = NULL),
  make_option("--splits", type = "integer", default = 5L), make_option("--n-groupings", type = "integer", default = 50L),
  make_option("--methods", type = "character", default = "all"), make_option("--out", type = "character", default = "results"),
  make_option("--tag", type = "character", default = "", help = "suffix for output files when re-running a subset of methods"),
  make_option("--master-seed", type = "integer", default = 1L), make_option("--min-prevalence", type = "double", default = 0.10),
  make_option("--bench-root", type = "character", default = NULL), make_option("--n-cores", type = "integer", default = 1L)))
opt <- parse_args(op)
if (!is.null(opt$`bench-root`)) Sys.setenv(PURSUE_BENCH_ROOT = opt$`bench-root`)
root <- Sys.getenv("PURSUE_BENCH_ROOT", unset = ".")
for (f in c("R/engine/templates.R", "R/engine/regimes.R", "R/engine/metrics.R", "R/engine/io.R",
            "R/methods/elementary.R", "R/methods/external.R", "R/methods/registry.R", "R/methods/pursue03.R")) source(file.path(root, f))
methods <- if (opt$methods == "all") method_registry()$id else strsplit(opt$methods, ",")[[1]]
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
tagsfx <- if (nzchar(opt$tag)) paste0("__", opt$tag) else ""
set.seed(opt$`master-seed`)
tpl <- load_template(opt$template)
gv <- if (!is.null(opt$group)) opt$group else tpl$group_var

prep <- function(ct, meta) { keep <- rowMeans(ct > 0) >= opt$`min-prevalence` & rowSums(ct > 0) >= 3; list(counts = ct[keep, , drop = FALSE], meta = meta) }
binary_meta <- function(meta, gv, levels = NULL) {
  v <- as.character(meta[[gv]]); lv <- if (!is.null(levels)) levels else sort(unique(v))[1:2]
  keep <- v %in% lv; data.frame(row.names = rownames(meta)[keep], group = factor(v[keep], levels = lv)) }

if (opt$axis == "D" && is.null(opt$spikein)) {
  lv <- if (!is.null(opt$levels)) strsplit(opt$levels, ",")[[1]] else NULL
  md <- binary_meta(tpl$meta, gv, lv); ct <- tpl$counts[, rownames(md), drop = FALSE]; d <- prep(ct, md)
  exp <- if (!is.null(opt$expected)) read.delim(file.path(root, opt$expected), stringsAsFactors = FALSE) else NULL
  rows <- list()
  for (mid in methods) {
    run <- run_method(mid, d$counts, d$meta, ~ group, "group", args = if (mid == "pursue") list(n_cores = opt$`n-cores`) else list())
    r <- run$result
    for (arm in unique(r$arm)) { ra <- r[r$arm == arm, ]; sig <- !is.na(ra$q) & ra$q <= 0.05
      row <- data.frame(axis = "D", template = opt$template, method = mid, arm = arm, n_tested = sum(is.finite(ra$p)), n_calls = sum(sig), runtime_s = run$runtime_s)
      if (!is.null(exp)) {
        e <- exp[match(ra$feature, exp$feature), ]; known <- !is.na(e$expected_direction)
        tab <- table(called = sig[known], expected = known[known] & TRUE)
        inexp <- sig & known; row$n_calls_annotated <- sum(inexp)
        row$frac_expected_direction <- if (any(inexp & is.finite(ra$estimate))) mean(sign(ra$estimate[inexp & is.finite(ra$estimate)]) == e$expected_direction[inexp & is.finite(ra$estimate)]) else NA
        row$enrichment_p <- tryCatch(stats::fisher.test(table(sig, known))$p.value, error = function(x) NA)
        row$enrichment_or <- tryCatch(stats::fisher.test(table(sig, known))$estimate, error = function(x) NA)
      }
      rows[[length(rows) + 1L]] <- row }
  }
  res <- do.call(rbind, rows); write.csv(res, file.path(opt$out, paste0("D__", opt$template, tagsfx, "__biotruth.csv")), row.names = FALSE); print(res, row.names = FALSE)
}

if (opt$axis == "D" && !is.null(opt$spikein)) {
  spk <- strsplit(opt$spikein, ",")[[1]]; ct <- tpl$counts; n <- ncol(ct)
  rows <- list()
  for (g in seq_len(opt$`n-groupings`)) {
    md <- data.frame(row.names = colnames(ct), group = factor(sample(rep(c("control", "case"), length.out = n)), levels = c("control", "case")))
    d <- prep(ct, md)
    for (mid in methods) { run <- run_method(mid, d$counts, d$meta, ~ group, "group"); r <- run$result
      for (arm in unique(r$arm)) { ra <- r[r$arm == arm, ]; sig <- !is.na(ra$q) & ra$q <= 0.05
        rows[[length(rows) + 1L]] <- data.frame(axis = "D", template = opt$template, grouping = g, method = mid, arm = arm,
          n_calls = sum(sig), spikein_calls = sum(sig & ra$feature %in% spk), spikein_tested = sum(ra$feature %in% spk & is.finite(ra$p)), runtime_s = run$runtime_s) } }
    cat("grouping", g, "done\n")
  }
  res <- do.call(rbind, rows); write.csv(res, file.path(opt$out, paste0("D__", opt$template, tagsfx, "__spikein.csv")), row.names = FALSE)
  print(aggregate(cbind(n_calls, spikein_calls) ~ method + arm, res, mean), row.names = FALSE)
}

if (opt$axis == "E") {
  lv <- if (!is.null(opt$levels)) strsplit(opt$levels, ",")[[1]] else NULL
  md <- binary_meta(tpl$meta, gv, lv); ct <- tpl$counts[, rownames(md), drop = FALSE]; n <- ncol(ct)
  rows <- list()
  for (s in seq_len(opt$splits)) {
    idx <- sample(n); h1 <- sort(idx[seq_len(floor(n / 2))]); h2 <- sort(idx[-seq_len(floor(n / 2))])
    d1 <- prep(ct[, h1], md[h1, , drop = FALSE]); d2 <- prep(ct[, h2], md[h2, , drop = FALSE])
    common <- intersect(rownames(d1$counts), rownames(d2$counts))
    for (mid in methods) {
      r1 <- run_method(mid, d1$counts[common, ], d1$meta, ~ group, "group")$result; r2 <- run_method(mid, d2$counts[common, ], d2$meta, ~ group, "group")$result
      for (arm in unique(r1$arm)) {
        a <- r1[r1$arm == arm, ]; b <- r2[r2$arm == arm, ]; b <- b[match(a$feature, b$feature), ]
        s1 <- !is.na(a$q) & a$q <= 0.05; s2 <- !is.na(b$q) & b$q <= 0.05
        dir_a <- sign(if (all(is.na(a$estimate))) -log(a$p) * 0 + 1 else a$estimate); dir_b <- sign(if (all(is.na(b$estimate))) -log(b$p) * 0 + 1 else b$estimate)
        hits <- sum(s1) + sum(s2); both <- s1 & s2
        rep_ <- sum(both & dir_a == dir_b); conf <- sum(both & dir_a != dir_b & !is.na(dir_a) & !is.na(dir_b))
        rows[[length(rows) + 1L]] <- data.frame(axis = "E", template = opt$template, split = s, method = mid, arm = arm, n_features = length(common),
          nhits = hits, replication_pct = if (hits) 100 * 2 * rep_ / hits else NA, conflict_pct = if (hits) 100 * 2 * conf / hits else NA)
      }
    }
    cat("split", s, "done\n")
  }
  res <- do.call(rbind, rows); write.csv(res, file.path(opt$out, paste0("E__", opt$template, tagsfx, "__replicability.csv")), row.names = FALSE)
  print(aggregate(cbind(nhits, replication_pct, conflict_pct) ~ method + arm, res, function(x) mean(x, na.rm = TRUE)), row.names = FALSE)
}
