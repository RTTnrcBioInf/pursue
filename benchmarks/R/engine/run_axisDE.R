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
#   --axis E --template crc_genus --group Diagnosis --levels Normal,Cancer --splits 5
#   for each of `splits` random 50/50 splits, runs every method on both halves and
#   reports replication% (same sign, q <= 0.05 in both), conflict% (opposite sign, both
#   significant) and NHits per method and arm.
#   --axis E --template A --template2 B --group g --levels c,k   (conf2, protocol 13): cross-study replication -- every
#   method on all of A and all of B, the same replication / conflict counts between the two.
# Axis D, measured load (conf2, protocol 13):
#   --axis D --qmp --template d2_vandeputte --group group --levels healthy,CD --case CD
#   truth on the load-scaled scale: abundance = count / depth x cells (meta$cell_count); per taxon the load-scaled log2
#   fold change (difference in group means of log2(abundance + c), c = half the taxon's smallest positive value) and a
#   reference set (BH q <= 0.05, Wilcoxon on load-scaled abundance). Scored per method: sign agreement of its calls with
#   the load-scaled fold change, recall of the reference set, calls outside it, Spearman of its estimates with the
#   load-scaled fold changes, and its direction on the taxa whose relative and load-scaled directions disagree (reference
#   taxa only: disagree_*; any tested taxon: disagree_any_*).
# -----------------------------------------------------------------------------
suppressMessages({ library(optparse); library(jsonlite) })
op <- OptionParser(option_list = list(
  make_option("--axis", type = "character"), make_option("--template", type = "character"),
  make_option("--group", type = "character", default = NULL), make_option("--levels", type = "character", default = NULL),
  make_option("--expected", type = "character", default = NULL), make_option("--spikein", type = "character", default = NULL),
  make_option("--case", type = "character", default = NULL, help = "regex for the level that expected_direction +1 refers to; it is made the tested (second) level"),
  make_option("--splits", type = "integer", default = 5L), make_option("--n-groupings", type = "integer", default = 50L),
  make_option("--methods", type = "character", default = "all"), make_option("--out", type = "character", default = "results"),
  make_option("--tag", type = "character", default = "", help = "suffix for output files when re-running a subset of methods"),
  make_option("--master-seed", type = "integer", default = 1L), make_option("--min-prevalence", type = "double", default = 0.10),
  make_option("--bench-root", type = "character", default = NULL), make_option("--n-cores", type = "integer", default = 1L),
  make_option("--template2", type = "character", default = NULL, help = "axis E: a second dataset for cross-study replication"),
  make_option("--qmp", action = "store_true", default = FALSE, help = "axis D: truth from measured load (meta$cell_count)")))
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
binary_meta <- function(meta, gv, levels = NULL, case = NULL) {
  if (!gv %in% names(meta)) { hit <- names(meta)[tolower(names(meta)) == tolower(gv)]      # e.g. crc metadata.csv has "Diagnosis"
    if (length(hit) != 1L) stop("group column '", gv, "' is not in the metadata; columns: ", paste(names(meta), collapse = ", "))
    gv <- hit }
  v <- as.character(meta[[gv]])
  if (!is.null(levels)) levels <- vapply(levels, function(l) { if (l %in% v) return(l)            # exact level, else a unique case-insensitive match
    h <- unique(v[!is.na(v) & tolower(v) == tolower(l)]); if (length(h) == 1L) h else l }, "", USE.NAMES = FALSE)
  lv <- if (!is.null(levels)) levels else sort(unique(v[!is.na(v)]))[1:2]
  if (!is.null(case)) { m <- grepl(case, lv, ignore.case = TRUE)
    if (sum(m) != 1L) stop("--case '", case, "' must match exactly one of the levels: ", paste(lv, collapse = ", "))
    lv <- c(lv[!m], lv[m]) }
  keep <- v %in% lv
  if (sum(v == lv[1], na.rm = TRUE) < 5L || sum(v == lv[2], na.rm = TRUE) < 5L)
    stop("group '", gv, "' levels ", paste(lv, collapse = " / "), " have too few samples; available: ", paste(unique(v), collapse = ", "))
  cat(">> contrast:", lv[2], "vs", lv[1], sprintf("(%d vs %d samples)\n", sum(v == lv[2], na.rm = TRUE), sum(v == lv[1], na.rm = TRUE)))
  data.frame(row.names = rownames(meta)[keep], group = factor(v[keep], levels = lv)) }

if (opt$axis == "D" && isTRUE(opt$qmp)) {
  lv <- if (!is.null(opt$levels)) strsplit(opt$levels, ",")[[1]] else NULL
  md <- binary_meta(tpl$meta, gv, lv, opt$case); ct <- tpl$counts[, rownames(md), drop = FALSE]; d <- prep(ct, md)
  cells <- tpl$meta[rownames(md), "cell_count"]; if (is.null(cells) || anyNA(cells)) stop("--qmp needs meta$cell_count for every sample")
  X <- d$counts; N <- colSums(ct); g <- md$group == levels(md$group)[2]
  A <- sweep(X, 2L, cells / N, "*"); Rl <- sweep(X, 2L, N, "/")
  lfc <- function(M) apply(M, 1L, function(v) { c0 <- min(v[v > 0]) / 2; mean(log2(v[g] + c0)) - mean(log2(v[!g] + c0)) })
  lfc_abs <- lfc(A); lfc_rel <- lfc(Rl)
  q_abs <- stats::p.adjust(apply(A, 1L, function(v) stats::wilcox.test(v[g], v[!g], exact = FALSE)$p.value), "BH")
  disall <- sign(lfc_rel) != sign(lfc_abs); ref <- q_abs <= 0.05; dis <- ref & disall
  tf <- file.path(opt$out, paste0("truth__", opt$template, "__qmp.csv"))
  if (!file.exists(tf)) { tmp <- paste0(tf, ".tmp", Sys.getpid())     # several methods run in parallel: write once, atomically
    write.csv(data.frame(feature = rownames(X), lfc_abs = lfc_abs, lfc_rel = lfc_rel, q_abs = q_abs, reference = ref, disagree = dis, disagree_any = disall), tmp, row.names = FALSE)
    file.rename(tmp, tf) }
  cat(sprintf(">> load-scaled truth: %d taxa tested, %d in the reference set, %d of them with opposite relative direction\n", nrow(X), sum(ref), sum(dis)))
  est_class <- method_registry()$estimand; names(est_class) <- method_registry()$id
  rows <- list()   # conf2 deviation (protocol 13, 2026-10-28): Spearman is NA for an arm without estimates (it stopped the job)
  for (mid in methods) {
    run <- run_method(mid, d$counts, d$meta, ~ group, "group"); r <- run$result
    for (arm in unique(r$arm)) { ra <- r[r$arm == arm, ]; i <- match(rownames(X), ra$feature); sig <- !is.na(ra$q[i]) & ra$q[i] <= 0.05
      est <- ra$estimate[i]; se <- sig & is.finite(est); sd_ <- se & dis; sa_ <- se & disall
      rows[[length(rows) + 1L]] <- data.frame(axis = "D", template = opt$template, contrast = paste(rev(levels(md$group)), collapse = " vs "),
        method = mid, estimand = unname(est_class[mid]), arm = arm, n_tested = sum(is.finite(ra$p[i])), n_calls = sum(sig),
        sign_agreement = if (any(se)) mean(sign(est[se]) == sign(lfc_abs[se])) else NA, recall_reference = sum(sig & ref) / max(sum(ref), 1),
        calls_outside_reference = sum(sig & !ref), spearman_lfc_abs = if (sum(is.finite(est) & is.finite(lfc_abs)) >= 3L) suppressWarnings(stats::cor(est, lfc_abs, method = "spearman", use = "complete.obs")) else NA,
        disagree_called = sum(sig & dis), disagree_absolute_direction = if (any(sd_)) mean(sign(est[sd_]) == sign(lfc_abs[sd_])) else NA,
        disagree_any_called = sum(sig & disall), disagree_any_absolute_direction = if (any(sa_)) mean(sign(est[sa_]) == sign(lfc_abs[sa_])) else NA,
        runtime_s = run$runtime_s) }
  }
  res <- do.call(rbind, rows); write.csv(res, file.path(opt$out, paste0("D__", opt$template, tagsfx, "__qmp.csv")), row.names = FALSE); print(res, row.names = FALSE)
}

if (opt$axis == "D" && is.null(opt$spikein) && !isTRUE(opt$qmp)) {
  lv <- if (!is.null(opt$levels)) strsplit(opt$levels, ",")[[1]] else NULL
  md <- binary_meta(tpl$meta, gv, lv, opt$case); ct <- tpl$counts[, rownames(md), drop = FALSE]; d <- prep(ct, md)
  contrast <- paste(rev(levels(md$group)), collapse = " vs ")
  exp <- if (!is.null(opt$expected)) read.delim(file.path(root, opt$expected), stringsAsFactors = FALSE) else NULL
  rows <- list()
  for (mid in methods) {
    run <- run_method(mid, d$counts, d$meta, ~ group, "group", args = if (mid == "pursue") list(n_cores = opt$`n-cores`) else list())
    r <- run$result
    for (arm in unique(r$arm)) { ra <- r[r$arm == arm, ]; sig <- !is.na(ra$q) & ra$q <= 0.05
      row <- data.frame(axis = "D", template = opt$template, contrast = contrast, method = mid, arm = arm, n_tested = sum(is.finite(ra$p)), n_calls = sum(sig), runtime_s = run$runtime_s)
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
  grp <- lapply(seq_len(opt$`n-groupings`), function(g) sample(rep(c("control", "case"), length.out = n)))   # drawn before any method runs
  rows <- list()
  for (g in seq_len(opt$`n-groupings`)) {
    md <- data.frame(row.names = colnames(ct), group = factor(grp[[g]], levels = c("control", "case")))
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

if (opt$axis == "E" && !is.null(opt$template2)) {
  lv <- if (!is.null(opt$levels)) strsplit(opt$levels, ",")[[1]] else NULL
  tpl2 <- load_template(opt$template2)
  m1 <- binary_meta(tpl$meta, gv, lv); m2 <- binary_meta(tpl2$meta, gv, lv)
  d1 <- prep(tpl$counts[, rownames(m1), drop = FALSE], m1); d2 <- prep(tpl2$counts[, rownames(m2), drop = FALSE], m2)
  common <- intersect(rownames(d1$counts), rownames(d2$counts)); cat(">> taxa tested in both:", length(common), "\n")
  rows <- list()
  for (mid in methods) {
    r1 <- run_method(mid, d1$counts[common, ], d1$meta, ~ group, "group")$result; r2 <- run_method(mid, d2$counts[common, ], d2$meta, ~ group, "group")$result
    for (arm in unique(r1$arm)) {
      a <- r1[r1$arm == arm, ]; b <- r2[r2$arm == arm, ]; b <- b[match(a$feature, b$feature), ]
      s1 <- !is.na(a$q) & a$q <= 0.05; s2 <- !is.na(b$q) & b$q <= 0.05
      dir_a <- sign(if (all(is.na(a$estimate))) -log(a$p) * 0 + 1 else a$estimate); dir_b <- sign(if (all(is.na(b$estimate))) -log(b$p) * 0 + 1 else b$estimate)
      hits <- sum(s1) + sum(s2); both <- s1 & s2
      rep_ <- sum(both & dir_a == dir_b); conf <- sum(both & dir_a != dir_b & !is.na(dir_a) & !is.na(dir_b))
      rows[[length(rows) + 1L]] <- data.frame(axis = "E", template = paste0(opt$template, "-vs-", opt$template2), split = 0L, method = mid, arm = arm,
        n_features = length(common), hits_1 = sum(s1), hits_2 = sum(s2), nhits = hits,
        replication_pct = if (hits) 100 * 2 * rep_ / hits else NA, conflict_pct = if (hits) 100 * 2 * conf / hits else NA)
    }
  }
  res <- do.call(rbind, rows)
  write.csv(res, file.path(opt$out, paste0("E__", opt$template, "-vs-", opt$template2, tagsfx, "__crossstudy.csv")), row.names = FALSE); print(res, row.names = FALSE)
}

if (opt$axis == "E" && is.null(opt$template2)) {
  lv <- if (!is.null(opt$levels)) strsplit(opt$levels, ",")[[1]] else NULL
  md <- binary_meta(tpl$meta, gv, lv); ct <- tpl$counts[, rownames(md), drop = FALSE]; n <- ncol(ct)
  perms <- lapply(seq_len(opt$splits), function(s) sample(n))                   # drawn before any method runs
  rows <- list()
  for (s in seq_len(opt$splits)) {
    idx <- perms[[s]]; h1 <- sort(idx[seq_len(floor(n / 2))]); h2 <- sort(idx[-seq_len(floor(n / 2))])
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
