#!/usr/bin/env Rscript
# diag_fp.R -- what ARE a candidate's false positives? Rebuilds dev-suite cells from the same seeds
# (devsuite.R: dev seed 9001, same cell_seed) and describes every rejected feature whose absolute
# truth is null: realised relative fold change, detection shift between groups, prevalence,
# abundance. Answers whether "false" positives are compositional (relative change, absolute null).
#   Rscript benchmarks/dev/diag_fp.R --candidates sharedlink,logistic_depth --settings house:R00,implant:B_ref
suppressMessages(library(optparse))
op <- OptionParser(option_list = list(make_option("--candidates", default = "sharedlink,logistic_depth"),
  make_option("--settings", default = "house:R00,implant:B_ref"), make_option("--templates", default = "hmp_tongue,twinsuk_stool"),
  make_option("--reps", type = "integer", default = 3L), make_option("--seed", type = "integer", default = 9001L),
  make_option("--bench-root", default = NULL)))
opt <- parse_args(op)
if (!is.null(opt$`bench-root`)) Sys.setenv(PURSUE_BENCH_ROOT = opt$`bench-root`)
root <- normalizePath(Sys.getenv("PURSUE_BENCH_ROOT", unset = "benchmarks")); Sys.setenv(PURSUE_BENCH_ROOT = root)
for (f in c("R/engine/templates.R", "R/engine/regimes.R", "R/engine/metrics.R", "R/simulators/sim_house.R",
            "R/simulators/implant.R", "R/simulators/dispatch.R", "R/methods/elementary.R")) source(file.path(root, f))
.cands <- new.env(); register_candidate <- function(name, fn, notes = "") assign(name, list(fn = fn), envir = .cands)
for (f in sort(list.files(file.path(root, "dev", "candidates"), pattern = "\\.R$", full.names = TRUE))) source(f)
regimes <- read.delim(file.path(root, "regimes.tsv"), stringsAsFactors = FALSE, comment.char = "#")
spec <- function(signal = "mixed", da = 0.10, balance = "balanced") list(n_per_group = 50, da_frac = da, effect = "medium",
  signal_type = signal, balance = balance, conf_phi = 0, exposure = "binary")
imp <- list(B_ref = spec(), B_null = spec(da = 0), B_prev = spec("prevalence"), B_abund = spec("abundance"), B_allup = spec(balance = "100_0"))
cell_seed <- function(...) { k <- paste(..., sep = "|"); as.integer((opt$seed * 1e5 + sum(utf8ToInt(k) * seq_along(utf8ToInt(k)))) %% .Machine$integer.max) }
cand <- strsplit(opt$candidates, ",")[[1]]
rows <- list()
for (sid in strsplit(opt$settings, ",")[[1]]) for (tid in strsplit(opt$templates, ",")[[1]]) for (rep in seq_len(opt$reps)) {
  tpl <- readRDS(file.path(root, "devdata", paste0(tid, ".rds")))
  sim_kind <- sub(":.*", "", sid); rg <- sub(".*:", "", sid)
  sim <- if (sim_kind == "implant") implant(tpl, imp[[rg]], cell_seed(sid, tid, rep)) else
    simulate_cell_data(sim_kind, tpl, regimes[regimes$regime_id == rg, ], cell_seed(sid, tid, rep))
  keep <- rowMeans(sim$counts > 0) >= 0.10 & rowSums(sim$counts > 0) >= 3; ct <- sim$counts[keep, , drop = FALSE]
  tr <- sim$truth[match(rownames(ct), sim$truth$feature), ]
  meta <- sim$meta[colnames(ct), , drop = FALSE]; if (is.null(meta$depth)) meta$depth <- colSums(sim$counts)
  g <- meta[[sim$tested_term]]; g <- if (is.factor(g)) as.integer(g) - 1L else g
  R <- sweep(sim$counts, 2, colSums(sim$counts), "/")[rownames(ct), , drop = FALSE]
  rel_lfc2 <- log2((rowMeans(R[, g == 1, drop = FALSE]) + 1e-9) / (rowMeans(R[, g == 0, drop = FALSE]) + 1e-9))
  det_shift <- rowMeans(ct[, g == 1, drop = FALSE] > 0) - rowMeans(ct[, g == 0, drop = FALSE] > 0)
  for (cn in cand) {
    r <- .cands[[cn]]$fn(ct, meta, sim$formula, sim$tested_term); p <- r$p[match(rownames(ct), r$feature)]
    q <- rep(NA, length(p)); ok <- is.finite(p); q[ok] <- p.adjust(p[ok], "BH"); rej <- !is.na(q) & q <= 0.05
    rows[[length(rows) + 1L]] <- data.frame(setting = sid, template = tid, rep = rep, candidate = cn, feature = rownames(ct),
      truth_abs = tr$truth_abs, truth_type = tr$truth_type, rel_lfc2 = rel_lfc2, det_shift = det_shift,
      prev = rowMeans(ct > 0), log10_mean_ra = log10(rowMeans(R) + 1e-9), rejected = rej, p = p, stringsAsFactors = FALSE)
  }
  cat(sprintf("  %s %s r%d\n", sid, tid, rep))
}
A <- do.call(rbind, rows)
fp <- A[A$rejected & A$truth_abs == 0, ]; tp <- A[A$rejected & A$truth_abs == 1, ]; nul <- A[A$truth_abs == 0 & is.finite(A$p), ]
cat("\n=== false positives vs other null features ===\n")
for (cn in cand) for (sid in unique(A$setting)) {
  f <- fp[fp$candidate == cn & fp$setting == sid, ]; n0 <- nul[nul$candidate == cn & nul$setting == sid, ]
  if (!nrow(n0)) next
  cat(sprintf("\n%s @ %s: %d false positives over %d cells\n", cn, sid, nrow(f), length(unique(paste(A$template, A$rep)))))
  if (!nrow(f)) next
  cmp <- function(v) sprintf("FP median %6.2f | other nulls median %6.2f", stats::median(f[[v]]), stats::median(n0[[v]]))
  cat("  |relative lfc2|  ", cmp_abs <- sprintf("FP median %6.2f | other nulls median %6.2f", median(abs(f$rel_lfc2)), median(abs(n0$rel_lfc2))), "\n")
  cat("  |detection shift|", sprintf("FP median %6.3f | other nulls median %6.3f", median(abs(f$det_shift)), median(abs(n0$det_shift))), "\n")
  cat("  prevalence       ", cmp("prev"), "\n")
  cat("  log10 rel. abund.", cmp("log10_mean_ra"), "\n")
  cat(sprintf("  FP with |relative lfc2| > 0.5: %d of %d   (other nulls: %.1f%%)\n",
              sum(abs(f$rel_lfc2) > 0.5), nrow(f), 100 * mean(abs(n0$rel_lfc2) > 0.5)))
}
saveRDS(A, file.path(root, "dev", "results", "diag_fp_all.rds"))
out <- file.path(root, "dev", "results", "diag_fp.csv"); utils::write.csv(A[A$rejected | A$truth_abs == 1, ], out, row.names = FALSE)
cat("\nwrote", out, "\n")
