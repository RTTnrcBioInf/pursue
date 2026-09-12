#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# aggregate.R -- collect every *.metrics.csv under a results directory into
#   results/summary/metrics_long.csv      one row per cell x method x arm x truth scale
#   results/summary/table_axisA.csv       per (simulator, regime, method, arm): mean +/- MC se
#   results/summary/table_axisB.csv, table_axisC.csv
#   results/summary/interaction.txt       method x simulator interaction (Axis A, lme4 if present)
#
#   Rscript R/analysis/aggregate.R results/
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE); res_dir <- if (length(args)) args[1] else "results"
out <- file.path(res_dir, "summary"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
fs <- list.files(res_dir, pattern = "\\.metrics\\.csv$", full.names = TRUE, recursive = TRUE)
if (!length(fs)) stop("no metrics files under ", res_dir)
M <- do.call(rbind, lapply(fs, read.csv, stringsAsFactors = FALSE))
# A method re-run at a new version (PURSUE iterations) must not be averaged together with the
# old one. Disambiguate only where more than one version is present, so ordinary runs are
# unaffected and the column set stays the same.
nv <- tapply(M$method_version, M$method, function(v) length(unique(v[!is.na(v)])))
multi <- names(nv)[!is.na(nv) & nv > 1]
if (length(multi)) {
  cat("multiple versions present for:", paste(multi, collapse = ", "), "-- labelling as method@version\n")
  i <- M$method %in% multi & !is.na(M$method_version)
  M$method[i] <- paste0(M$method[i], "@", M$method_version[i])
}
write.csv(M, file.path(out, "metrics_long.csv"), row.names = FALSE)
cat(sprintf("%d cells, %d rows\n", length(fs), nrow(M)))

mse <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) NA else sd(x) / sqrt(length(x)) }
summ <- function(d, by) {
  vals <- c("fdr_05", "tpr_05", "fdr_10", "tpr_10", "pauc10", "fpr05", "ks_null", "shadow", "est_bias", "est_rmse", "ci_cover", "runtime_s")
  vals <- intersect(vals, names(d))
  a <- aggregate(d[vals], d[by], function(x) mean(x, na.rm = TRUE)); s <- aggregate(d[vals], d[by], mse)
  names(s)[names(s) %in% vals] <- paste0(vals, "_se"); n <- aggregate(list(n_cells = d[[vals[1]]]), d[by], length)
  merge(merge(a, s, by = by), n, by = by)
}
for (ax in unique(M$axis)) {
  d <- M[M$axis == ax & M$truth_scale == "abs", ]
  by <- c("simulator", "regime_id", "method", "arm"); if (ax == "C") by <- c("regime_id", "method", "arm")
  t <- summ(d, by); write.csv(t, file.path(out, paste0("table_axis", ax, ".csv")), row.names = FALSE)
  cat("axis", ax, ":", nrow(t), "summary rows\n")
}

# --- method x simulator interaction on Axis A (protocol 4.5) ---
A <- M[M$axis == "A" & M$truth_scale == "abs" & M$arm %in% c("single", "combined"), ]
if (nrow(A) && length(unique(A$simulator)) >= 2 && requireNamespace("lme4", quietly = TRUE)) {
  A$y_fdr <- qlogis(pmin(pmax(A$fdr_05, 0.005), 0.995)); A$y_tpr <- qlogis(pmin(pmax(A$tpr_05, 0.005), 0.995))
  sink(file.path(out, "interaction.txt"))
  for (y in c("y_fdr", "y_tpr")) {
    cat("\n=====", y, "=====\n")
    f1 <- lme4::lmer(as.formula(paste(y, "~ method * simulator + method * regime_id + (1 | template)")), data = A, REML = FALSE)
    f0 <- lme4::lmer(as.formula(paste(y, "~ method + simulator + method * regime_id + (1 | template)")), data = A, REML = FALSE)
    print(anova(f0, f1))
    # per-method fragility: SD across simulators of the method's simulator-specific mean
    em <- aggregate(A[[y]], A[c("method", "simulator")], mean, na.rm = TRUE)
    fr <- aggregate(em$x, list(method = em$method), sd); names(fr)[2] <- "fragility_sd_across_simulators"
    print(fr[order(fr$fragility_sd_across_simulators), ], row.names = FALSE)
  }
  sink()
  cat("interaction analysis -> summary/interaction.txt\n")
} else cat("interaction analysis skipped (need >= 2 simulators on Axis A and lme4)\n")
