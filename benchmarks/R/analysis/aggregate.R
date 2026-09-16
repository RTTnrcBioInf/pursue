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
# --- cell-level design audit (protocol section 4) -----------------------------------------
# Every simulator clamps the requested feature count to the template's own: `min(regime$m,
# nrow(counts))`. So R04 / R00 / R05 -- the m = 200 / 500 / 1000 levels of the feature-count
# factor -- can deliver the SAME m. On hmp_gingiva (364 features) they are 200 / 364 / 364, and
# an m effect fitted across all templates would be reading a difference that is not there.
# Cells written before 2026-09-14 carry no m_clamped field in their manifest, so DERIVE it from
# n_features_total and regimes.tsv rather than trusting the field; that covers every cell ever
# written and needs nothing re-run.
script <- sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])
bench <- tryCatch(normalizePath(file.path(dirname(script), "..", "..")), error = function(e) getwd())
mf <- list.files(res_dir, pattern = "\\.manifest\\.json$", full.names = TRUE, recursive = TRUE)
cells <- NULL
if (length(mf) && requireNamespace("jsonlite", quietly = TRUE)) {
  one <- function(f) {
    j <- tryCatch(jsonlite::fromJSON(f, simplifyVector = TRUE), error = function(e) NULL)
    if (is.null(j) || is.null(j$cell)) return(NULL)
    g <- function(k) { v <- j[[k]]; if (is.null(v) || !length(v) || is.na(v[[1]])) NA else v[[1]] }
    data.frame(axis = j$cell$axis, simulator = j$cell$simulator, template = j$cell$template,
               regime_id = j$cell$regime_id, replicate = as.integer(j$cell$replicate),
               n_samples = as.integer(g("n_samples")), n_features_total = as.integer(g("n_features_total")),
               n_features_tested = as.integer(g("n_features_tested")), stringsAsFactors = FALSE)
  }
  cells <- do.call(rbind, lapply(mf, one))
}
if (!is.null(cells) && nrow(cells)) {
  rp <- file.path(bench, "regimes.tsv")
  if (file.exists(rp)) {
    reg <- read.delim(rp, comment.char = "#", stringsAsFactors = FALSE)
    cells$m_requested <- ifelse(cells$axis == "A", reg$m[match(cells$regime_id, reg$regime_id)], NA)
    cells$m_clamped <- !is.na(cells$m_requested) & cells$n_features_total < cells$m_requested
  }
  write.csv(cells, file.path(out, "cells.csv"), row.names = FALSE)
  sub <- cells[cells$axis == "A" & cells$regime_id %in% c("R04", "R00", "R05"), ]
  if (nrow(sub) && "m_clamped" %in% names(cells)) {
    # BY SIMULATOR, not averaged over them: each simulator clamps differently (sd2 ignores m
    # entirely and uses the template's features), so a mean across simulators shows numbers no
    # cell ever had. The 2026-09-14 run printed hmp_gingiva R05 = 523 that way, on a template
    # with 364 features.
    sub$sim_tpl <- paste(sub$simulator, sub$template, sep = " / ")
    mt <- tapply(sub$n_features_total, list(sub$sim_tpl, sub$regime_id), function(x) round(mean(x, na.rm = TRUE)))
    mt <- mt[, intersect(c("R04", "R00", "R05"), colnames(mt)), drop = FALSE]
    cat("\nrealised feature count per simulator/template, m factor (requested 200/500/1000):\n"); print(mt)
    collapsed <- rownames(mt)[apply(mt, 1, function(r) { r <- r[!is.na(r)]; length(r) > 1 && length(unique(r)) < length(r) })]
    if (length(collapsed))
      cat("  !! m is NOT estimable on: ", paste(collapsed, collapse = ", "),
          "\n     (levels collapse to the same feature count). Drop these from any m contrast.\n", sep = "")
    cat(sprintf("  %d of %d axis-A cells ran at a clamped m.\n",
                sum(cells$m_clamped[cells$axis == "A"], na.rm = TRUE), sum(cells$axis == "A")))
  }
  key <- c("axis", "simulator", "template", "regime_id", "replicate")
  M <- merge(M, cells, by = key, all.x = TRUE, sort = FALSE)
}

# --- design coverage: which cells are structurally absent, and why -------------------------
# Not every simulator can produce every regime (msq cannot model depth confounding; only the
# in-house simulator does bloom/hetero/repeated designs; msq cannot draw more samples than the
# template has). Those cells write <id>.skipped.json instead of a manifest. A reader needs to
# see them: an empty cell in the Axis A table is a design fact, not a missing result.
sk <- list.files(res_dir, pattern = "\\.skipped\\.json$", full.names = TRUE, recursive = TRUE)
if (length(sk) && requireNamespace("jsonlite", quietly = TRUE)) {
  ones <- function(f) {
    j <- tryCatch(jsonlite::fromJSON(f, simplifyVector = TRUE), error = function(e) NULL)
    if (is.null(j) || is.null(j$cell)) return(NULL)
    data.frame(axis = j$cell$axis, simulator = j$cell$simulator, template = j$cell$template,
               regime_id = j$cell$regime_id, replicate = as.integer(j$cell$replicate),
               reason = paste(unlist(j$unsupported), collapse = "+"), stringsAsFactors = FALSE)
  }
  S <- do.call(rbind, lapply(sk, ones))
  if (!is.null(S) && nrow(S)) {
    write.csv(S, file.path(out, "coverage_gaps.csv"), row.names = FALSE)
    cat(sprintf("\n%d cells not producible, by reason:\n", nrow(S)))
    print(sort(table(S$reason), decreasing = TRUE))
    cat("\nsimulators x regimes NOT producible (axis A, cell = replicates missing):\n")
    A <- S[S$axis == "A", ]
    if (nrow(A)) print(table(A$simulator, A$regime_id))
  }
}

write.csv(M, file.path(out, "metrics_long.csv"), row.names = FALSE)
cat(sprintf("%d cells, %d rows\n", length(fs), nrow(M)))

mse <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) NA else sd(x) / sqrt(length(x)) }
summ <- function(d, by) {
  vals <- c("fdr_05", "tpr_05", "fdr_10", "tpr_10", "pauc10", "fpr05", "ks_null", "shadow", "est_bias", "est_rmse", "ci_cover", "runtime_s",
            "n_pos", "n_tested", "n_rej_05", "tp_05", "fp_05")
  vals <- intersect(vals, names(d))
  a <- aggregate(d[vals], d[by], function(x) mean(x, na.rm = TRUE)); s <- aggregate(d[vals], d[by], mse)
  names(s)[names(s) %in% vals] <- paste0(vals, "_se"); n <- aggregate(list(n_cells = d[[vals[1]]]), d[by], length)
  merge(merge(a, s, by = by), n, by = by)
}
# Per-cell counts, averaged by summ() below. TPR x n_pos is the number of true positives the
# arm actually found; n_rej x FDR is the number of false ones. Rates mislead across arms.
if (all(c("tpr_05", "n_pos") %in% names(M))) M$tp_05 <- M$tpr_05 * M$n_pos
if (all(c("n_rej_05", "fdr_05") %in% names(M))) M$fp_05 <- M$n_rej_05 * M$fdr_05

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
