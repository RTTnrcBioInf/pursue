# -----------------------------------------------------------------------------
# metrics.R -- performance measures (protocol section 10) computed from one cell's
# long result table (one row per method x arm x feature, joined with truth).
# -----------------------------------------------------------------------------

#' Truth column that matches an arm and a scale
truth_for <- function(truth, arm, scale = c("abs", "rel")) {
  scale <- match.arg(scale)
  switch(arm,
    presence  = as.integer(truth$truth_type %in% c("prevalence", "both")),
    abundance = as.integer(truth$truth_type %in% c("abundance", "both", "bloom")),
    if (scale == "abs") as.integer(truth$truth_abs) else as.integer(truth$truth_rel))
}

pauc <- function(p, truth, fpr_max = 0.10) {
  ok <- is.finite(p) & !is.na(truth); p <- p[ok]; truth <- truth[ok]
  if (sum(truth == 1) == 0 || sum(truth == 0) == 0) return(NA_real_)
  o <- order(p); t <- truth[o]; tpr <- cumsum(t == 1) / sum(t == 1); fpr <- cumsum(t == 0) / sum(t == 0)
  keep <- fpr <= fpr_max; if (!any(keep)) return(0)
  x <- c(0, fpr[keep]); y <- c(0, tpr[keep])
  sum(diff(x) * (y[-1] + y[-length(y)]) / 2) / fpr_max
}

cell_metrics <- function(res, truth, alpha = c(0.05, 0.10)) {
  # res: feature, arm, p, q, estimate, status ; truth: feature-level truth table
  tr <- truth[match(res$feature, truth$feature), ]
  out <- list()
  for (arm in unique(res$arm)) {
    r <- res[res$arm == arm, ]; t <- tr[res$arm == arm, ]
    for (scale in c("abs", "rel")) {
      y <- truth_for(t, arm, scale); ok <- is.finite(r$p)
      row <- data.frame(arm = arm, truth_scale = scale, n_tested = sum(ok), n_pos = sum(y[ok] == 1), stringsAsFactors = FALSE)
      for (a in alpha) { rej <- ok & !is.na(r$q) & r$q <= a; tp <- sum(rej & y == 1); fp <- sum(rej & y == 0)
        row[[sprintf("fdr_%02d", a * 100)]] <- if (tp + fp > 0) fp / (tp + fp) else 0
        row[[sprintf("tpr_%02d", a * 100)]] <- if (sum(y[ok] == 1) > 0) tp / sum(y[ok] == 1) else NA
        row[[sprintf("n_rej_%02d", a * 100)]] <- tp + fp }
      row$pauc10 <- pauc(r$p, y); pn <- r$p[ok & y == 0]
      row$fpr05 <- if (length(pn)) mean(pn < 0.05) else NA; row$fpr01 <- if (length(pn)) mean(pn < 0.01) else NA
      row$ks_null <- if (length(pn) >= 10) suppressWarnings(stats::ks.test(pn, "punif")$statistic) else NA
      # cross-arm shadow
      if (arm == "presence") { s <- ok & t$truth_type == "abundance"; row$shadow <- if (any(s)) mean(r$p[s] < 0.05) else NA }
      else if (arm == "abundance") { s <- ok & t$truth_type == "prevalence"; row$shadow <- if (any(s)) mean(r$p[s] < 0.05) else NA }
      else row$shadow <- NA
      # effect recovery (abundance-type truths, log2 scale) where estimates exist
      if (arm %in% c("abundance", "single") && "truth_lfc2" %in% names(t) && any(is.finite(r$estimate))) {
        e <- ok & is.finite(r$estimate) & t$truth_type %in% c("abundance", "both", "bloom", "none")
        tv <- if (scale == "abs") t$truth_lfc2 else if ("truth_rel_lfc2" %in% names(t)) t$truth_rel_lfc2 else t$truth_lfc2
        row$est_bias <- mean(r$estimate[e] - tv[e]); row$est_rmse <- sqrt(mean((r$estimate[e] - tv[e])^2))
        ci <- e & is.finite(r$ci_lo) & is.finite(r$ci_hi)
        row$ci_cover <- if (any(ci)) mean(r$ci_lo[ci] <= tv[ci] & r$ci_hi[ci] >= tv[ci]) else NA
      } else { row$est_bias <- NA; row$est_rmse <- NA; row$ci_cover <- NA }
      out[[length(out) + 1L]] <- row
    }
  }
  do.call(rbind, out)
}
