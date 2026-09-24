# -----------------------------------------------------------------------------
# score.R -- the charter's scoring, shared by devsuite.R and compare.R so they cannot drift.
#   summarise_cells(D)  per setting x candidate, pooled over templates and replicates
#   scoreboard(sm)      the gate (calibration), then power
#   print_report(sm, sb)
# Gate, in three parts:
#  * FPR among truly null features at p < 0.05 <= 0.075 in EVERY setting. Thousands of null
#    features per setting make this precise (Monte Carlo SE ~0.004), so 1.5x nominal is a real
#    signal of miscalibration.
#  * FDR is the mean over cells of the false discovery proportion (0 where nothing is rejected),
#    as BH defines it. Per setting it is noisy, and false discoveries cluster within cells, so a
#    pooled binomial test on FP out of rejections is wrong: it flagged 2 FP in 2 rejections over 10
#    cells, which is a FWER of 2/10 and entirely compatible with 0.05 (it1, logistic_depth, R19).
#    FDP is bounded in [0,1], and for a fixed mean the Bernoulli is the most dispersed such variable,
#    so the test treats each cell's FDP as at worst a Bernoulli(0.05) draw: p = P(Binomial(cells, 0.05)
#    >= sum of FDPs, rounded up). Clustering within a cell cannot inflate it.
#    A setting fails if mean FDP > 0.10 AND p < 0.01.
#  * Pooled over every setting with signal: mean FDP > 0.075 AND p < 0.01, same test. This catches
#    systematic mild excess that no single setting can resolve (it1: sharedlink ~0.14 in 7 of 10).
# Power: TP per cell relative to the best candidate in each setting, averaged over settings, so
# easy settings do not dominate. Relative power depends on which candidates are compared.
# -----------------------------------------------------------------------------
GATE_FPR <- 0.075; GATE_FDR <- 0.10; GATE_FDR_POOLED <- 0.075
.fdp_p <- function(fdp) if (length(fdp)) stats::pbinom(ceiling(sum(fdp) - 1e-9) - 1, length(fdp), 0.05, lower.tail = FALSE) else 1

summarise_cells <- function(D) {
  structure(do.call(rbind, lapply(split(D, list(D$setting, D$candidate), drop = TRUE), function(x) data.frame(
    setting = x$setting[1], candidate = x$candidate[1], cells = nrow(x),
    n_pos = mean(x$n_pos), TP = mean(x$tp), FP = mean(x$fp),
    FDR = mean(x$fp / pmax(x$rej, 1)), FDR_pool = sum(x$fp) / max(sum(x$rej), 1),
    FPR = sum(x$null_p05) / max(sum(x$null_n), 1), FDR_p = .fdp_p(x$fp / pmax(x$rej, 1)),
    secs = mean(x$secs), errors = sum(!is.na(x$error)), stringsAsFactors = FALSE))), cells = D)
}

scoreboard <- function(sm) {
  sig <- sm$n_pos > 0
  best <- vapply(split(sm$TP[sig], sm$setting[sig]), max, numeric(1))
  sm$rel <- ifelse(sig, sm$TP / pmax(unname(best[sm$setting]), 1e-9), NA)
  sb <- do.call(rbind, lapply(split(sm, sm$candidate), function(x) {
    s <- x$n_pos > 0
    data.frame(candidate = x$candidate[1],
      worst_FPR = max(x$FPR, na.rm = TRUE), worst_FPR_at = x$setting[which.max(x$FPR)],
      worst_FDR = if (any(s)) max(x$FDR[s]) else NA, worst_FDR_at = if (any(s)) x$setting[s][which.max(x$FDR[s])] else NA,
      mean_TP = mean(x$TP[s]), rel_power = mean(x$rel, na.rm = TRUE),
      secs_per_cell = mean(x$secs), errors = sum(x$errors), settings = nrow(x), stringsAsFactors = FALSE) }))
  bad_fdr <- tapply(sm$n_pos > 0 & sm$FDR > GATE_FDR & sm$FDR_p < 0.01, sm$candidate, any)
  sb$fdr_flags <- vapply(sb$candidate, function(c) { x <- sm[sm$candidate == c & sm$n_pos > 0 & sm$FDR > GATE_FDR & sm$FDR_p < 0.01, ]
    if (nrow(x)) paste(x$setting, collapse = ",") else "" }, character(1))
  # pooled FDP over all signal cells needs the cells themselves: summarise_cells keeps them as an attribute
  cl <- attr(sm, "cells")
  pooled <- vapply(sb$candidate, function(c) { x <- cl[cl$candidate == c & cl$n_pos > 0, ]; fdp <- x$fp / pmax(x$rej, 1)
    c(mean(fdp), .fdp_p(fdp)) }, numeric(2))
  sb$FDR_all <- pooled[1, ]; sb$FDR_all_p <- pooled[2, ]
  bad_pool <- sb$FDR_all > GATE_FDR_POOLED & sb$FDR_all_p < 0.01
  sb$fdr_flags <- ifelse(bad_pool, ifelse(nzchar(sb$fdr_flags), paste0("ALL,", sb$fdr_flags), "ALL"), sb$fdr_flags)
  sb$calibrated <- sb$worst_FPR <= GATE_FPR & !bad_fdr[sb$candidate] & !bad_pool & sb$errors == 0
  sb[order(!sb$calibrated, -sb$rel_power), ]
}

print_report <- function(sm, sb, per_setting = TRUE) {
  fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
  cat("=== scoreboard: calibration gate, then power ===\n")
  print(data.frame(candidate = sb$candidate, gate = ifelse(sb$calibrated, "PASS", "fail"),
                   worst_FPR = fmt(sb$worst_FPR), at = sb$worst_FPR_at, FDR_all = fmt(sb$FDR_all), sig_FDR_excess = ifelse(nzchar(sb$fdr_flags), sb$fdr_flags, "-"),
                   rel_power = fmt(sb$rel_power, 2), mean_TP = fmt(sb$mean_TP, 2), s_cell = fmt(sb$secs_per_cell, 1)),
        row.names = FALSE)
  if (!per_setting) return(invisible())
  cat("\n=== per setting: TP / FDR / FPR ===\n")
  for (s in sort(unique(sm$setting))) {
    x <- sm[sm$setting == s, ]; x <- x[order(-x$TP), ]
    cat(sprintf("\n%s  (DA per cell %.1f)\n", s, max(x$n_pos)))
    print(data.frame(candidate = x$candidate, TP = fmt(x$TP, 2), FP = fmt(x$FP, 2), FDR = fmt(x$FDR), FPR = fmt(x$FPR),
                     s = fmt(x$secs, 1)), row.names = FALSE)
  }
}
