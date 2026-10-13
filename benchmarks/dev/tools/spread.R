# -----------------------------------------------------------------------------
# spread.R -- spread-only implants into real samples (2026-10-13, after the real-data transfer check).
#
# R25 (30 / 70 groups, 3x log-scale variance in the cases) is the only regime where groups differ in spread, and only
# `house` generates it; every method fails it (c10 FDR 0.19, the best). Real case groups are often more variable than
# controls. This implants a spread change into real tuning-template samples, with no simulator, and leaves the
# location of every taxon unchanged, so that under the benchmark's location truth every call is a false one.
#
# For an altered taxon j, on the samples where it is present: l = log(count / depth) and mu = the median of l over all
# drawn samples; l' = mu + k_i (l - mu), with k_i <= 1 the sample's shrink factor. The spread difference is made by
# shrinking, never by stretching: stretching a taxon's log-scale spread by sqrt(3) turns its upper tail into blooms
# (a taxon at 30% where its median is 5% goes to ~110% before rescaling; twinsuk_stool's long-tailed taxa swamp whole
# samples), which is a bloom test, not a spread test. So, binary: the controls are shrunk by 1 / k and the cases are
# real samples, untouched -- the cases have k^2 times the controls' log-scale variance (R25's 3x with k = sqrt(3), in
# the small group). Continuous: k_i = k^((z_i - max z) / 2) for standardised z (variance ratio k^2 between z = +1 and
# z = -1, as before; the highest-z sample untouched). Zeros stay zeros. Each sample is then rescaled to its original
# depth (a per-sample factor: the compositional side effect every centred method must absorb) and either rounded with
# every present cell kept >= 1 (redraw = FALSE: presence and the median of the present values unchanged) or
# multinomially redrawn as implant() does (redraw = TRUE: a low value can fall to zero, so presence drops a little in
# the high-spread group, as it would in real data).
#
# Present log abundances are right-skewed (a pile of low counts, a long upper tail), so shrinking around the median
# also moves the mean of log abundance, by (1 - 1 / k)(mean - median). Measured over 12 draws per template (case -
# control, sd units, altered taxa): H_unb mean +0.03 / -0.02, median +0.00 / -0.02 (twinsuk / tongue; single draws
# vary by +-0.07 because the taxa share samples). centre = "mean" shrinks around the mean instead (median then -0.07 /
# -0.05). Both are small against the spread change, so the default stays the median and no mean-centred setting is run.
#
# implant_spread(template, spec, seed); spec = list(n_case, n_ctrl, frac, k, exposure = "binary"|"continuous",
#   redraw, centre = "median"|"mean"). Altered taxa are drawn among those present in >= 10% of the drawn samples (the dev suite's filter), so
#   `frac` is the altered share of tested taxa. truth_abs is 0 everywhere; `spread` marks the altered taxa.
# -----------------------------------------------------------------------------

implant_spread <- function(template, spec, seed = 1L) {
  set.seed(seed)
  ct <- template$counts; m <- nrow(ct); n_avail <- ncol(ct)
  continuous <- identical(spec$exposure, "continuous")
  n <- spec$n_case + spec$n_ctrl
  if (n > n_avail) return(list(unsupported = sprintf("n_exceeds_template[%d>%d]", n, n_avail)))
  s <- sample.int(n_avail, n); X <- ct[, s, drop = FALSE]; N <- colSums(X)
  if (continuous) { x <- as.numeric(scale(stats::rnorm(n))); kf <- spec$k^((x - max(x)) / 2)
  } else { x <- sample(c(rep(1L, spec$n_case), rep(0L, spec$n_ctrl))); kf <- ifelse(x == 1L, 1, 1 / spec$k) }
  elig <- which(rowMeans(X > 0) >= 0.10 & rowSums(X > 0) >= 3)
  alt <- if (spec$frac >= 1) elig else sort(sample(elig, round(spec$frac * length(elig))))
  Xf <- X * 1.0
  for (j in alt) {
    pos <- X[j, ] > 0; l <- log(X[j, pos] / N[pos]); mu <- if (identical(spec$centre, "mean")) mean(l) else stats::median(l)
    Xf[j, pos] <- N[pos] * exp(mu + kf[pos] * (l - mu))
  }
  # location diagnostic on the absolute (pre-rescaling) scale, altered taxa, present cells: the case - control
  # difference (binary) or the correlation with z (continuous) of the mean and the median of log abundance, in units
  # of the taxon's sd. The median is preserved by construction up to the random split; the mean moves by
  # (k - 1)(mean - median) for a skewed taxon.
  diag <- sapply(alt, function(j) { pos <- Xf[j, ] > 0; l <- log(Xf[j, pos] / N[pos]); sdl <- stats::sd(l); g <- x[pos]
    if (sdl <= 0 || !is.finite(sdl)) return(c(NA, NA))
    if (continuous) return(c(stats::cor(l, g), stats::cor(l, g, method = "spearman")))
    if (sum(g == 1L) < 2 || sum(g == 0L) < 2) return(c(NA, NA))
    c((mean(l[g == 1L]) - mean(l[g == 0L])) / sdl, (stats::median(l[g == 1L]) - stats::median(l[g == 0L])) / sdl) })
  Xs <- sweep(Xf, 2L, N / colSums(Xf), "*")
  if (isTRUE(spec$redraw)) {
    Xn <- sapply(seq_len(n), function(i) stats::rmultinom(1L, N[i], Xs[, i] / sum(Xs[, i]))[, 1])
  } else { Xn <- round(Xs); Xn[X > 0 & Xn < 1] <- 1 }
  dimnames(Xn) <- list(rownames(ct), paste0("S", seq_len(n))); storage.mode(Xn) <- "integer"
  sp <- integer(m); sp[alt] <- 1L
  truth <- data.frame(feature = rownames(ct), truth_type = "none", truth_abs = 0L, truth_lfc2 = 0, truth_rel = 0L,
                      spread = sp, stringsAsFactors = FALSE)
  meta <- data.frame(row.names = colnames(Xn),
                     group = if (continuous) NA else factor(ifelse(x == 1L, "case", "control"), levels = c("control", "case")),
                     exposure = x, depth = colSums(Xn), stringsAsFactors = FALSE)
  list(counts = Xn, meta = meta, truth = truth, tested_term = if (continuous) "exposure" else "group",
       formula = stats::as.formula(paste("~", if (continuous) "exposure" else "group")),
       loc_diag = c(mean_shift = mean(diag[1, ], na.rm = TRUE), median_shift = mean(diag[2, ], na.rm = TRUE),
                    abs_mean_shift = mean(abs(diag[1, ]), na.rm = TRUE)))
}

# The settings. H_unb is R25 on real samples restricted to 20% of taxa; H_all alters every eligible taxon, as R25 does;
# H_bal separates the imbalance from the spread; H_rd lets presence follow the spread; H_cont is the continuous analogue
# (spread trending with the exposure), which it54's plain-correlation statistic was never checked against.
SPREAD_SPECS <- list(
  H_unb  = list(n_case = 30, n_ctrl = 70, frac = 0.2, k = sqrt(3), exposure = "binary", redraw = FALSE),
  H_bal  = list(n_case = 50, n_ctrl = 50, frac = 0.2, k = sqrt(3), exposure = "binary", redraw = FALSE),
  H_all  = list(n_case = 30, n_ctrl = 70, frac = 1.0, k = sqrt(3), exposure = "binary", redraw = FALSE),
  H_rd   = list(n_case = 30, n_ctrl = 70, frac = 0.2, k = sqrt(3), exposure = "binary", redraw = TRUE),
  H_cont = list(n_case = 50, n_ctrl = 50, frac = 0.2, k = sqrt(3), exposure = "continuous", redraw = FALSE))
