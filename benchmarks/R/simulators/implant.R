# -----------------------------------------------------------------------------
# implant.R -- Axis B: signal implantation into a real template (SIMBA semantics,
# Wirbel et al. 2024). Two random groups of n are drawn from the template; effects are
# implanted into a fraction of features by count scaling (abundance), by moving
# non-zero cells across groups (prevalence), or both; an optional confounder signal is
# implanted into a disjoint feature set correlated with the exposure at phi; depth is
# re-drawn by multinomial resampling so counts stay integers.
#
# implant(template, spec, seed) with spec = list(n_per_group, da_frac, effect, signal_type,
#   balance, conf_phi, exposure = "binary"|"continuous")
# -----------------------------------------------------------------------------

implant <- function(template, spec, seed = 1L) {
  set.seed(seed)
  ct <- template$counts; m <- nrow(ct); n_avail <- ncol(ct)
  n <- 2L * spec$n_per_group
  if (n > n_avail) stop("template has ", n_avail, " samples; need ", n)
  s <- sample.int(n_avail, n); X <- ct[, s, drop = FALSE]
  u <- stats::rnorm(n)
  continuous <- identical(spec$exposure, "continuous")
  x <- if (continuous) as.numeric(scale(u)) else as.integer(u > stats::median(u))
  phi <- if (is.null(spec$conf_phi)) 0 else spec$conf_phi
  conf <- if (phi > 0) phi * as.numeric(scale(u)) + sqrt(1 - phi^2) * stats::rnorm(n) else stats::rnorm(n)
  es <- effect_sizes(spec$effect)
  n_da <- round(spec$da_frac * m)
  typ <- rep("none", m); lfc2 <- numeric(m); shift <- numeric(m)
  Xf <- X * 1.0
  if (n_da > 0) {
    da <- sample.int(m, n_da)
    kind <- switch(spec$signal_type, abundance = rep("abundance", n_da), prevalence = rep("prevalence", n_da),
                   mixed = sample(c("abundance", "prevalence", "both"), n_da, TRUE))
    p_up <- switch(spec$balance, balanced = 0.5, "80_20" = 0.8, "100_0" = 1.0)
    sgn <- ifelse(stats::runif(n_da) < p_up, 1, -1)
    mag_a <- if (es$graded) stats::runif(n_da, es$lfc2[1], es$lfc2[2]) else rep(es$lfc2, n_da)
    mag_p <- if (es$graded) stats::runif(n_da, 0.2, 0.5) else rep(switch(spec$effect, small = 0.2, medium = 0.35, large = 0.5), n_da)
    typ[da] <- kind
    for (k in seq_len(n_da)) {
      j <- da[k]
      if (kind[k] %in% c("abundance", "both")) {
        lfc2[j] <- sgn[k] * mag_a[k]
        fac <- if (continuous) 2^(lfc2[j] * x) else ifelse(x == 1L, 2^lfc2[j], 1)
        Xf[j, ] <- Xf[j, ] * fac
      }
      if (kind[k] %in% c("prevalence", "both") && !continuous) {
        # move a fraction of non-zero cells from the "down" group to zero, and seed the same
        # number of zero cells in the "up" group with a typical non-zero value
        shift[j] <- sgn[k] * mag_p[k]
        up <- if (sgn[k] > 0) x == 1L else x == 0L
        nz_down <- which(!up & Xf[j, ] > 0); z_up <- which(up & Xf[j, ] == 0)
        k_move <- min(length(nz_down), length(z_up), round(mag_p[k] * length(nz_down)))
        if (k_move > 0) {
          med_nz <- stats::median(Xf[j, Xf[j, ] > 0])
          Xf[j, sample(nz_down, k_move)] <- 0
          Xf[j, sample(z_up, k_move)] <- med_nz
        }
      }
    }
  }
  gconf <- numeric(m)
  if (phi > 0) { cf <- sample(which(typ == "none"), round(0.10 * m)); gconf[cf] <- sample(c(-1, 1), length(cf), TRUE)
    for (j in cf) Xf[j, ] <- Xf[j, ] * 2^(gconf[j] * conf) }
  # re-draw depth by multinomial resampling: keeps integers, depth distribution of the template
  N <- colSums(X)
  Xn <- sapply(seq_len(n), function(i) { p <- Xf[, i]; if (sum(p) <= 0) rep(0L, m) else stats::rmultinom(1L, N[i], p / sum(p))[, 1] })
  dimnames(Xn) <- list(rownames(ct), paste0("S", seq_len(n)))
  storage.mode(Xn) <- "integer"
  truth <- data.frame(feature = rownames(ct), truth_type = typ, truth_abs = as.integer(typ != "none"),
                      truth_lfc2 = lfc2, truth_prev_shift = shift, truth_rel = as.integer(typ != "none"),
                      conf_affected = as.integer(gconf != 0), stringsAsFactors = FALSE)
  meta <- data.frame(row.names = colnames(Xn),
                     group = if (continuous) NA else factor(ifelse(x == 1L, "case", "control"), levels = c("control", "case")),
                     exposure = x, confounder = conf, depth = N, stringsAsFactors = FALSE)
  list(counts = Xn, meta = meta, truth = truth,
       tested_term = if (continuous) "exposure" else "group",
       formula = stats::as.formula(paste("~", if (continuous) "exposure" else "group", if (phi > 0) "+ confounder" else "")))
}
