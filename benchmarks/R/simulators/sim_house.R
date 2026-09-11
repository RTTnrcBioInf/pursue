# -----------------------------------------------------------------------------
# sim_house.R -- the in-house generator ("house"): structural presence x lognormal
# abundance x Poisson (or NB) sampling with multinomial closure, parameterised from a
# real template's per-feature prevalence and log relative abundance. It is the
# reference implementation of every regime factor in regimes.R; the external
# simulators map the same factors onto their own parameters as far as they can.
#
# simulate_house(template, regime, seed) -> list(counts [features x samples], meta, truth)
# truth: feature, truth_abs (0/1), truth_type (none/abundance/prevalence/both/bloom),
#        truth_lfc2 (abundance effect, log2), truth_logor (presence effect), truth_rel (0/1)
# -----------------------------------------------------------------------------

simulate_house <- function(template, regime, seed = 1L, sampling = c("poisson", "nb"), nb_size = 5) {
  sampling <- match.arg(sampling); set.seed(seed)
  prof <- template_profile(template)
  m <- regime$m; n_per <- regime$n_per_group
  # ---- samples, exposure, confounder, design ----
  hetero <- identical(regime$hetero, "yes")
  n <- 2L * n_per
  if (identical(regime$design, "repeated")) {
    n_subj <- max(4L, round(n / 5)); visits <- 5L; n <- n_subj * visits
    subject <- rep(seq_len(n_subj), each = visits)
    u_lat <- rep(stats::rnorm(n_subj), each = visits)          # exposure is subject-level
  } else { subject <- seq_len(n); u_lat <- stats::rnorm(n) }
  if (identical(regime$exposure, "continuous")) { x <- as.numeric(scale(u_lat)); x_bin <- NULL }
  else { cut <- if (hetero) stats::qnorm(0.70) else stats::median(u_lat); x <- as.integer(u_lat > cut); x_bin <- x }
  phi <- regime$conf_phi
  conf <- if (phi > 0) phi * as.numeric(scale(u_lat)) + sqrt(1 - phi^2) * stats::rnorm(n) else stats::rnorm(n)
  # ---- depth ----
  N <- sample(prof$depth, n, replace = TRUE)
  if (!is.null(x_bin)) N <- round(N * ifelse(x_bin == 1L, regime$depth_conf, 1)) else N <- round(N * regime$depth_conf^(x > 0))
  N <- pmax(N, 200L)
  # ---- features drawn from the template profile ----
  idx <- sample(seq_len(prof$n_features), m, replace = m > prof$n_features)
  psi <- pmin(pmax(prof$prevalence[idx], 0.05), 0.995)
  mu  <- prof$mean_log_relab[idx]; mu[!is.finite(mu)] <- stats::median(mu, na.rm = TRUE)
  sg  <- prof$sd_log_relab[idx]; sg[!is.finite(sg) | sg < 0.3] <- stats::median(sg[is.finite(sg) & sg >= 0.3], na.rm = TRUE); sg <- pmin(sg, 2.5)
  # ---- effects ----
  es <- effect_sizes(regime$effect)
  n_da <- round(regime$da_frac * m)
  typ <- rep("none", m); lfc2 <- numeric(m); logor <- numeric(m)
  if (n_da > 0) {
    da <- sample.int(m, n_da)
    kind <- switch(regime$signal_type,
      abundance  = rep("abundance", n_da),
      prevalence = rep("prevalence", n_da),
      mixed      = sample(c("abundance", "prevalence", "both"), n_da, replace = TRUE))
    p_up <- switch(regime$balance, balanced = 0.5, "80_20" = 0.8, "100_0" = 1.0)
    sgn <- ifelse(stats::runif(n_da) < p_up, 1, -1)
    mag_a <- if (es$graded) stats::runif(n_da, es$lfc2[1], es$lfc2[2]) else rep(es$lfc2, n_da)
    mag_p <- if (es$graded) stats::runif(n_da, es$logor[1], es$logor[2]) else rep(es$logor, n_da)
    typ[da] <- kind
    lfc2[da]  <- ifelse(kind %in% c("abundance", "both"), sgn * mag_a, 0)
    logor[da] <- ifelse(kind %in% c("prevalence", "both"), sgn * mag_p, 0)
  }
  if (identical(regime$bloom, "x20")) { b <- which(typ == "none")[1]; typ[b] <- "bloom"; lfc2[b] <- log2(20)
    mu[b] <- stats::quantile(mu, 0.95) }        # a bloom of an already-common feature
  # confounder affects 10% of non-DA features
  gconf <- numeric(m)
  if (phi > 0) { cf <- sample(which(typ == "none"), round(0.10 * m)); gconf[cf] <- log(2) * sample(c(-1, 1), length(cf), TRUE) }
  # ---- generate ----
  beta <- lfc2 * log(2)
  icc <- if (identical(regime$design, "repeated")) 0.3 else 0
  A <- matrix(0, n, m); S <- matrix(0L, n, m)
  for (j in seq_len(m)) {
    eta_psi <- stats::qlogis(psi[j]) + logor[j] * x
    S[, j] <- stats::rbinom(n, 1L, stats::plogis(eta_psi))
    u_subj <- if (icc > 0) rep(stats::rnorm(max(subject), 0, sqrt(icc / (1 - icc)) * sg[j]), each = 5L)[seq_len(n)] else 0
    eps <- stats::rnorm(n, 0, sg[j]); if (hetero && !is.null(x_bin)) eps <- eps * ifelse(x_bin == 1L, sqrt(3), 1)
    A[, j] <- exp(mu[j] + beta[j] * x + gconf[j] * conf + u_subj + eps) * S[, j]
  }
  tot <- rowSums(A); tot[tot <= 0] <- 1
  R <- A / tot
  X <- matrix(0L, n, m)
  for (i in seq_len(n)) if (sum(R[i, ]) > 0) {
    if (sampling == "poisson") X[i, ] <- stats::rpois(m, N[i] * R[i, ])
    else X[i, ] <- stats::rnbinom(m, size = nb_size, mu = N[i] * R[i, ])
  }
  feat <- sprintf("F%04d", seq_len(m)); colnames(X) <- feat; rownames(X) <- sprintf("S%04d", seq_len(n))
  # ---- truth on the relative scale: realised mean relative abundance by exposure ----
  if (!is.null(x_bin)) {
    r1 <- colMeans(R[x_bin == 1L, , drop = FALSE]); r0 <- colMeans(R[x_bin == 0L, , drop = FALSE])
    rel_lfc2 <- log2((r1 + 1e-12) / (r0 + 1e-12))
  } else { rel_lfc2 <- apply(R, 2L, function(r) stats::coef(stats::lm(log(r + 1e-12) ~ x))[2]) / log(2) }
  truth <- data.frame(feature = feat, truth_type = typ,
                      truth_abs = as.integer(typ != "none"),
                      truth_lfc2 = lfc2, truth_logor = logor, truth_rel_lfc2 = rel_lfc2,
                      truth_rel = as.integer(abs(rel_lfc2) > 0.25 & (typ != "none" | abs(rel_lfc2) > 0.25)),
                      conf_affected = as.integer(gconf != 0), stringsAsFactors = FALSE)
  meta <- data.frame(row.names = rownames(X),
                     group = if (!is.null(x_bin)) factor(ifelse(x_bin == 1L, "case", "control"), levels = c("control", "case")) else NA,
                     exposure = x, confounder = conf, subject = factor(subject), depth = N, stringsAsFactors = FALSE)
  list(counts = t(X), meta = meta, truth = truth,
       tested_term = if (!is.null(x_bin)) "group" else "exposure",
       formula = stats::as.formula(paste("~", if (!is.null(x_bin)) "group" else "exposure", if (phi > 0) "+ confounder" else "")))
}
