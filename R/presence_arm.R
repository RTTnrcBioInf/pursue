# -----------------------------------------------------------------------------
# presence_arm.R -- structural presence via a stabilized zero-inflated NB
#
# Per feature j:
#   logit psi_i = Xz_i' alpha                 (structural presence; Xz = design incl. tested)
#   log  mu_i   = log N_i + Xc_i' gamma       (abundance when present; NB with dispersion theta)
#   P(X=0) = (1 - psi) + psi * NB(0 | mu, theta);   P(X = x > 0) = psi * NB(x | mu, theta)
#
# Stabilizers (see design-experiments/step4_joint):
#   * theta_j from the zero-truncated NB on detected cells, shrunk on the log scale toward a
#     loess trend across features (prior weight d0_prior pseudo-observations), then fixed.
#   * Gaussian ridge on the presence coefficients (sd tau_ridge) against separation.
#   * Likelihood-ratio test of the tested term's presence coefficients (df = #columns),
#     constrained model refitted from the full solution, several starts.
# The count part exists to model detection; fold changes come from the abundance arm.
# -----------------------------------------------------------------------------

#' Fit the presence arm for every feature
#' @keywords internal
fit_presence_arm <- function(otu, depth, design, d0_prior = 10, tau_ridge = 5,
                             min_detected = 3L, n_cores = 1L) {
  m <- ncol(otu); n <- nrow(otu); lN <- log(depth)
  Xz <- design$X; Xc <- design$X
  tcols <- design$tested_cols; k <- length(tcols); cn <- design$coef_names[tcols]
  pz <- ncol(Xz); pc <- ncol(Xc)

  # ---- stage 1: per-feature dispersion from the zero-truncated NB on detected cells ----
  ztnb_theta <- function(j) {
    xj <- otu[, j]; det <- xj > 0; n1 <- sum(det)
    if (n1 < 3L) return(c(NA_real_, NA_real_))
    xd <- xj[det]; ld <- lN[det]
    Xd <- if (n1 >= pc + 3L && qr(Xc[det, , drop = FALSE])$rank == pc) Xc[det, , drop = FALSE] else matrix(1, n1, 1L)
    q <- ncol(Xd)
    nll <- function(par) { mu <- exp(ld + as.vector(Xd %*% par[seq_len(q)])); th <- exp(par[q + 1L])
      lp0 <- th * (log(th) - log(th + mu))
      -sum(stats::dnbinom(xd, size = th, mu = mu, log = TRUE) - log1p(-pmin(exp(lp0), 1 - 1e-12))) }
    st <- c(mean(log(xd / depth[det])), rep(0, q - 1L), 0)
    f <- tryCatch(stats::nlminb(st, nll, lower = c(-30, rep(-10, q - 1L), -5), upper = c(10, rep(10, q - 1L), 8)),
                  error = function(e) NULL)
    if (is.null(f) || !is.finite(f$objective)) return(c(NA_real_, mean(log(xd / depth[det]))))
    c(f$par[q + 1L], mean(log(xd / depth[det])))
  }
  th <- if (n_cores > 1L) parallel::mclapply(seq_len(m), ztnb_theta, mc.cores = n_cores) else lapply(seq_len(m), ztnb_theta)
  th <- do.call(rbind, th); log_theta <- th[, 1L]; amean <- th[, 2L]
  n_det <- colSums(otu > 0)
  ok_th <- is.finite(log_theta) & is.finite(amean)
  trend <- rep(if (any(ok_th)) stats::median(log_theta[ok_th]) else 0, m)
  if (sum(ok_th) >= 20L) {
    lo <- tryCatch(suppressWarnings(stats::loess(log_theta[ok_th] ~ amean[ok_th], span = 0.75, degree = 1L)), error = function(e) NULL)
    if (!is.null(lo)) { pr <- suppressWarnings(stats::predict(lo, newdata = data.frame(amean = amean)))
      pr <- as.numeric(pr); pr[!is.finite(pr)] <- trend[!is.finite(pr)]; trend <- pr }
  }
  w <- n_det / (n_det + d0_prior)
  log_theta_star <- ifelse(ok_th, w * log_theta + (1 - w) * trend, trend)
  theta_star <- exp(pmin(pmax(log_theta_star, -4), 6))

  # ---- stage 2: joint fit with theta fixed ----
  fit_one <- function(j) {
    xj <- otu[, j]; D <- xj > 0; n1 <- sum(D); theta <- theta_star[j]
    out <- list(status = "ok", p = NA_real_, lrt = NA_real_, alpha = rep(NA_real_, k), se = rep(NA_real_, k),
                theta = theta, n_det = n1, converged = FALSE)
    if (n1 < min_detected || n1 > n - min_detected) { out$status <- "no_variation_in_detection"; return(out) }
    if (qr(Xz[D, , drop = FALSE])$rank < pz && qr(Xz[!D, , drop = FALSE])$rank < pz) {
      # both detection groups rank-deficient for the design: still fine for the joint likelihood; continue
    }
    nll <- function(par, free_z) {
      a <- numeric(pz); a[free_z] <- par[seq_len(sum(free_z))]; g <- par[sum(free_z) + seq_len(pc)]
      psi <- stats::plogis(as.vector(Xz %*% a)); mu <- exp(lN + as.vector(Xc %*% g))
      l0 <- exp(theta * (log(theta) - log(theta + mu)))
      P0 <- (1 - psi) + psi * l0
      ll <- sum(log(pmax(P0[!D], 1e-300))) +
            sum(log(pmax(psi[D], 1e-300)) + stats::dnbinom(xj[D], size = theta, mu = mu[D], log = TRUE))
      -ll + sum(a^2) / (2 * tau_ridge^2)
    }
    # starts
    a0 <- rep(0, pz); a0[1L] <- stats::qlogis(min(max(mean(D), 0.05), 0.95))
    g0 <- rep(0, pc); g0[1L] <- if (n1 > 0) mean(log(xj[D] / depth[D])) else -10
    if (n1 >= pc + 3L) { cf <- tryCatch(stats::lm.fit(Xc[D, , drop = FALSE], log(xj[D] / depth[D]))$coefficients, error = function(e) NULL)
      if (!is.null(cf) && all(is.finite(cf))) g0 <- unname(cf) }
    bounds <- function(free_z) { nz <- sum(free_z)
      list(lo = c(rep(-12, nz), -30, rep(-10, pc - 1L)), hi = c(rep(12, nz), 10, rep(10, pc - 1L))) }
    fit <- function(free_z, start) {
      b <- bounds(free_z)
      best <- tryCatch(stats::nlminb(start, nll, free_z = free_z, lower = b$lo, upper = b$hi), error = function(e) NULL)
      best
    }
    free_full <- rep(TRUE, pz); free_null <- free_full; free_null[tcols] <- FALSE
    st_full <- c(a0, g0)
    full <- fit(free_full, st_full)
    if (is.null(full) || !is.finite(full$objective)) { out$status <- "fit_failed"; return(out) }
    # extra starts: perturb the first tested presence coefficient and the count intercept
    for (dz in c(-1.5, 1.5)) { st <- st_full; st[tcols[1L]] <- dz
      f2 <- fit(free_full, st); if (!is.null(f2) && is.finite(f2$objective) && f2$objective < full$objective) full <- f2 }
    pf <- full$par; a_full <- numeric(pz); a_full[] <- pf[seq_len(pz)]; g_full <- pf[pz + seq_len(pc)]
    st_null <- c(a_full[free_null], g_full)
    null <- fit(free_null, st_null)
    if (is.null(null) || !is.finite(null$objective)) { out$status <- "null_fit_failed"; return(out) }
    for (dg in c(-1, 1)) { st <- st_null; st[sum(free_null) + 1L] <- st[sum(free_null) + 1L] + dg
      f2 <- fit(free_null, st); if (!is.null(f2) && is.finite(f2$objective) && f2$objective < null$objective) null <- f2 }
    lrt <- 2 * (null$objective - full$objective)
    if (is.finite(lrt) && lrt < 0 && lrt > -1e-6) lrt <- 0
    if (!is.finite(lrt) || lrt < 0) { out$status <- "lrt_invalid"; return(out) }
    out$lrt <- lrt; out$p <- stats::pchisq(lrt, df = k, lower.tail = FALSE); out$converged <- TRUE
    out$alpha <- a_full[tcols]
    H <- tryCatch(stats::optimHess(pf, nll, free_z = free_full), error = function(e) NULL)
    if (!is.null(H)) { V <- tryCatch(solve(H), error = function(e) NULL)
      if (!is.null(V)) { d <- diag(V)[tcols]; out$se <- ifelse(is.finite(d) & d > 0, sqrt(d), NA_real_) } }
    out
  }
  fits <- if (n_cores > 1L) parallel::mclapply(seq_len(m), fit_one, mc.cores = n_cores) else lapply(seq_len(m), fit_one)

  status <- vapply(fits, function(f) f$status, character(1))
  alpha <- t(vapply(fits, function(f) f$alpha, numeric(k))); if (k == 1L) alpha <- matrix(alpha, ncol = 1L)
  se <- t(vapply(fits, function(f) f$se, numeric(k))); if (k == 1L) se <- matrix(se, ncol = 1L)
  colnames(alpha) <- colnames(se) <- cn
  out <- data.frame(pres_n_detected = n_det, pres_status = status,
                    pres_theta = theta_star, pres_lrt = vapply(fits, function(f) f$lrt, numeric(1)),
                    pres_p = vapply(fits, function(f) f$p, numeric(1)), stringsAsFactors = FALSE)
  for (i in seq_len(k)) {
    out[[paste0("pres_logor_", cn[i])]] <- alpha[, i]
    out[[paste0("pres_se_",    cn[i])]] <- se[, i]
  }
  if (k == 1L) { names(out) <- sub(paste0("_", cn[1L], "$"), "", names(out))
    out$pres_ci_lo <- out$pres_logor - 1.96 * out$pres_se; out$pres_ci_hi <- out$pres_logor + 1.96 * out$pres_se }
  attr(out, "theta_trend") <- trend; attr(out, "log_theta_raw") <- log_theta
  out
}
