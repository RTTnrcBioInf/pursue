# -----------------------------------------------------------------------------
# centering.R -- empirical-null mixture centring of effect estimates across features
#
# For each tested coefficient, the per-feature estimates b_j with standard errors s_j
# are modelled as a two-component mixture:
#     null:  b_j ~ N(delta, s_j^2)
#     alt:   b_j ~ N(eta, omega^2 + s_j^2)
# fitted by EM. delta is the community-wide shift (the "typical" effect); the
# reported effect is b_j - delta, with SE^2 = s_j^2 + SE(delta)^2. pi0 is the
# estimated fraction of null features. This is Efron's empirical null applied to
# compositional centring; it replaces the median/mode/reference-set choices.
# -----------------------------------------------------------------------------

#' Empirical-null centre of a vector of estimates with known SEs
#' @param b estimates; @param se standard errors (same length); NA allowed.
#' @return list(delta, se_delta, pi0, eta, omega2, w = posterior null probability, n_used)
#' @keywords internal
enull_center <- function(b, se, max_iter = 500L, tol = 1e-8) {
  ok <- is.finite(b) & is.finite(se) & se > 0
  bb <- b[ok]; ss <- se[ok]; n <- length(bb)
  if (n < 10L) return(list(delta = 0, se_delta = 0, pi0 = NA_real_, eta = NA_real_, omega2 = NA_real_,
                           w = rep(NA_real_, length(b)), n_used = n, converged = FALSE,
                           note = "fewer than 10 usable features: no centring applied"))
  # start: density mode as null centre, moments for the alternative
  dens <- tryCatch(stats::density(bb, bw = "SJ"), error = function(e) stats::density(bb))
  delta <- dens$x[which.max(dens$y)]; eta <- mean(bb); omega2 <- max(stats::var(bb), 1e-4); pi0 <- 0.8
  converged <- FALSE
  for (it in seq_len(max_iter)) {
    f0 <- stats::dnorm(bb, delta, ss); f1 <- stats::dnorm(bb, eta, sqrt(omega2 + ss^2))
    w <- pi0 * f0 / (pi0 * f0 + (1 - pi0) * f1 + 1e-300)
    pi0_new <- min(max(mean(w), 0.02), 0.995)
    delta_new <- sum(w * bb / ss^2) / sum(w / ss^2)
    v1 <- omega2 + ss^2
    eta_new <- sum((1 - w) * bb / v1) / sum((1 - w) / v1)
    omega2_new <- max(sum((1 - w) * ((bb - eta_new)^2 - ss^2) / v1) / sum((1 - w) / v1), 1e-4)
    done <- abs(delta_new - delta) < tol && abs(pi0_new - pi0) < tol
    delta <- delta_new; eta <- eta_new; omega2 <- omega2_new; pi0 <- pi0_new
    if (done) { converged <- TRUE; break }
  }
  w_full <- rep(NA_real_, length(b)); w_full[ok] <- w
  list(delta = delta, se_delta = sqrt(1 / sum(w / ss^2)), pi0 = pi0, eta = eta, omega2 = omega2,
       w = w_full, n_used = n, converged = converged, note = NULL)
}

#' Centre every tested column of the abundance arm
#' @keywords internal
center_effects <- function(arm) {
  k <- arm$k
  delta <- numeric(k); se_delta <- numeric(k); pi0 <- numeric(k); notes <- character(k); w <- matrix(NA_real_, nrow(arm$beta), k)
  for (i in seq_len(k)) {
    ce <- enull_center(arm$beta[, i], arm$se[, i])
    delta[i] <- ce$delta; se_delta[i] <- ce$se_delta; pi0[i] <- ce$pi0; w[, i] <- ce$w
    notes[i] <- if (is.null(ce$note)) "" else ce$note
  }
  names(delta) <- names(se_delta) <- names(pi0) <- arm$coef_names
  list(delta = delta, se_delta = se_delta, pi0 = pi0, w_null = w, notes = notes,
       median_delta = apply(arm$beta, 2L, stats::median, na.rm = TRUE))
}
