# -----------------------------------------------------------------------------
# 05_enull.R -- shared calibration layer: empirical null with an estimated null SCALE
#
# PURSUE 0.2's enull_center() fits   null: b_j ~ N(delta, s_j^2),  alt: b_j ~ N(eta, omega^2 + s_j^2)
# and centres at delta. It trusts s_j. Every count model tried so far (corncob, the ZINB mu-test,
# the NB GLM) is anti-conservative on realistic counts because its s_j is too SMALL. Most features
# are null, so their spread measures how wrong s_j is: fit the null as N(delta, sigma0^2 s_j^2) with
# sigma0 >= 1 estimated too (Efron's empirical null, mixture form) and test against it. Same
# cross-feature principle as the shared link, applied to variance rather than to depth.
#
# enull2(b, se) -> list(delta, sigma0, pi0, p = two-sided p against N(delta, sigma0^2 se^2))
# sigma0 is floored at 1: the layer may only make a test MORE conservative (calibration first).
# -----------------------------------------------------------------------------
enull2 <- function(b, se, max_iter = 500L, tol = 1e-8, floor1 = TRUE, scale = TRUE) {
  ok <- is.finite(b) & is.finite(se) & se > 0 & se < 10 * stats::median(se[is.finite(se) & se > 0])
  bb <- b[ok]; ss <- se[ok]; n <- length(bb); p <- rep(NA_real_, length(b))
  if (n < 20L) return(list(delta = 0, sigma0 = 1, pi0 = NA, p = p, note = "too few features"))
  z <- bb / ss
  dens <- tryCatch(stats::density(bb, bw = "SJ"), error = function(e) stats::density(bb))
  delta <- dens$x[which.max(dens$y)]; s0 <- 1; eta <- mean(bb); om2 <- max(stats::var(bb), 1e-4); pi0 <- 0.85
  for (it in seq_len(max_iter)) {
    f0 <- stats::dnorm(bb, delta, s0 * ss); f1 <- stats::dnorm(bb, eta, sqrt(om2 + ss^2))
    w <- pi0 * f0 / (pi0 * f0 + (1 - pi0) * f1 + 1e-300)
    pi0n <- min(max(mean(w), 0.5), 0.995)                 # pi0 >= 0.5: the null is the majority
    dn <- sum(w * bb / ss^2) / sum(w / ss^2)
    s0n <- if (scale) sqrt(sum(w * (bb - dn)^2 / ss^2) / sum(w)) else 1; if (floor1) s0n <- max(s0n, 1)
    v1 <- om2 + ss^2; en <- sum((1 - w) * bb / v1) / sum((1 - w) / v1)
    om2n <- max(sum((1 - w) * ((bb - en)^2 - ss^2) / v1) / sum((1 - w) / v1), 1e-4)
    done <- abs(dn - delta) < tol && abs(s0n - s0) < tol && abs(pi0n - pi0) < tol
    delta <- dn; s0 <- s0n; eta <- en; om2 <- om2n; pi0 <- pi0n
    if (done) break
  }
  p[ok] <- 2 * stats::pnorm(-abs(bb - delta) / (s0 * ss))
  list(delta = delta, sigma0 = s0, pi0 = pi0, p = p)
}

# Robust centre only (scale fixed at 1). PURSUE 0.2's enull_center lets pi0 fall to 0.02, and then the
# "null" component can swap roles with the alternative: it1, house:R08, twinsuk_stool rep 3 -- pi0 0.07,
# delta -0.80 against a true-null median of -0.04, not converged, 79 false positives in one cell.
# Requiring pi0 >= 0.5 is the assumption the centring rests on anyway (most features are null).
enull_robust <- function(b, se) enull2(b, se, scale = FALSE)
