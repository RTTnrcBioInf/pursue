# -----------------------------------------------------------------------------
# 10_sharedlink.R -- iteration 1: the shared learned depth link  (claude/rd-charter.md, backlog)
#
#   P(detected_ij) = f( a_j + log N_i + x_i' b_j ),   f monotone, ONE curve for all features
#
# Why: under strong depth confounding a single feature cannot tell a depth effect from a group
# effect -- a free per-feature depth slope (logistic_depth) is calibrated but powerless there, and
# a fixed parametric shape (PURSUE 0.2's NB, cloglog) keeps power only when the shape happens to be
# right (house) and makes false positives when it is not (mid). But depth acts on every taxon
# through the same sampling mechanism and the exposure on a sparse few, so the SHAPE of the
# detection-depth curve is learnable from all features at once. log N enters with coefficient 1
# (expected reads scale with depth: sampling, not an assumption); only f is learned. With f known,
# b_j stays identified even when group and depth are collinear within the feature.
#
# The assumption that remains: detection depends on abundance level and depth only through their
# sum (a shared mechanism). Weaker than a fixed parametric link, and it is what this tests.
#
# f = plogis(g), g(e) = c0 + sum_l d_l T_l(e), d_l >= 0, T_l = tail sums of a cubic B-spline
# basis -- monotone by construction. Alternates: link from intercept-only fits -> per-feature
# full fits -> link refitted -> per-feature full fits with Fisher-information SEs -> tested
# coefficients centred at PURSUE 0.2's empirical-null centre -> LRT against the centre.
# -----------------------------------------------------------------------------
.sl_memo <- new.env()

.sl_link <- function(E, Dm, K = 8L, max_pts = 60000L) {
  e <- as.vector(E); d <- as.vector(Dm)
  if (length(e) > max_pts) { i <- sample.int(length(e), max_pts); e <- e[i]; d <- d[i] }
  rng <- range(e); inner <- stats::quantile(e, probs = seq(0, 1, length.out = K - 2L)[-c(1L, K - 2L)], names = FALSE)
  kn <- c(rep(rng[1], 4L), inner, rep(rng[2], 4L))
  basis <- function(x, deriv = 0L) {
    Bm <- splines::splineDesign(kn, pmin(pmax(x, rng[1]), rng[2]), ord = 4L, derivs = rep(deriv, length(x)), outer.ok = TRUE)
    U <- upper.tri(diag(ncol(Bm)), diag = TRUE); storage.mode(U) <- "double"
    (Bm %*% t(U))[, -1L, drop = FALSE]           # tail sums T_2..T_K (monotone 0 -> 1)
  }
  Tm <- basis(e)
  nll <- function(p) { g <- p[1] + drop(Tm %*% p[-1]); P <- stats::plogis(g); -sum(d * log(P + 1e-12) + (1 - d) * log(1 - P + 1e-12)) }
  gr  <- function(p) { g <- p[1] + drop(Tm %*% p[-1]); r <- d - stats::plogis(g); -c(sum(r), drop(crossprod(Tm, r))) }
  fit <- stats::nlminb(c(stats::qlogis(min(max(mean(d), 0.01), 0.99)), rep(0.5, ncol(Tm))), nll, gr,
                       lower = c(-Inf, rep(0, ncol(Tm))))
  eg <- seq(rng[1], rng[2], length.out = 2001L)
  g <- fit$par[1] + drop(basis(eg) %*% fit$par[-1]); gp <- drop(basis(eg, 1L) %*% fit$par[-1]); gpp <- drop(basis(eg, 2L) %*% fit$par[-1])
  P <- pmin(pmax(stats::plogis(g), 1e-9), 1 - 1e-9)
  list(fP = function(x) stats::approx(eg, P, xout = x, rule = 2)$y,
       fgp = function(x) stats::approx(eg, gp, xout = x, yleft = 0, yright = 0)$y,
       fgpp = function(x) stats::approx(eg, gpp, xout = x, yleft = 0, yright = 0)$y,
       range = rng, par = fit$par)
}

.sl_feature <- function(d, lN, X, L, start, fixed = integer(0), fixed_val = numeric(0)) {
  free <- setdiff(seq_len(ncol(X)), fixed)
  off <- lN + if (length(fixed)) drop(X[, fixed, drop = FALSE] %*% fixed_val) else 0
  Xf <- X[, free, drop = FALSE]
  nll <- function(b) { P <- L$fP(off + drop(Xf %*% b)); -sum(d * log(P) + (1 - d) * log(1 - P)) }
  gr  <- function(b) { e <- off + drop(Xf %*% b); -drop(crossprod(Xf, (d - L$fP(e)) * L$fgp(e))) }
  f <- tryCatch(stats::nlminb(start[free], nll, gr), error = function(e) NULL)
  if (is.null(f) || !is.finite(f$objective)) return(NULL)
  b <- numeric(ncol(X)); b[free] <- f$par; if (length(fixed)) b[fixed] <- fixed_val
  list(b = b, ll = -f$objective)
}

.sl_fit <- function(counts, meta, formula, tested_term) {
  key <- list(counts, meta, formula, tested_term)
  if (!is.null(.sl_memo$key) && identical(.sl_memo$key, key)) return(.sl_memo$val)
  X <- stats::model.matrix(formula, meta); asg <- attr(X, "assign")
  tc <- which(asg == which(attr(stats::terms(formula), "term.labels") == tested_term)); k <- length(tc)
  depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts)
  lN <- log(depth); D <- counts > 0; m <- nrow(D); n <- ncol(D)
  ok <- which(rowSums(D) >= 3 & rowSums(D) <= n - 3)
  pbar <- rowMeans(D[ok, , drop = FALSE])
  B <- matrix(0, m, ncol(X)); B[ok, 1] <- log(-log(1 - pmin(pmax(pbar, 0.01), 0.99))) - mean(lN)
  eta <- function(rows) B[rows, , drop = FALSE] %*% t(X) + matrix(lN, length(rows), n, byrow = TRUE)
  full_pass <- function(L) for (j in ok) { r <- .sl_feature(as.numeric(D[j, ]), lN, X, L, B[j, ]); if (!is.null(r)) B[j, ] <<- r$b }
  L <- .sl_link(eta(ok), D[ok, , drop = FALSE])          # round 1: intercept-only eta
  full_pass(L)
  L <- .sl_link(eta(ok), D[ok, , drop = FALSE])          # round 2: eta from full fits
  full_pass(L)
  ll_full <- rep(NA_real_, m); se <- matrix(NA_real_, m, k)
  for (j in ok) {
    e <- lN + drop(X %*% B[j, ]); P <- L$fP(e); w <- P * (1 - P) * L$fgp(e)^2
    ll_full[j] <- sum(D[j, ] * log(P) + (1 - D[j, ]) * log(1 - P))
    I <- crossprod(X * sqrt(w)); V <- tryCatch(solve(I), error = function(e) NULL)
    if (!is.null(V)) se[j, ] <- sqrt(pmax(diag(V)[tc], 0))
  }
  delta <- vapply(seq_len(k), function(c) {
    ce <- tryCatch(get("enull_center", envir = asNamespace("PURSUE"))(B[, tc[c]], se[, c]), error = function(e) NULL)
    if (is.null(ce)) stats::median(B[ok, tc[c]]) else ce$delta }, numeric(1))
  lrt <- function(val) vapply(seq_len(m), function(j) {
    if (!(j %in% ok) || !is.finite(ll_full[j])) return(NA_real_)
    r <- .sl_feature(as.numeric(D[j, ]), lN, X, L, B[j, ], fixed = tc, fixed_val = val)
    if (is.null(r)) return(NA_real_)
    stats::pchisq(max(2 * (ll_full[j] - r$ll), 0), k, lower.tail = FALSE) }, numeric(1))
  delta_rc <- vapply(seq_len(k), function(c) enull_robust(B[, tc[c]], se[, c])$delta, numeric(1))
  val <- list(feature = rownames(counts), p = lrt(delta), p_raw = lrt(rep(0, k)), p_rc = lrt(delta_rc), delta_rc = delta_rc,
              L = L, B = B, ok = ok, X = X, tc = tc,
              est = B[, tc[1]] - delta[1], delta = delta, link = L$par,
              b = B[, tc[1]], se = se[, 1])            # raw estimate + Fisher SE, for calibration layers
  .sl_memo$key <- key; .sl_memo$val <- val; val
}

register_candidate("sharedlink", function(counts, meta, formula, tested_term) {
  v <- .sl_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p, estimate = v$est)
}, notes = "it1: detection through ONE monotone link learned across features, log-depth at slope 1, effect tested against the empirical-null centre")

register_candidate("sharedlink_raw", function(counts, meta, formula, tested_term) {
  v <- .sl_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_raw)
}, notes = "it1: as sharedlink, effect tested against 0 (no centring) -- isolates what the centring does")

register_candidate("sharedlink_rc", function(counts, meta, formula, tested_term) {
  v <- .sl_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_rc, estimate = v$b - v$delta_rc[1])
}, notes = "it3: sharedlink with the robust empirical-null centre (pi0 >= 0.5)")
