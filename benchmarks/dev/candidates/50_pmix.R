# -----------------------------------------------------------------------------
# 50_pmix.R -- iteration 4: one count likelihood with a SHARED LEARNED mixing distribution
#
#   y_ij | u_ij ~ Poisson( N_i * exp(u_ij) ),   u_ij = a_j + x_i' b_j + s_j * z_ij,   z_ij ~ G
#
# G is ONE distribution for all features, learned nonparametrically on a grid (NPMLE-style EM
# across features); each feature has its own location a_j, effect b_j and spread s_j.
#
# Why this shape (it1-it3 evidence): the detection-only models waste the counts, and the parametric
# count models (NB GLM, ZINB mu-test, corncob) are anti-conservative because the per-feature
# variance model is wrong. Here depth enters only through Poisson thinning -- the one part of the
# data-generating process that is known (reads are sampled in proportion to depth) -- and
# everything biological, including zero inflation (a heavy left tail of G), is learned from all
# features at once. It is the shared-link idea of it1 extended from detection to the whole count
# distribution: the detection curve and the count distribution are both implied by (G, s_j).
#
# pln  : G fixed at a discretised standard normal (Poisson-lognormal) -- the ablation.
# pmix : G learned.
# Both: LRT on b_j against the robust empirical-null centre (enull_robust); pmix_raw against 0.
# -----------------------------------------------------------------------------
# Each grid atom is a negative-binomial kernel, not a point: a point grid makes the likelihood a
# sum of spikes once counts are large (Poisson resolution in log space ~1/sqrt(y) << grid spacing),
# which gave multiple optima and a null FPR of 0.21 on data simulated from the model itself. The NB
# kernel is a log-gamma smear of width ~ the grid spacing in u (kappa = 1/(s*dz)^2), so the
# effective mixing distribution is G convolved with a narrow log-gamma: smooth in every parameter,
# still closed form, no quadrature.
.pm_ll_parts <- function(y, eta, s, G) {
  n <- length(y); mu <- exp(outer(eta, s * G$z, "+")); kap <- min(max(1 / (s * G$dz)^2, 0.5), 1e5)
  C <- lgamma(y + kap) - lgamma(kap) + kap * log(kap)
  lkm <- log(kap + mu)
  ll <- C - kap * lkm + y * log(mu) - y * lkm + rep(log(G$w), each = n)
  mx <- ll[cbind(seq_len(n), max.col(ll, ties.method = "first"))]; E <- exp(ll - mx); rs <- rowSums(E)
  list(nll = -sum(mx + log(rs)), R = E / rs, mu = mu, kap = kap, lkm = lkm)
}
.pm_feature <- function(y, lN, X, G, start, fixed = integer(0), fixed_val = numeric(0), hess = FALSE) {
  n <- length(y); free <- setdiff(seq_len(ncol(X)), fixed); nf <- length(free)
  off <- lN + if (length(fixed)) drop(X[, fixed, drop = FALSE] %*% fixed_val) else 0
  Xf <- X[, free, drop = FALSE]; cache <- new.env()
  core <- function(th) {
    if (!is.null(cache$th) && identical(cache$th, th)) return(cache$v)
    v <- .pm_ll_parts(y, off + drop(Xf %*% th[seq_len(nf)]), exp(th[nf + 1L]), G); v$s <- exp(th[nf + 1L])
    cache$th <- th; cache$v <- v; v }
  obj <- function(th) { v <- core(th)$nll; if (is.finite(v)) v else 1e12 }
  grad <- function(th) { v <- core(th); k <- v$kap
    dl <- y - (y + k) * v$mu / (k + v$mu)                               # d ll / d log mu, n x G
    dk <- digamma(y + k) - digamma(k) + log(k) + 1 - v$lkm - (k + y) / (k + v$mu)   # d ll / d kappa
    dk_ds <- if (k > 0.5 && k < 1e5) -2 * k else 0
    -c(drop(crossprod(Xf, rowSums(v$R * dl))), sum(v$R * (dl * rep(v$s * G$z, each = n) + dk * dk_ds))) }
  st <- c(start[free], start[length(start)])
  f <- tryCatch(stats::nlminb(st, obj, grad, lower = c(rep(-Inf, nf), -4), upper = c(rep(Inf, nf), 2.5)), error = function(e) NULL)
  if (is.null(f) || f$objective >= 1e12) return(NULL)
  b <- numeric(ncol(X) + 1L); b[free] <- f$par[seq_len(nf)]; if (length(fixed)) b[fixed] <- fixed_val; b[ncol(X) + 1L] <- f$par[nf + 1L]
  out <- list(par = b, ll = -f$objective, grad = grad, obj = obj)
  if (hess) {                                                    # observed information by differencing the gradient
    h <- 1e-4; H <- vapply(seq_along(f$par), function(k) { e <- numeric(length(f$par)); e[k] <- h
      (grad(f$par + e) - grad(f$par - e)) / (2 * h) }, numeric(length(f$par)))
    out$V <- tryCatch(solve((H + t(H)) / 2), error = function(e) NULL)
  }
  out
}

# posterior grid weights of one fitted feature (for the G update)
.pm_post <- function(y, lN, X, G, par) {
  p <- ncol(X); colSums(.pm_ll_parts(y, lN + drop(X %*% par[seq_len(p)]), exp(par[p + 1L]), G)$R)
}

# per-sample score contributions at par (n x (p+1)): for the sandwich and the robust score test
.pm_scores <- function(y, lN, X, G, par) {
  p <- ncol(X); s <- exp(par[p + 1L]); v <- .pm_ll_parts(y, lN + drop(X %*% par[seq_len(p)]), s, G); k <- v$kap; n <- length(y)
  dl <- y - (y + k) * v$mu / (k + v$mu)
  dk <- digamma(y + k) - digamma(k) + log(k) + 1 - v$lkm - (k + y) / (k + v$mu); dk_ds <- if (k > 0.5 && k < 1e5) -2 * k else 0
  cbind(X * rowSums(v$R * dl), rowSums(v$R * (dl * rep(s * G$z, each = n) + dk * dk_ds)))
}

.pm_G0 <- function(ng = 41L) { z <- seq(-7, 4, length.out = ng); w <- stats::dnorm(z); list(z = z, w = w / sum(w), dz = z[2] - z[1]) }
.pm_standardise <- function(G) { m <- sum(G$w * G$z); v <- sqrt(sum(G$w * (G$z - m)^2)); list(z = (G$z - m) / v, w = G$w, dz = G$dz / v, m = m, v = v) }

.pm_memo <- new.env()
.pm_fit <- function(counts, meta, formula, tested_term, learn = TRUE, rounds = 3L, n_learn = 150L) {
  key <- list(counts, meta, formula, tested_term, learn)
  if (!is.null(.pm_memo[[as.character(learn)]]$key) && identical(.pm_memo[[as.character(learn)]]$key, key)) return(.pm_memo[[as.character(learn)]]$val)
  X <- stats::model.matrix(formula, meta); asg <- attr(X, "assign")
  tc <- which(asg == which(attr(stats::terms(formula), "term.labels") == tested_term)); k <- length(tc); p <- ncol(X)
  lN <- log(if (!is.null(meta$depth)) meta$depth else colSums(counts)); Y <- as.matrix(counts); storage.mode(Y) <- "double"
  m <- nrow(Y); n <- ncol(Y); ok <- which(rowSums(Y > 0) >= 3)
  P <- matrix(NA_real_, m, p + 1L)
  P[ok, 1] <- log((rowSums(Y[ok, , drop = FALSE]) + 0.5) / sum(exp(lN))); P[ok, -1] <- 0; P[ok, p + 1L] <- log(1.5)
  if (p > 1) P[ok, 2:p] <- 0
  G <- .pm_G0()
  fit_rows <- function(rows, G) for (j in rows) { r <- .pm_feature(Y[j, ], lN, X, G, P[j, ]); if (!is.null(r)) P[j, ] <<- r$par }
  if (learn) {
    set.seed(1L); L <- if (length(ok) > n_learn) sort(sample(ok, n_learn)) else ok
    for (r in seq_len(rounds)) {
      fit_rows(L, G)
      for (em in 1:25) { w <- Reduce(`+`, lapply(L, function(j) .pm_post(Y[j, ], lN, X, G, P[j, ])));
        G$w <- pmax(w / sum(w), 1e-10); G$w <- G$w / sum(G$w) }
      Gs <- .pm_standardise(G)                                   # keep G at mean 0, sd 1: move location/scale into a_j, s_j
      P[ok, 1] <- P[ok, 1] + exp(P[ok, p + 1L]) * Gs$m; P[ok, p + 1L] <- P[ok, p + 1L] + log(Gs$v); G <- Gs[c("z", "w", "dz")]
    }
  }
  ll_full <- rep(NA_real_, m); B <- rep(NA_real_, m); se <- rep(NA_real_, m)
  for (j in ok) { r <- .pm_feature(Y[j, ], lN, X, G, P[j, ], hess = TRUE); if (is.null(r)) next
    P[j, ] <- r$par; ll_full[j] <- r$ll; B[j] <- r$par[tc[1]]
    if (!is.null(r$V)) se[j] <- sqrt(max(r$V[tc[1], tc[1]], 0)) }
  delta <- if (k == 1L) enull_robust(B, se)$delta else rep(0, k)
  # sandwich Wald (model fit as estimating equation) and Boos' generalised (robust) score test
  n_ <- ncol(Y); cf <- n_ / (n_ - p - 1)
  sw <- rep(NA_real_, m); rs <- rep(NA_real_, m); rs_raw <- rep(NA_real_, m)
  for (j in ok) { if (!is.finite(ll_full[j])) next
    Hf <- tryCatch({ r <- .pm_feature(Y[j, ], lN, X, G, P[j, ], hess = TRUE); r$V }, error = function(e) NULL); if (is.null(Hf)) next
    S <- .pm_scores(Y[j, ], lN, X, G, P[j, ]); Vs <- Hf %*% crossprod(S) %*% Hf * cf
    sw[j] <- 2 * stats::pnorm(-abs(B[j] - delta[1]) / sqrt(max(Vs[tc[1], tc[1]], 1e-12)))
    for (w in 1:2) { val <- if (w == 1) delta else rep(0, k)
      r0 <- .pm_feature(Y[j, ], lN, X, G, P[j, ], fixed = tc, fixed_val = val); if (is.null(r0)) next
      S0 <- .pm_scores(Y[j, ], lN, X, G, r0$par); I0 <- crossprod(S0)             # outer-product information
      nu <- setdiff(seq_len(p + 1L), tc)
      H0 <- tryCatch({ h <- 1e-4; th <- r0$par; g <- function(t) -colSums(.pm_scores(Y[j, ], lN, X, G, t))
        vapply(seq_len(p + 1L), function(q) { e <- numeric(p + 1L); e[q] <- h; (g(th + e) - g(th - e)) / (2 * h) }, numeric(p + 1L)) }, error = function(e) NULL)
      if (is.null(H0)) next
      H0 <- (H0 + t(H0)) / 2; A <- H0[tc, nu, drop = FALSE] %*% solve(H0[nu, nu, drop = FALSE])
      Ueff <- S0[, tc, drop = FALSE] - S0[, nu, drop = FALSE] %*% t(A)                  # efficient score contributions
      U <- colSums(Ueff); Vu <- crossprod(Ueff) * cf
      st <- tryCatch(drop(t(U) %*% solve(Vu, U)), error = function(e) NA_real_)
      pv <- stats::pchisq(st, k, lower.tail = FALSE); if (w == 1) rs[j] <- pv else rs_raw[j] <- pv }
  }
  lrt <- function(val) vapply(seq_len(m), function(j) {
    if (!is.finite(ll_full[j])) return(NA_real_)
    r <- .pm_feature(Y[j, ], lN, X, G, P[j, ], fixed = tc, fixed_val = val); if (is.null(r)) return(NA_real_)
    stats::pchisq(max(2 * (ll_full[j] - r$ll), 0), k, lower.tail = FALSE) }, numeric(1))
  val <- list(feature = rownames(counts), p = lrt(delta), p_raw = lrt(rep(0, k)), p_sw = sw, p_rs = rs, p_rs_raw = rs_raw,
              est = B - delta[1], b = B, se = se, G = G, s = exp(P[, p + 1L]))
  .pm_memo[[as.character(learn)]] <- list(key = key, val = val); val
}

register_candidate("pmix", function(counts, meta, formula, tested_term) {
  v <- .pm_fit(counts, meta, formula, tested_term, learn = TRUE); data.frame(feature = v$feature, p = v$p, estimate = v$est)
}, notes = "it4: Poisson sampling x shared learned mixing distribution G (NPMLE across features), per-feature location/effect/spread; LRT vs robust centre")
register_candidate("pmix_raw", function(counts, meta, formula, tested_term) {
  v <- .pm_fit(counts, meta, formula, tested_term, learn = TRUE); data.frame(feature = v$feature, p = v$p_raw)
}, notes = "it4: pmix tested against 0")
register_candidate("pln", function(counts, meta, formula, tested_term) {
  v <- .pm_fit(counts, meta, formula, tested_term, learn = FALSE); data.frame(feature = v$feature, p = v$p, estimate = v$est)
}, notes = "it4 ablation: Poisson-lognormal (G fixed standard normal), same test")
register_candidate("pln_raw", function(counts, meta, formula, tested_term) {
  v <- .pm_fit(counts, meta, formula, tested_term, learn = FALSE); data.frame(feature = v$feature, p = v$p_raw)
}, notes = "it4: Poisson-lognormal LRT against 0")
register_candidate("pln_sw", function(counts, meta, formula, tested_term) {
  v <- .pm_fit(counts, meta, formula, tested_term, learn = FALSE); data.frame(feature = v$feature, p = v$p_sw, estimate = v$est)
}, notes = "it4: Poisson-lognormal, sandwich Wald against the robust centre")
register_candidate("pln_rs", function(counts, meta, formula, tested_term) {
  v <- .pm_fit(counts, meta, formula, tested_term, learn = FALSE); data.frame(feature = v$feature, p = v$p_rs, estimate = v$est)
}, notes = "it4: Poisson-lognormal, robust (Boos) score test against the robust centre")
register_candidate("pln_rs_raw", function(counts, meta, formula, tested_term) {
  v <- .pm_fit(counts, meta, formula, tested_term, learn = FALSE); data.frame(feature = v$feature, p = v$p_rs_raw)
}, notes = "it4: Poisson-lognormal, robust score test against 0")
