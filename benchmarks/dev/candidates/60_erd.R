# -----------------------------------------------------------------------------
# 60_erd.R -- iteration 5: expected rarefied detection (ERD)
#
# Rarefying every sample to a common depth D removes depth from detection exactly (for read
# sampling), but throws reads away at random. Its EXPECTATION does not: the probability that
# feature j would still be seen after rarefying sample i to D reads is hypergeometric,
#     f_ij(D) = 1 - C(N_i - y_ij, D) / C(N_i, D),
# a deterministic, bounded transform of the counts. For every sample with N_i >= D its expectation
# depends only on the sample's true composition, not on N_i, so depth confounding is removed by
# construction rather than modelled -- no depth covariate, no link to learn. Deeper samples just
# give a less noisy f. f is ~ D*p for rare features (a relative-abundance scale) and saturates at 1
# for common ones (a detection scale): a presence/abundance compromise set by D.
#
# D = the smallest library in the cell (every sample usable). Tests on f:
#   erd_lm  : linear model on f, HC3 sandwich t-test against 0 (relative change; no centring)
#   erd_qb  : quasi-binomial logit GLM on f, HC3 Wald; coefficient centred at the robust empirical
#             null (on the logit scale a compositional factor c is ~ a common shift log c for the
#             many rare features)
#   erd_qb_raw : as erd_qb, against 0
# -----------------------------------------------------------------------------
.erd_f <- function(counts, depth, D = min(depth)) {
  Y <- as.matrix(counts); N <- matrix(depth, nrow(Y), ncol(Y), byrow = TRUE)
  lr <- (lchoose(N - Y, D) - lchoose(N, D)); lr[!is.finite(lr)] <- -Inf     # N - y < D: always detected
  1 - exp(lr)
}
.hc3 <- function(X, r, h, bread) { meat <- crossprod(X * (r / (1 - pmin(h, 0.99)))); bread %*% meat %*% bread }

.erd_memo <- new.env()
.erd_fit <- function(counts, meta, formula, tested_term) {
  key <- list(counts, meta, formula, tested_term)
  if (!is.null(.erd_memo$key) && identical(.erd_memo$key, key)) return(.erd_memo$val)
  X <- stats::model.matrix(formula, meta); asg <- attr(X, "assign")
  tc <- which(asg == which(attr(stats::terms(formula), "term.labels") == tested_term))[1]
  depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts)
  Fm <- .erd_f(counts, depth); m <- nrow(Fm); n <- ncol(Fm); p <- ncol(X)
  XtXi <- solve(crossprod(X)); H <- rowSums((X %*% XtXi) * X)
  b_lm <- se_lm <- b_qb <- se_qb <- b_qp <- se_qp <- rep(NA_real_, m)
  for (j in seq_len(m)) {
    f <- Fm[j, ]; if (sum(f > 0) < 3 || stats::var(f) == 0) next
    beta <- XtXi %*% crossprod(X, f); r <- f - drop(X %*% beta)
    V <- .hc3(X, r, H, XtXi); b_lm[j] <- beta[tc]; se_lm[j] <- sqrt(max(V[tc, tc], 0))
    g <- tryCatch(suppressWarnings(stats::glm.fit(X, f, family = stats::quasibinomial())), error = function(e) NULL)
    if (is.null(g) || !g$converged) next
    mu <- g$fitted.values; w <- mu * (1 - mu); Xw <- X * sqrt(w)
    Bi <- tryCatch(solve(crossprod(Xw)), error = function(e) NULL); if (is.null(Bi)) next
    hq <- rowSums((Xw %*% Bi) * Xw); meat <- crossprod(X * ((f - mu) / (1 - pmin(hq, 0.99))))
    Vq <- Bi %*% meat %*% Bi; b_qb[j] <- g$coefficients[tc]; se_qb[j] <- sqrt(max(Vq[tc, tc], 0))
    # log link: a compositional factor c moves E[f] by ~ log c for every rare feature -> centrable
    h <- tryCatch(suppressWarnings(stats::glm.fit(X, f, family = stats::quasipoisson())), error = function(e) NULL)
    if (is.null(h) || !h$converged) next
    mu <- h$fitted.values; Xw <- X * sqrt(mu); Bi <- tryCatch(solve(crossprod(Xw)), error = function(e) NULL); if (is.null(Bi)) next
    hq <- rowSums((Xw %*% Bi) * Xw); meat <- crossprod(X * ((f - mu) / (1 - pmin(hq, 0.99))))
    Vp <- Bi %*% meat %*% Bi; b_qp[j] <- h$coefficients[tc]; se_qp[j] <- sqrt(max(Vp[tc, tc], 0))
  }
  dq <- enull_robust(b_qb, se_qb)$delta; dp <- enull_robust(b_qp, se_qp)$delta; df <- n - p
  val <- list(feature = rownames(counts),
              p_lm = 2 * stats::pt(-abs(b_lm / se_lm), df),
              p_qb = 2 * stats::pt(-abs((b_qb - dq) / se_qb), df), p_qb_raw = 2 * stats::pt(-abs(b_qb / se_qb), df),
              p_qp = 2 * stats::pt(-abs((b_qp - dp) / se_qp), df), est_qp = b_qp - dp,
              est = b_qb - dq, D = min(depth))
  .erd_memo$key <- key; .erd_memo$val <- val; val
}
register_candidate("erd_lm", function(counts, meta, formula, tested_term) {
  v <- .erd_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_lm)
}, notes = "it5: expected rarefied detection at the smallest library, LM + HC3 t-test (relative)")
register_candidate("erd_qb", function(counts, meta, formula, tested_term) {
  v <- .erd_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_qb, estimate = v$est)
}, notes = "it5: expected rarefied detection, quasi-binomial logit + HC3, robust empirical-null centre")
register_candidate("erd_qp", function(counts, meta, formula, tested_term) {
  v <- .erd_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_qp, estimate = v$est_qp)
}, notes = "it5: expected rarefied detection, quasi-Poisson log link + HC3, robust empirical-null centre (absolute)")
register_candidate("erd_qb_raw", function(counts, meta, formula, tested_term) {
  v <- .erd_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_qb_raw)
}, notes = "it5: erd_qb against 0")

# --- generalisation: expected rarefied h(count), E[h(Y_D)], Y_D ~ Hypergeometric(N_i, y_ij, D) ------
# h(y) = 1{y > 0} is ERD above. h(y) = log(1 + y) -- the expected rarefied log count (ERL) -- keeps
# ERD's depth invariance but does not saturate: ~ log(1 + D p) for common features (abundance
# information that detection throws away), ~ detection for rare ones. Exact hypergeometric sums
# where D p <= 60; second-order delta method beyond (relative error < 1e-4 there).
.erl_f <- function(counts, depth, D = min(depth), kmax = 150L) {
  Y <- as.matrix(counts); N <- matrix(depth, nrow(Y), ncol(Y), byrow = TRUE)
  mu <- D * Y / N; out <- matrix(0, nrow(Y), ncol(Y))
  big <- mu > 60; if (any(big)) { pp <- (Y / N)[big]; v <- D * pp * (1 - pp) * (N[big] - D) / pmax(N[big] - 1, 1)
    out[big] <- log1p(mu[big]) - v / (2 * (1 + mu[big])^2) }
  sm <- which(!big & Y > 0)
  if (length(sm)) { y <- Y[sm]; nn <- N[sm] - Y[sm]; acc <- numeric(length(sm))
    for (k in 1:kmax) { live <- k <= y & k <= D; if (!any(live)) break
      acc[live] <- acc[live] + stats::dhyper(k, y[live], nn[live], D) * log1p(k) }
    out[sm] <- acc }
  out
}
# Clustered samples (repeated measures: meta$subject with repeats) get a cluster-robust variance
# instead of HC3: CR1 (G/(G-1)) with t on G-1 df, G = number of clusters. Every benchmark method
# ignores subject in R23 and is anti-conservative there (full benchmark, house: FDR 0.07-0.17).
.cluster_of <- function(meta) { s <- meta$subject; if (is.null(s)) return(NULL); s <- as.factor(s)
  if (nlevels(s) < length(s) && nlevels(s) >= 6L) s else NULL }
.lm_cr_fit <- function(Fm, X, tc, cl) {
  XtXi <- solve(crossprod(X)); a <- drop(XtXi[tc, , drop = FALSE] %*% t(X))          # influence weights of the tested coefficient
  B <- Fm %*% t(XtXi %*% t(X)); b <- B[, tc]; R <- Fm - B %*% t(X)
  Z <- stats::model.matrix(~ cl - 1); G <- ncol(Z)
  S <- (R * rep(a, each = nrow(R))) %*% Z; v <- rowSums(S^2) * G / (G - 1)
  ok <- rowSums(Fm > 0) >= 3 & apply(Fm, 1, stats::var) > 0
  b[!ok] <- NA; list(b = b, se = ifelse(ok, sqrt(v), NA_real_), df = G - 1L)
}
.lm_hc3_fit <- function(Fm, X, tc, cl = NULL) {
  if (!is.null(cl)) return(.lm_cr_fit(Fm, X, tc, cl))
  XtXi <- solve(crossprod(X)); H <- rowSums((X %*% XtXi) * X); m <- nrow(Fm); b <- se <- rep(NA_real_, m)
  for (j in seq_len(m)) { f <- Fm[j, ]; if (sum(f > 0) < 3 || stats::var(f) == 0) next
    beta <- XtXi %*% crossprod(X, f); r <- f - drop(X %*% beta); V <- .hc3(X, r, H, XtXi)
    b[j] <- beta[tc]; se[j] <- sqrt(max(V[tc, tc], 0)) }
  list(b = b, se = se, df = ncol(Fm) - ncol(X))
}
.erl_memo <- new.env()
.erl_fit <- function(counts, meta, formula, tested_term) {
  key <- list(counts, meta, formula, tested_term)
  if (!is.null(.erl_memo$key) && identical(.erl_memo$key, key)) return(.erl_memo$val)
  X <- stats::model.matrix(formula, meta); asg <- attr(X, "assign")
  tc <- which(asg == which(attr(stats::terms(formula), "term.labels") == tested_term))[1]
  depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts)
  r <- .lm_hc3_fit(.erl_f(counts, depth), X, tc); d <- enull_robust(r$b, r$se)$delta
  val <- list(feature = rownames(counts), p = 2 * stats::pt(-abs((r$b - d) / r$se), r$df), p_raw = 2 * stats::pt(-abs(r$b / r$se), r$df), est = r$b - d)
  .erl_memo$key <- key; .erl_memo$val <- val; val
}
register_candidate("erl_lm", function(counts, meta, formula, tested_term) {
  v <- .erl_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p, estimate = v$est)
}, notes = "it5: expected rarefied log count at the smallest library, LM + HC3, robust empirical-null centre")
register_candidate("erl_lm_raw", function(counts, meta, formula, tested_term) {
  v <- .erl_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_raw)
}, notes = "it5: erl_lm against 0 (relative)")

# rank test on ERD for a lone binary factor (otherwise falls back to erd_lm)
register_candidate("erd_wx", function(counts, meta, formula, tested_term) {
  g <- meta[[tested_term]]; tl <- attr(stats::terms(formula), "term.labels")
  if (length(tl) != 1L || length(unique(g)) != 2L) { v <- .erd_fit(counts, meta, formula, tested_term); return(data.frame(feature = v$feature, p = v$p_lm)) }
  depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts); Fm <- .erd_f(counts, depth); g1 <- g == sort(unique(g))[2]
  p <- apply(Fm, 1, function(f) if (sum(f > 0) < 3) NA_real_ else suppressWarnings(stats::wilcox.test(f[g1], f[!g1], exact = FALSE)$p.value))
  data.frame(feature = rownames(counts), p = p)
}, notes = "it5: expected rarefied detection, Wilcoxon rank-sum (binary designs)")

# --- adaptive depth per feature (label-blind) ------------------------------------------------------
# One D for all features wastes most of them: common features saturate (f ~ 1 everywhere, so an
# abundance change is invisible), rare ones barely register at the smallest library. The transform's
# sensitivity to a proportional change in p, df/dlog p = D p (1-p)^(D-1), peaks at D p ~ 1. So each
# feature gets its own D_j = 1 / (median relative abundance among samples where it is seen), clamped
# to [10, D_max]: an abundance change in a common taxon becomes a detection change at a small D; a
# rare taxon is read at the deepest D the cell allows. D_j uses pooled data only, never the labels,
# so selecting it cannot bias the test (the same argument as independent filtering).
#   erd_a   : D_max = smallest library (every sample kept)
#   erd_a10 : D_max = 10th percentile of depth; samples below D_j dropped for that feature
.erd_adapt <- function(counts, meta, formula, tested_term, qmax = 0) {
  X <- stats::model.matrix(formula, meta); asg <- attr(X, "assign")
  tc <- which(asg == which(attr(stats::terms(formula), "term.labels") == tested_term))[1]
  depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts)
  Y <- as.matrix(counts); m <- nrow(Y); Dmax <- as.numeric(stats::quantile(depth, qmax, type = 1)); P <- sweep(Y, 2, depth, "/")
  b <- se <- dfs <- rep(NA_real_, m)
  for (j in seq_len(m)) {
    pj <- P[j, ]; if (sum(pj > 0) < 3) next
    Dj <- floor(min(max(1 / stats::median(pj[pj > 0]), 10), Dmax)); keep <- depth >= Dj
    f <- 1 - exp(lchoose(depth[keep] - Y[j, keep], Dj) - lchoose(depth[keep], Dj)); f[!is.finite(f)] <- 1
    Xk <- X[keep, , drop = FALSE]; if (sum(f > 0) < 3 || stats::var(f) == 0 || qr(Xk)$rank < ncol(X)) next
    XtXi <- solve(crossprod(Xk)); beta <- XtXi %*% crossprod(Xk, f); r <- f - drop(Xk %*% beta)
    V <- .hc3(Xk, r, rowSums((Xk %*% XtXi) * Xk), XtXi); b[j] <- beta[tc]; se[j] <- sqrt(max(V[tc, tc], 0)); dfs[j] <- sum(keep) - ncol(X)
  }
  data.frame(feature = rownames(counts), p = 2 * stats::pt(-abs(b / se), dfs))
}
register_candidate("erd_a", function(counts, meta, formula, tested_term) .erd_adapt(counts, meta, formula, tested_term, 0),
  notes = "it5: ERD with a per-feature depth D_j ~ 1/typical relative abundance, capped at the smallest library; LM + HC3")
register_candidate("erd_a10", function(counts, meta, formula, tested_term) .erd_adapt(counts, meta, formula, tested_term, 0.10),
  notes = "it5: as erd_a, capped at the 10th depth percentile (shallower samples dropped per feature)")

# --- it6: a rank test that tolerates unequal spreads, and a design-chosen D ---------------------------
# it5: erd_wx (Wilcoxon) had the most power of any calibrated candidate (rel. 0.82) but FPR 0.075 at
# R19 -- depth confounding makes the two groups' spreads of f differ, and Wilcoxon assumes they do
# not. Brunner-Munzel is the rank test built for exactly that (Behrens-Fisher) case.
.bm_test <- function(x, y) {
  n1 <- length(x); n2 <- length(y); r <- rank(c(x, y)); r1 <- r[seq_len(n1)]; r2 <- r[n1 + seq_len(n2)]
  i1 <- rank(x); i2 <- rank(y); m1 <- mean(r1); m2 <- mean(r2)
  v1 <- sum((r1 - i1 - m1 + (n1 + 1) / 2)^2) / (n1 - 1); v2 <- sum((r2 - i2 - m2 + (n2 + 1) / 2)^2) / (n2 - 1)
  den <- n1 * v1 + n2 * v2; if (!is.finite(den) || den <= 0) return(NA_real_)
  W <- n1 * n2 * (m2 - m1) / ((n1 + n2) * sqrt(den))
  df <- den^2 / ((n1 * v1)^2 / (n1 - 1) + (n2 * v2)^2 / (n2 - 1))
  2 * stats::pt(-abs(W), df)
}
# D chosen from the design alone (labels and depths, never counts), so the choice cannot bias any
# feature's test. For a rare feature with rate p, detections at depth D are ~ D p per kept sample,
# and a two-group contrast's information scales with the harmonic size of the kept groups:
# H(D) = D * n1(D) n0(D) / (n1(D) + n0(D)), samples with N_i < D dropped. Unconfounded, raising D
# above the smallest library pays; under confounding it empties the shallow group and does not.
.erd_Dopt <- function(depth, g) {
  cand <- sort(unique(stats::quantile(depth, c(0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50), type = 1)))
  H <- vapply(cand, function(D) { k <- depth >= D; n1 <- sum(k & g); n0 <- sum(k & !g)
    if (n1 < 5 || n0 < 5) return(-Inf); D * n1 * n0 / (n1 + n0) }, numeric(1))
  cand[which.max(H)]
}
.erd_two <- function(counts, meta, formula, tested_term, opt = FALSE, test = c("lm", "bm")) {
  test <- match.arg(test)
  g <- meta[[tested_term]]; tl <- attr(stats::terms(formula), "term.labels")
  binary <- length(tl) == 1L && length(unique(g)) == 2L
  depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts)
  if (!binary) { if (test == "bm" && !opt) { v <- .erd_fit(counts, meta, formula, tested_term); return(data.frame(feature = v$feature, p = v$p_lm)) }
    D <- min(depth) } else { g1 <- g == sort(unique(g))[2]; D <- if (opt) .erd_Dopt(depth, g1) else min(depth) }
  keep <- depth >= D; Fm <- .erd_f(counts[, keep, drop = FALSE], depth[keep], D)
  if (test == "bm" && binary) {
    gk <- g1[keep]; p <- apply(Fm, 1, function(f) if (sum(f > 0) < 3) NA_real_ else .bm_test(f[!gk], f[gk]))
  } else {
    X <- stats::model.matrix(formula, meta[keep, , drop = FALSE]); asg <- attr(X, "assign")
    tc <- which(asg == which(tl == tested_term))[1]; r <- .lm_hc3_fit(Fm, X, tc); p <- 2 * stats::pt(-abs(r$b / r$se), r$df)
  }
  data.frame(feature = rownames(counts), p = p)
}
register_candidate("erd_bm", function(counts, meta, formula, tested_term) .erd_two(counts, meta, formula, tested_term, FALSE, "bm"),
  notes = "it6: ERD at the smallest library, Brunner-Munzel (binary designs; else erd_lm)")
register_candidate("erd_opt", function(counts, meta, formula, tested_term) .erd_two(counts, meta, formula, tested_term, TRUE, "lm"),
  notes = "it6: ERD at a design-chosen D (max D*harmonic kept group size; labels+depth only), LM + HC3")
register_candidate("erd_opt_bm", function(counts, meta, formula, tested_term) .erd_two(counts, meta, formula, tested_term, TRUE, "bm"),
  notes = "it6: ERD at the design-chosen D, Brunner-Munzel")

# --- it7: multi-depth ERD -- the whole expected rarefaction curve, invariance kept -----------------
# ERD at the smallest library is invariant but wastes deep samples' rare detections (it5: B_prev
# 25 vs 38 TP for logistic_depth). For every D, the mean of f(D) over the samples deep enough to be
# rarefied to D is depth-invariant; so the effect on E[f(D)] can be estimated at several depths,
# each from its own eligible samples, and combined. Each depth's LM coefficient has influence
# function psi_ik = [(X_k'X_k)^-1 x_i r_ik / (1 - h_ik)] (HC3-type) for eligible i, 0 otherwise; the
# joint covariance of the coefficients across depths is sum_i psi_i psi_i'. Combination: GLS
# was tried first and was anti-conservative everywhere (R00 FPR 0.084; R19 0.82 with thin layers):
# weights estimated from S exploit its noise. Now: equal weights on the coefficients, one variance
# estimated, and layers with any leverage > 0.1 skipped. Depth grid from
# design quantiles (0, 10, 25, 50%) -- labels and counts not used to choose it.
.erd_md <- function(counts, meta, formula, tested_term, qs = c(0, 0.10, 0.25, 0.50), min_n = 10L) {
  X <- stats::model.matrix(formula, meta); asg <- attr(X, "assign")
  tc <- which(asg == which(attr(stats::terms(formula), "term.labels") == tested_term))[1]
  depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts)
  Ds <- sort(unique(floor(stats::quantile(depth, qs, type = 1)))); n <- ncol(counts); m <- nrow(counts)
  layers <- list()
  for (D in Ds) { k <- depth >= D; Xk <- X[k, , drop = FALSE]
    if (sum(k) < max(min_n, ncol(X) + 3) || qr(Xk)$rank < ncol(X)) next
    XtXi <- solve(crossprod(Xk)); A <- (XtXi %*% t(Xk))[tc, ]                 # row of (X'X)^-1 X' for the tested coefficient
    if (max(rowSums((Xk %*% XtXi) * Xk)) > 0.1) next                           # a layer where any sample has leverage > 0.1 (e.g. < 10 per group) is too thin for HC3
    layers[[length(layers) + 1L]] <- list(D = D, k = which(k), A = A, Xk = Xk, XtXi = XtXi, h = rowSums((Xk %*% XtXi) * Xk),
                                          F = .erd_f(counts[, k, drop = FALSE], depth[k], D)) }
  L <- length(layers); z <- rep(NA_real_, m)
  for (j in seq_len(m)) {
    Psi <- matrix(0, n, L); b <- rep(NA_real_, L)
    for (l in seq_len(L)) { f <- layers[[l]]$F[j, ]; if (sum(f > 0) < 3 || stats::var(f) == 0) next
      beta <- layers[[l]]$XtXi %*% crossprod(layers[[l]]$Xk, f); r <- f - drop(layers[[l]]$Xk %*% beta)
      b[l] <- beta[tc]; Psi[layers[[l]]$k, l] <- layers[[l]]$A * r / (1 - pmin(layers[[l]]$h, 0.99)) }
    ok <- which(is.finite(b)); if (!length(ok)) next
    S <- crossprod(Psi[, ok, drop = FALSE]); v <- sum(S); if (!is.finite(v) || v <= 0) next
    z[j] <- sum(b[ok]) / sqrt(v)                                               # fixed equal weights: only one variance is estimated
  }
  data.frame(feature = rownames(counts), p = 2 * stats::pt(-abs(z), n - ncol(X)))
}
register_candidate("erd_md", function(counts, meta, formula, tested_term) .erd_md(counts, meta, formula, tested_term),
  notes = "it7: ERD at depths {q0,q10,q25,q50}, each from eligible samples, GLS-combined with a joint HC3-type covariance")

# --- it8: the absolute estimand -- compositional centring by exposure-specific rarefaction depth ---
# erd_lm tests the relative scale: a x20 bloom of one common taxon dilutes every other taxon's
# relative abundance by a common factor c, moving E[f] for all of them. For f = 1 - (1-p)^D, a
# common factor on p is (to first order in p, exactly the regime where f is informative) a factor
# on depth: 1 - (1 - c p)^(D/c) ~ 1 - (1 - p)^D. So rarefying the exposed samples to D/c instead of D
# undoes the compositional shift INSIDE the transform, and the test stays an HC3 t-test on means,
# with its depth invariance. Generally: sample i is rarefied to D_i = D0 exp(-gamma x_i) (x the tested
# covariate, binary or continuous), D0 as large as every sample allows, and gamma is the value that
# puts the median feature's t-statistic at 0 -- the "most features are null" assumption, applied
# where it belongs. gamma is found on the whole cell (all features), before any feature is tested.
.erd_f_var <- function(counts, depth, Dvec) {                 # per-sample rarefaction depth
  Y <- as.matrix(counts); N <- matrix(depth, nrow(Y), ncol(Y), byrow = TRUE); Dm <- matrix(Dvec, nrow(Y), ncol(Y), byrow = TRUE)
  lr <- lchoose(N - Y, Dm) - lchoose(N, Dm); lr[!is.finite(lr)] <- -Inf; 1 - exp(lr)
}
.erd_c <- function(counts, meta, formula, tested_term, sub = 400L) {
  X <- stats::model.matrix(formula, meta); asg <- attr(X, "assign")
  tc <- which(asg == which(attr(stats::terms(formula), "term.labels") == tested_term))[1]
  depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts); x <- X[, tc]
  tstat <- function(gam, rows) {
    D0 <- min(depth * exp(gam * x)); Dv <- pmax(floor(D0 * exp(-gam * x)), 1)
    r <- .lm_hc3_fit(.erd_f_var(counts[rows, , drop = FALSE], depth, Dv), X, tc, .cluster_of(meta)); list(t = r$b / r$se, r = r) }
  set.seed(7L); rows <- if (nrow(counts) > sub) sort(sample.int(nrow(counts), sub)) else seq_len(nrow(counts))
  med <- function(g) stats::median(tstat(g, rows)$t, na.rm = TRUE)
  lo <- -2; hi <- 2; flo <- med(lo); fhi <- med(hi)
  gam <- if (is.finite(flo) && is.finite(fhi) && sign(flo) != sign(fhi))
    stats::uniroot(function(g) med(g), c(lo, hi), f.lower = flo, f.upper = fhi, tol = 1e-3)$root else 0
  fin <- tstat(gam, seq_len(nrow(counts)))
  data.frame(feature = rownames(counts), p = 2 * stats::pt(-abs(fin$t), fin$r$df), estimate = fin$r$b, gamma = gam)
}
register_candidate("erd_c", function(counts, meta, formula, tested_term) .erd_c(counts, meta, formula, tested_term)[, 1:3],
  notes = "it8: ERD with compositional centring -- exposed samples rarefied to D*exp(-gamma x), gamma setting the median t to 0 (absolute estimand)")

# --- it9: depth-stratified ERD -- rarefy within strata, compare within strata ---------------------
# ERD at the smallest library pays for invariance even when nothing is confounded (B_prev 25 vs 38
# TP for logistic_depth): every deep sample is thinned to the shallowest one. Invariance only needs
# samples COMPARED with each other to share a D. So: cut samples into K depth strata, rarefy each
# stratum to its own smallest library D_s, and estimate the exposure effect within strata (stratum
# fixed effects in the LM, HC3). Under H0, E[f] is equal across exposure levels inside every stratum,
# so the test stays valid for any K. K is chosen from the design alone: the information for the
# exposure contrast at depth D scales ~ D x (within-stratum sum of squares of x), so
# K = argmax_K sum_s D_s * SS_s(x), K in 1..6 (K = 1 is erd_lm). Unconfounded, strata hold both
# groups and D_s rises; confounded, strata go group-pure, SS_s -> 0, and K = 1 wins.
# Compositional centring (it8) carries over: D_i = D_s(i) exp(-gamma x_i).
.erd_strata <- function(depth, x, Kmax = 6L, X0 = NULL, hmax = 0.15) {
  best <- list(H = -Inf, s = rep(1L, length(depth)), K = 1L)
  for (K in seq_len(Kmax)) {
    br <- unique(stats::quantile(depth, seq(0, 1, length.out = K + 1L), type = 1)); if (length(br) < K + 1L) next
    s <- cut(depth, br, include.lowest = TRUE, labels = FALSE)
    if (any(tabulate(s) < 8L)) next
    # every stratum must keep at least half the design's exposure variance per sample: strata that
    # are nearly exposure-pure (depth confounding) carry a tiny, leverage-dominated contrast (smoke:
    # R19 picked K = 3 on D_s alone and lost 11 of 14 TP with FPR 0.069)
    if (K > 1L && any(vapply(split(x, s), function(v) mean((v - mean(v))^2), numeric(1)) < 0.5 * mean((x - mean(x))^2))) next
    if (!is.null(X0) && K > 1L) {             # HC3 needs every sample to have modest leverage: a stratum holding a
      Xs <- cbind(X0, stats::model.matrix(~ factor(s))[, -1, drop = FALSE])   # handful of one group would not
      if (qr(Xs)$rank < ncol(Xs)) next
      if (max(rowSums((Xs %*% solve(crossprod(Xs))) * Xs)) > hmax) next }
    H <- sum(vapply(split(seq_along(depth), s), function(i) min(depth[i]) * sum((x[i] - mean(x[i]))^2), numeric(1)))
    if (H > best$H) best <- list(H = H, s = s, K = K)
  }
  best
}
.erd_s <- function(counts, meta, formula, tested_term, centre = FALSE, sub = 400L) {
  X0 <- stats::model.matrix(formula, meta); asg <- attr(X0, "assign")
  tc0 <- which(asg == which(attr(stats::terms(formula), "term.labels") == tested_term))[1]
  depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts); x <- X0[, tc0]
  st <- .erd_strata(depth, x, X0 = X0); s <- st$s
  X <- if (max(s) > 1L) cbind(X0, stats::model.matrix(~ factor(s))[, -1, drop = FALSE]) else X0; tc <- tc0
  if (qr(X)$rank < ncol(X)) { X <- X0; s[] <- 1L }
  Dvec <- function(gam) { d <- numeric(length(depth))
    for (k in unique(s)) { i <- s == k; D0 <- min(depth[i] * exp(gam * x[i])); d[i] <- pmax(floor(D0 * exp(-gam * x[i])), 1) }; d }
  tstat <- function(gam, rows) { r <- .lm_hc3_fit(.erd_f_var(counts[rows, , drop = FALSE], depth, Dvec(gam)), X, tc, .cluster_of(meta)); list(t = r$b / r$se, r = r) }
  gam <- 0
  if (centre) {
    set.seed(7L); rows <- if (nrow(counts) > sub) sort(sample.int(nrow(counts), sub)) else seq_len(nrow(counts))
    med <- function(g) stats::median(tstat(g, rows)$t, na.rm = TRUE); flo <- med(-2); fhi <- med(2)
    if (is.finite(flo) && is.finite(fhi) && sign(flo) != sign(fhi))
      gam <- stats::uniroot(med, c(-2, 2), f.lower = flo, f.upper = fhi, tol = 1e-3)$root
  }
  fin <- tstat(gam, seq_len(nrow(counts)))
  data.frame(feature = rownames(counts), p = 2 * stats::pt(-abs(fin$t), fin$r$df), estimate = fin$r$b)
}
register_candidate("erd_s", function(counts, meta, formula, tested_term) .erd_s(counts, meta, formula, tested_term, FALSE),
  notes = "it9: depth-stratified ERD (K strata from the design, each rarefied to its own minimum), stratum fixed effects, HC3")
register_candidate("erd_sc", function(counts, meta, formula, tested_term) .erd_s(counts, meta, formula, tested_term, TRUE),
  notes = "it9: erd_s + compositional centring by exposure-specific depth (absolute estimand)")

# --- it11: one test across the transform family -- detection AND log, depth-centred, max-combined ---
# srv2 (msq/mid): the family holds -- erd_c worst FPR 0.055, erl_lm 0.058 -- but the two transforms
# win on different simulators. house/implant signal is mostly prevalence (ERD >> ERL: R16 8.7 vs
# 2.3); msq/mid signal is abundance (ERL >> ERD: msq R00 10.0 vs 3.4). Which one carries the signal
# is not knowable from the pooled data, so test both and pay only the multiplicity of their
# correlation: max(|t_det|, |t_log|) referred to the bivariate normal with the correlation estimated
# from the two coefficients' HC3 influence functions (not ACAT: no dilution by an uninformative arm
# beyond the max's own price, which shrinks as the two correlate). A 2-df Wald on the pair is the
# alternative. The it8 depth-scaling centring holds for any h(Y_D): Y_D ~ Bin(D, p) ~ Pois(Dp) depends
# on D and p only through Dp, so rarefying exposed samples to D/c undoes a compositional c for the
# log transform as well.
.erl_f_var <- function(counts, depth, Dvec, kmax = 70L, big_mu = 20) {
  # exact hypergeometric sum where D p <= big_mu; above it a THIRD-order moment expansion, which is
  # exact at D = N (no thinning: variance and skew vanish) and within ~1e-4 elsewhere, so a thinned
  # and an unthinned sample are treated alike (the invariance the U-statistic relies on)
  Y <- as.matrix(counts); N <- matrix(depth, nrow(Y), ncol(Y), byrow = TRUE); Dm <- matrix(Dvec, nrow(Y), ncol(Y), byrow = TRUE)
  mu <- Dm * Y / N; out <- matrix(0, nrow(Y), ncol(Y))
  big <- mu > big_mu; if (any(big)) { pp <- (Y / N)[big]; Nb <- N[big]; Db <- Dm[big]; fpc <- (Nb - Db) / pmax(Nb - 1, 1)
    v <- Db * pp * (1 - pp) * fpc; k3 <- v * (1 - 2 * pp) * (Nb - 2 * Db) / pmax(Nb - 2, 1); m1 <- 1 + mu[big]
    out[big] <- log(m1) - v / (2 * m1^2) + k3 / (3 * m1^3) }
  sm <- which(!big & Y > 0)
  if (length(sm)) { y <- Y[sm]; nn <- N[sm] - Y[sm]; dd <- Dm[sm]; acc <- numeric(length(sm))
    for (k in 1:kmax) { live <- k <= y & k <= dd; if (!any(live)) break
      acc[live] <- acc[live] + stats::dhyper(k, y[live], nn[live], dd[live]) * log1p(k) }
    out[sm] <- acc }
  out
}
# HC3 influence functions of the tested coefficient, one row per feature (features x samples)
.lm_infl <- function(Fm, X, tc, cl = NULL) {
  XtXi <- solve(crossprod(X)); a <- drop(XtXi[tc, , drop = FALSE] %*% t(X)); H <- rowSums((X %*% XtXi) * X)
  B <- Fm %*% t(XtXi %*% t(X)); R <- Fm - B %*% t(X)
  Psi <- R * rep(a, each = nrow(R)); if (is.null(cl)) Psi <- Psi / rep(1 - pmin(H, 0.99), each = nrow(R))
  list(b = B[, tc], Psi = Psi, cl = cl)
}
.infl_cov <- function(P1, P2, cl) { if (is.null(cl)) return(rowSums(P1 * P2))
  Z <- stats::model.matrix(~ cl - 1); G <- ncol(Z); rowSums((P1 %*% Z) * (P2 %*% Z)) * G / (G - 1) }
.gauss_legendre <- function(n) {                                  # Golub-Welsch; no package dependency
  k <- seq_len(n - 1L); b <- k / sqrt(4 * k^2 - 1); J <- matrix(0, n, n); J[cbind(k, k + 1L)] <- b; J[cbind(k + 1L, k)] <- b
  e <- eigen(J, symmetric = TRUE); o <- order(e$values); list(nodes = e$values[o], weights = 2 * e$vectors[1, o]^2)
}
.pmax2 <- function(m, rho) {                                      # P(max(|Z1|,|Z2|) >= m), corr rho, vectorised
  gl <- .gauss_legendre(48L); out <- numeric(length(m))
  for (i in seq_along(m)) { if (!is.finite(m[i]) || !is.finite(rho[i])) { out[i] <- NA; next }
    r <- max(min(rho[i], 0.999), -0.999); s <- sqrt(1 - r^2); z <- m[i] * gl$nodes
    inner <- stats::pnorm((m[i] - r * z) / s) - stats::pnorm((-m[i] - r * z) / s)
    out[i] <- 1 - m[i] * sum(gl$weights * stats::dnorm(z) * inner) }
  pmin(pmax(out, 0), 1)
}
.erdl_memo <- new.env()
.erdl_fit <- function(counts, meta, formula, tested_term, sub = 400L) {
  key <- list(counts, meta, formula, tested_term)
  if (!is.null(.erdl_memo$key) && identical(.erdl_memo$key, key)) return(.erdl_memo$val)
  X <- stats::model.matrix(formula, meta); asg <- attr(X, "assign")
  tc <- which(asg == which(attr(stats::terms(formula), "term.labels") == tested_term))[1]
  depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts); x <- X[, tc]; cl <- .cluster_of(meta)
  Dv <- function(gam) { D0 <- min(depth * exp(gam * x)); pmax(floor(D0 * exp(-gam * x)), 1) }
  fmat <- function(h, gam, rows) if (h == "det") .erd_f_var(counts[rows, , drop = FALSE], depth, Dv(gam)) else .erl_f_var(counts[rows, , drop = FALSE], depth, Dv(gam))
  set.seed(7L); rows <- if (nrow(counts) > sub) sort(sample.int(nrow(counts), sub)) else seq_len(nrow(counts))
  gam_for <- function(h) { med <- function(g) { r <- .lm_hc3_fit(fmat(h, g, rows), X, tc, cl); stats::median(r$b / r$se, na.rm = TRUE) }
    flo <- med(-2); fhi <- med(2)
    if (is.finite(flo) && is.finite(fhi) && sign(flo) != sign(fhi)) stats::uniroot(med, c(-2, 2), f.lower = flo, f.upper = fhi, tol = 1e-3)$root else 0 }
  all <- seq_len(nrow(counts)); g1 <- gam_for("det"); g2 <- gam_for("log")
  F1 <- fmat("det", g1, all); F2 <- fmat("log", g2, all)
  I1 <- .lm_infl(F1, X, tc, cl); I2 <- .lm_infl(F2, X, tc, cl)
  v1 <- .infl_cov(I1$Psi, I1$Psi, cl); v2 <- .infl_cov(I2$Psi, I2$Psi, cl); c12 <- .infl_cov(I1$Psi, I2$Psi, cl)
  ok <- rowSums(counts > 0) >= 3 & v1 > 0 & v2 > 0
  # t -> z through the reference t (df = n - p, or G - 1 clusters) so the max/Wald respect small-df tails
  df <- if (is.null(cl)) ncol(counts) - ncol(X) else nlevels(cl) - 1L
  tz <- function(t) stats::qnorm(stats::pt(t, df))
  z1 <- ifelse(ok, tz(I1$b / sqrt(v1)), NA); z2 <- ifelse(ok, tz(I2$b / sqrt(v2)), NA); rho <- ifelse(ok, c12 / sqrt(v1 * v2), NA)
  pm <- .pmax2(pmax(abs(z1), abs(z2)), rho)
  w2 <- (z1^2 - 2 * rho * z1 * z2 + z2^2) / (1 - pmin(rho^2, 0.998)); pw <- stats::pchisq(w2, 2, lower.tail = FALSE)
  val <- list(feature = rownames(counts), p_max = pm, p_w2 = pw, p_log = 2 * stats::pnorm(-abs(z2)), est_log = I2$b, gam = c(g1, g2), rho = rho)
  .erdl_memo$key <- key; .erdl_memo$val <- val; val
}
register_candidate("erdl_max", function(counts, meta, formula, tested_term) {
  v <- .erdl_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_max, estimate = v$est_log)
}, notes = "it11: depth-centred ERD and ERL, max(|t|) against their bivariate normal (HC3 / cluster-robust influence correlation)")
register_candidate("erdl_w2", function(counts, meta, formula, tested_term) {
  v <- .erdl_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_w2, estimate = v$est_log)
}, notes = "it11: depth-centred ERD and ERL, joint 2-df Wald")
register_candidate("erl_c", function(counts, meta, formula, tested_term) {
  v <- .erdl_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_log, estimate = v$est_log)
}, notes = "it11: expected rarefied log count with depth-scaling compositional centring (it8 trick applied to ERL)")

# --- it12: pairwise common-depth rarefaction -- a two-sample U-statistic -----------------------------
# Outer loop (p03): the family is the best-calibrated in the benchmark but gives up power without
# confounding (B_prev 4.2 vs 8.2 TP for logistic). Diagnosis on dev cells: a linear-probability test
# on RAW detection matches logistic (B_prev twinsuk 74 vs 74), ERD at the smallest library gets 57-65
# -- the loss is the thinning, not the statistic. Invariance does not need one global D: it needs the
# two samples being COMPARED to share a depth. So compare every exposed sample i with every control j
# at their own common depth D_ij = min(N_i c, N_j) (c = exp(gamma), the it8 compositional factor):
#     h_ij = f_i(D_ij / c) - f_j(D_ij),   U = mean_ij h_ij.
# Under H0, E h_ij = 0 for every pair whatever the depths, so E U = 0 exactly. Unconfounded, most pairs
# have similar depths and are barely thinned; confounded, a deep/shallow pair is thinned to the shallow
# one -- ERD's behaviour, pair by pair instead of globally. Variance: two-sample U-statistic (the
# DeLong/Hoeffding form) -- var U = var_i(a_i)/n1 + var_j(b_j)/n0 with a_i = mean_j h_ij,
# b_j = mean_i h_ij -- Welch-Satterthwaite df. Binary exposure without other terms; otherwise it falls
# back to erd_c (LM). gamma comes from erd_c's own search (same compositional factor).
# Large n: each exposed sample is compared with at most K controls, spread evenly across the control
# list (an incomplete U-statistic; cost n1*K*m instead of n1*n0*m -- 128 s -> ~30 s at n = 400).
# The DeLong estimate from partner means then carries the partner-sampling term automatically, so
# it errs conservative. With n0 <= K (the dev suite's n = 100) nothing changes.
# rho < 1 (it13): both members of every pair are thinned, to rho x their common depth. srv4: with
# rho = 1 the shallower sample of a pair is compared UNTHINNED against a thinned deep one, which is
# equivalent only if the simulator makes low-depth zeros by read thinning; MIDASim decides presence
# from library size its own way, and the detection U-statistic went to FPR 0.087 / 0.094 at mid
# R17 / R19 (the global-minimum ERD, which thins nearly every sample, held at 0.050 / 0.054).
.u_stat <- function(counts, depth, g1, gam, h = c("det", "log"), K = 50L, rho = 1) {
  h <- match.arg(h); Y <- as.matrix(counts); storage.mode(Y) <- "double"; c <- exp(gam)
  i1 <- which(g1); i0 <- which(!g1); n1 <- length(i1); n0 <- length(i0); m <- nrow(Y); Kk <- min(K, n0)
  step <- max(1L, n0 %/% Kk)
  A <- matrix(0, m, n1); Bs <- matrix(0, m, n0); cnt <- numeric(n0)
  fcol <- function(Ysub, Nvec, Dvec) if (h == "det") .erd_f_var(Ysub, Nvec, Dvec) else .erl_f_var(Ysub, Nvec, Dvec)
  for (k in seq_len(n1)) { i <- i1[k]
    J <- if (Kk == n0) seq_len(n0) else ((k - 1L + (seq_len(Kk) - 1L) * step) %% n0) + 1L
    jj <- i0[J]
    D <- pmax(floor(rho * pmin(depth[i] * c, depth[jj])), 1)               # (rho x) common depth with each partner control
    f0 <- fcol(Y[, jj, drop = FALSE], depth[jj], D)                        # controls at D
    f1 <- fcol(Y[, rep(i, length(jj)), drop = FALSE], rep(depth[i], length(jj)), pmax(floor(D / c), 1))   # exposed at D/c
    H <- f1 - f0; A[, k] <- rowMeans(H); Bs[, J] <- Bs[, J] + H; cnt[J] <- cnt[J] + 1 }
  keepb <- cnt > 0; B <- sweep(Bs[, keepb, drop = FALSE], 2, cnt[keepb], "/"); n0 <- sum(keepb); U <- rowMeans(A)
  va <- apply(A, 1, stats::var) / n1; vb <- apply(B, 1, stats::var) / n0; v <- va + vb
  df <- v^2 / (va^2 / (n1 - 1) + vb^2 / (n0 - 1))
  list(U = U, v = v, df = df, a = A - U, b = B - U)                        # centred influence pieces
}
.u_design <- function(meta, formula, tested_term) {
  tl <- attr(stats::terms(formula), "term.labels"); g <- meta[[tested_term]]
  if (length(tl) != 1L || length(unique(g)) != 2L || !is.null(.cluster_of(meta))) return(NULL)
  g == sort(unique(g))[2]
}
register_candidate("erd_u", function(counts, meta, formula, tested_term) {
  g1 <- .u_design(meta, formula, tested_term)
  if (is.null(g1)) { r <- .erd_c(counts, meta, formula, tested_term); return(r[, 1:3]) }
  depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts)
  gam <- .erd_c(counts, meta, formula, tested_term)$gamma[1]
  u <- .u_stat(counts, depth, g1, gam, "det"); ok <- rowSums(counts > 0) >= 3 & u$v > 0
  data.frame(feature = rownames(counts), p = ifelse(ok, 2 * stats::pt(-abs(u$U / sqrt(u$v)), u$df), NA), estimate = u$U)
}, notes = "it12: pairwise common-depth expected detection, two-sample U-statistic (DeLong variance), compositional centring from erd_c")

# gamma estimated ON the U-statistic (median z = 0 over a feature subsample): with pairs thinned to
# their common depth rather than the global minimum, D_ij p is larger and the first-order
# depth-scaling identity is less exact, so a gamma borrowed from erd_c left a residual compositional
# bias (hard bloom, tongue: null FPR 0.058-0.071 vs erd_c 0.041-0.049).
.u_gamma <- function(counts, depth, g1, h = "det", sub = 300L, rho = 1) {
  set.seed(7L); rows <- if (nrow(counts) > sub) sort(sample.int(nrow(counts), sub)) else seq_len(nrow(counts))
  ct <- counts[rows, , drop = FALSE]; keep <- rowSums(ct > 0) >= 3; ct <- ct[keep, , drop = FALSE]
  med <- function(g) { u <- .u_stat(ct, depth, g1, g, h, rho = rho); stats::median(u$U / sqrt(u$v), na.rm = TRUE) }
  flo <- med(-2); fhi <- med(2)
  if (is.finite(flo) && is.finite(fhi) && sign(flo) != sign(fhi)) stats::uniroot(med, c(-2, 2), f.lower = flo, f.upper = fhi, tol = 2e-3)$root else 0
}
.eu_memo <- new.env()
.eu_fit <- function(counts, meta, formula, tested_term, rho = 1) {
  key <- list(counts, meta, formula, tested_term, rho); slot <- paste0("r", rho)
  if (!is.null(.eu_memo[[slot]]) && identical(.eu_memo[[slot]]$key, key)) return(.eu_memo[[slot]]$val)
  g1 <- .u_design(meta, formula, tested_term)
  if (is.null(g1)) { v <- .erdl_fit(counts, meta, formula, tested_term)             # covariates / continuous / clusters
    val <- list(feature = v$feature, p_det = NA, p_log = v$p_log, p_max = v$p_max, est = v$est_log, fallback = TRUE)
    val$p_det <- .erd_c(counts, meta, formula, tested_term)$p
  } else {
    depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts)
    # one compositional factor for both transforms, estimated on the LOG U-statistic: under depth
    # confounding the detection U-statistic's median-z curve is nearly flat in gamma and its root
    # wandered (house R19 twinsuk r1: -1.21, which put the log test at null FPR 0.251); the log
    # scale moves by ~log c for every common taxon and pins c down (0.07 there, FPR 0.064).
    gl <- .u_gamma(counts, depth, g1, "log", rho = rho); gd <- gl
    ud <- .u_stat(counts, depth, g1, gd, "det", rho = rho); ul <- .u_stat(counts, depth, g1, gl, "log", rho = rho)
    ok <- rowSums(counts > 0) >= 3 & ud$v > 0 & ul$v > 0
    n1 <- sum(g1); n0 <- sum(!g1)
    cv <- (rowSums(ud$a * ul$a) / (n1 - 1)) / n1 + (rowSums(ud$b * ul$b) / (n0 - 1)) / n0   # joint U-stat covariance
    tz <- function(t, df) stats::qnorm(stats::pt(t, df))
    z1 <- ifelse(ok, tz(ud$U / sqrt(ud$v), ud$df), NA); z2 <- ifelse(ok, tz(ul$U / sqrt(ul$v), ul$df), NA)
    rho <- ifelse(ok, cv / sqrt(ud$v * ul$v), NA)
    val <- list(feature = rownames(counts), p_det = 2 * stats::pnorm(-abs(z1)), p_log = 2 * stats::pnorm(-abs(z2)),
                p_max = .pmax2(pmax(abs(z1), abs(z2)), rho), est = ul$U / log(2), gam = c(gd, gl), fallback = FALSE)
  }
  .eu_memo[[slot]] <- list(key = key, val = val); val
}
register_candidate("erd_u2", function(counts, meta, formula, tested_term) {
  v <- .eu_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_det)
}, notes = "it12: pairwise common-depth detection U-statistic, gamma estimated on the log U-statistic")
register_candidate("erl_u", function(counts, meta, formula, tested_term) {
  v <- .eu_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_log, estimate = v$est)
}, notes = "it12: pairwise common-depth expected log count U-statistic")
register_candidate("erdl_u", function(counts, meta, formula, tested_term) {
  v <- .eu_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_max, estimate = v$est)
}, notes = "it12: pairwise common-depth detection AND log U-statistics, max(|z|) with their joint U-statistic covariance")

# it13: both pair members thinned (rho < 1)
for (.rho in c(0.5, 0.7)) local({ r <- .rho; tag <- sprintf("r%02d", round(10 * r))
  register_candidate(paste0("erdl_u_", tag), function(counts, meta, formula, tested_term) {
    v <- .eu_fit(counts, meta, formula, tested_term, rho = r); data.frame(feature = v$feature, p = v$p_max, estimate = v$est)
  }, notes = sprintf("it13: erdl_u with both pair members thinned to %.1f x their common depth", r))
  register_candidate(paste0("erd_u_", tag), function(counts, meta, formula, tested_term) {
    v <- .eu_fit(counts, meta, formula, tested_term, rho = r); data.frame(feature = v$feature, p = v$p_det)
  }, notes = sprintf("it13: detection U-statistic, both pair members thinned to %.1f x their common depth", r))
})
