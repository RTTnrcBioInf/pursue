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
.erl_f_var <- function(counts, depth, Dvec, kmax = 70L, big_mu = 20, off = 1) {   # E[log(off + Y_D)]; off = 0: E[log Y_D; Y_D > 0]
  # exact hypergeometric sum where D p <= big_mu; above it a THIRD-order moment expansion, which is
  # exact at D = N (no thinning: variance and skew vanish) and within ~1e-4 elsewhere, so a thinned
  # and an unthinned sample are treated alike (the invariance the U-statistic relies on)
  Y <- as.matrix(counts); N <- matrix(depth, nrow(Y), ncol(Y), byrow = TRUE); Dm <- matrix(Dvec, nrow(Y), ncol(Y), byrow = TRUE)
  mu <- Dm * Y / N; out <- matrix(0, nrow(Y), ncol(Y))
  big <- mu > big_mu; if (any(big)) { pp <- (Y / N)[big]; Nb <- N[big]; Db <- Dm[big]; fpc <- (Nb - Db) / pmax(Nb - 1, 1)
    v <- Db * pp * (1 - pp) * fpc; k3 <- v * (1 - 2 * pp) * (Nb - 2 * Db) / pmax(Nb - 2, 1); m1 <- off + mu[big]
    out[big] <- log(m1) - v / (2 * m1^2) + k3 / (3 * m1^3)
    if (off == 0) out[big] <- out[big] - (v * (1 - 6 * pp * (1 - pp)) + 3 * v^2) / (4 * m1^4) }   # 4th order: log y has no +1 cushion
  sm <- which(!big & Y > 0)
  if (length(sm)) { y <- Y[sm]; nn <- N[sm] - Y[sm]; dd <- Dm[sm]; acc <- numeric(length(sm))
    for (k in 1:kmax) { live <- k <= y & k <= dd; if (!any(live)) break
      acc[live] <- acc[live] + stats::dhyper(k, y[live], nn[live], dd[live]) * log(off + k) }
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

# --- it14: thinning factor from the design ---------------------------------------------------------
# srv4/srv5: rho = 1 is anti-conservative on mid ONLY under depth confounding (R17/R19); at mid
# R00/R06/R11 it is calibrated (0.047-0.053), because when depth is unrelated to the exposure, which
# member of a pair is the shallower (unthinned) one is random with respect to group and the
# thinning-inconsistency bias cancels. So rho is set from the design alone -- labels and depths,
# never counts, hence no effect on validity under H0: d = |mean log N (exposed) - mean log N
# (controls)| / pooled sd of log N; rho = 1 for d <= 0.25 (randomisation noise at n = 50 per group
# has sd ~0.2), 0.5 for d >= 0.5, linear between. Unconfounded cells keep rho = 1's power (B_prev
# twinsuk 65 vs 55 TP at rho 0.5); confounded cells get srv5's protection.
.rho_design <- function(depth, g1, lo = 0.25, hi = 0.5, rmin = 0.5) {
  l <- log(depth); s <- sqrt(((sum(g1) - 1) * stats::var(l[g1]) + (sum(!g1) - 1) * stats::var(l[!g1])) / (length(l) - 2))
  d <- abs(mean(l[g1]) - mean(l[!g1])) / max(s, 1e-8)
  1 - (1 - rmin) * min(max((d - lo) / (hi - lo), 0), 1)
}
.eu_fit_ad <- function(counts, meta, formula, tested_term) {
  g1 <- .u_design(meta, formula, tested_term)
  r <- if (is.null(g1)) 1 else .rho_design(if (!is.null(meta$depth)) meta$depth else colSums(counts), g1)
  r <- round(r, 2); v <- .eu_fit(counts, meta, formula, tested_term, rho = r); v$rho <- r; v
}
register_candidate("erdl_u_ad", function(counts, meta, formula, tested_term) {
  v <- .eu_fit_ad(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_max, estimate = v$est)
}, notes = "it14: pairwise U-statistics, detection + log max, thinning factor rho chosen from the depth-exposure imbalance (1 -> 0.5)")
register_candidate("erd_u_ad", function(counts, meta, formula, tested_term) {
  v <- .eu_fit_ad(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_det)
}, notes = "it14: detection U-statistic with design-chosen rho")

# --- it15: covariates inside the pairwise test -- pairwise-difference regression ---------------------
# p06: erdlu_ad exceeds FDR 0.10 only at R21 (a confounder correlated with the exposure, phi 0.7,
# that also moves 10% of null taxa), where it falls back to the global-minimum LM. Keeping the
# pairwise thinning: for each case i / control j pair, h_ij = f_i(D_ij/c) - f_j(D_ij) is regressed on
# x_ij = (1, z_i - z_j) over all pairs: h_ij = theta + beta'(z_i - z_j) + e_ij. theta is the
# covariate-adjusted exposure effect (with no covariates it is exactly the U-statistic). The
# estimating function sum_ij x_ij e_ij is a two-sample U-statistic, so its variance is the same
# Hoeffding projection as before -- per-sample partner means a_i, b_j of x e -- and
# Var(gamma) = M^-1 (Var a / n1 + Var b / n0) M^-1 with M = mean x x'. Accumulated per sample
# (sums of x h and x x'), so no pair-level storage.
# vcorr (it16): the DeLong/projection estimate var(a)/n1 + var(b)/n0 counts the pair-level residual
# variance twice (each partner mean carries sigma_e^2 / #partners), while the true variance carries
# it once: Var(theta) ~ s10/n1 + s01/n0 + sigma_e^2 [M^-1]_11 / P. It is second-order in general,
# but here h_ij depends on the PAIR's common depth, so on the log scale sigma_e^2 is large and the
# test runs conservative in the tail (null FPR 0.0003 at 0.001 for the log half). vcorr subtracts
# one copy, sigma_e^2 estimated from the two-way decomposition of the pair residuals e_ij.
# raw (it17): no thinning at all -- each sample at its own depth -- with the pair's log-depth
# difference as a covariate, offset by the compositional factor for the exposed sample
# ((log N_i + gamma) - log N_j). Used only when the design is depth-balanced (the choice uses labels
# and depths only), where it is logistic_depth's information inside the pairwise machinery.
# cl (it18): subject labels for repeated measures. The pairwise estimator is unchanged; its
# projection is summed within subject before the variance (cluster-robust, G/(G-1)), df from the
# numbers of exposed and control subjects.
# wk (it21): pair weights w_ij = D_ij^wk (design-only, so any value keeps the test exact); the pairwise
# regression becomes weighted least squares and each sample's projection is scaled by its share of the
# total weight (W_i / W). wk = 0 is the unweighted estimator exactly.
.u_stat_cov <- function(counts, depth, g1, Z, gam, h = c("det", "log", "pi", "cen", "pre", "sqrt", "pos", "q4"), K = 50L, rho = 1, vcorr = FALSE, raw = FALSE, cl = NULL, wk = 0, tob = NULL, wfun = NULL, dslope = NULL, S = NULL, coef = 1L, sf = NULL) {
  h <- match.arg(h); Y <- as.matrix(counts); storage.mode(Y) <- "double"; c <- exp(gam)
  if (raw) { Fall <- if (h == "det") (Y > 0) * 1 else log1p(Y); ld <- log(depth)
    Z <- cbind(Z, rep(0, length(depth))) }                                    # placeholder column, filled per pair
  i1 <- which(g1); i0 <- which(!g1); n1 <- length(i1); n0 <- length(i0); m <- nrow(Y); Kk <- min(K, n0)
  step <- max(1L, n0 %/% Kk); q <- if (is.null(Z)) 0L else ncol(Z); p <- 1L + q
  Axh <- array(0, c(m, p, n1)); Axx <- array(0, c(p, p, n1)); Kc <- numeric(n1); SSh <- numeric(m)
  Bxh <- array(0, c(m, p, n0)); Bxx <- array(0, c(p, p, n0)); cnt <- numeric(n0)
  fcol <- function(Ysub, Nvec, Dvec) if (h == "det") .erd_f_var(Ysub, Nvec, Dvec) else if (h == "sqrt") .ers_f_var(Ysub, Nvec, Dvec) else if (h == "q4") .erp_f_var(Ysub, Nvec, Dvec, 0.25) else .erl_f_var(Ysub, Nvec, Dvec)
  for (k in seq_len(n1)) { i <- i1[k]
    J <- if (Kk == n0) seq_len(n0) else ((k - 1L + (seq_len(Kk) - 1L) * step) %% n0) + 1L; jj <- i0[J]; nk <- length(jj)
    if (raw) {
      H <- Fall[, i] - Fall[, jj, drop = FALSE]
      X <- cbind(1, -sweep(Z[jj, , drop = FALSE], 2, Z[i, ], "-")); X[, p] <- (ld[i] + gam) - ld[jj]
    } else {
    D <- pmax(floor(rho * pmin(depth[i] * c, depth[jj])), 1)
    H <- if (h == "pre") S[, i] - S[, jj, drop = FALSE] else if (h == "pi") .pi_kernel(Y[, i], depth[i], pmax(floor(D / c), 1), Y[, jj, drop = FALSE], depth[jj], D) else
      if (h == "pos") .pos_kernel(Y[, rep(i, nk), drop = FALSE], rep(depth[i], nk), pmax(floor(D / c), 1), Y[, jj, drop = FALSE], depth[jj], D) else
      if (h == "cen") .cen_kernel(Y[, rep(i, nk), drop = FALSE], rep(depth[i], nk), pmax(floor(D / c), 1), Y[, jj, drop = FALSE], depth[jj], D, tob) else
      if (!is.null(sf)) { bi <- c * exp(sf[i]); bj <- exp(sf[jj]); De <- rho * pmin(depth[i] * bi, depth[jj] * bj)   # it34: per-sample size factors
        fcol(Y[, rep(i, nk), drop = FALSE], rep(depth[i], nk), pmax(floor(De / bi), 1)) - fcol(Y[, jj, drop = FALSE], depth[jj], pmax(floor(De / bj), 1)) } else
      fcol(Y[, rep(i, nk), drop = FALSE], rep(depth[i], nk), pmax(floor(D / c), 1)) - fcol(Y[, jj, drop = FALSE], depth[jj], D)
    X <- if (q) cbind(1, -sweep(Z[jj, , drop = FALSE], 2, Z[i, ], "-")) else matrix(1, nk, 1)   # (1, z_i - z_j)
    if (!is.null(dslope)) H <- H - outer(dslope, log(depth[i]) - log(depth[jj]))   # it25: pooled depth slope removed
    }
    w <- if (!is.null(wfun)) wfun(i, jj) else if (wk == 0 || raw) rep(1, nk) else { Dw <- pmax(floor(rho * pmin(depth[i] * c, depth[jj])), 1); (Dw / stats::median(depth))^wk }
    Hw <- H * rep(w, each = m)
    Axh[, , k] <- Hw %*% X; Axx[, , k] <- crossprod(X * sqrt(w)); Kc[k] <- sum(w); SSh <- SSh + rowSums(H^2)
    for (cc in seq_len(p)) Bxh[, cc, J] <- Bxh[, cc, J] + Hw * rep(X[, cc], each = m)
    for (t in seq_len(nk)) Bxx[, , J[t]] <- Bxx[, , J[t]] + w[t] * tcrossprod(X[t, ])
    cnt[J] <- cnt[J] + w }
  P <- sum(Kc); M <- apply(Axx, c(1, 2), sum) / P; Mi <- solve(M)
  G <- (apply(Axh, c(1, 2), sum) / P) %*% Mi                                  # m x p: gamma per feature (rows)
  # partner-mean estimating functions per sample, mapped through M^-1; component 1 = theta's influence
  phiA <- vapply(seq_len(n1), function(k) ((Axh[, , k, drop = TRUE] - G %*% Axx[, , k]) / Kc[k]) %*% Mi[, coef], numeric(m))
  keep <- cnt > 0
  phiB <- vapply(which(keep), function(j) ((matrix(Bxh[, , j], m, p) - G %*% Bxx[, , j]) / cnt[j]) %*% Mi[, coef], numeric(m))
  phiA <- matrix(phiA, m); phiB <- matrix(phiB, m); n0k <- ncol(phiB)
  if ((wk != 0 || !is.null(wfun)) && !raw) {                                   # weight shares: W_i / W (reduces to 1/n1 unweighted)
    phiA <- phiA * rep(Kc / P * n1, each = m); phiB <- phiB * rep(cnt[keep] / P * n0k, each = m) }
  va <- apply(phiA, 1, stats::var) / n1; vb <- apply(phiB, 1, stats::var) / n0k; v <- va + vb
  df <- v^2 / (va^2 / (n1 - 1) + vb^2 / (n0k - 1))
  vind <- v; dfind <- df
  if (!is.null(cl)) {
    ca <- factor(cl[i1]); cb <- factor(cl[i0[keep]]); Ga <- nlevels(ca); Gb <- nlevels(cb)
    Sa <- (phiA - rowMeans(phiA)) %*% stats::model.matrix(~ ca - 1); Sb <- (phiB - rowMeans(phiB)) %*% stats::model.matrix(~ cb - 1)
    va <- rowSums(Sa^2) / n1^2 * Ga / (Ga - 1); vb <- rowSums(Sb^2) / n0k^2 * Gb / (Gb - 1); v <- va + vb
    df <- v^2 / (va^2 / (Ga - 1) + vb^2 / (Gb - 1))
    attr(v, "clusters") <- c(Ga, Gb)
  }
  if (vcorr) {
    Sxx <- apply(Axx, c(1, 2), sum); Sxh <- apply(Axh, c(1, 2), sum)
    SSe <- SSh - 2 * rowSums(G * Sxh) + rowSums((G %*% Sxx) * G)                 # sum of e_ij^2 over pairs
    ea <- vapply(seq_len(n1), function(k) (Axh[, 1, k] - G %*% Axx[, 1, k]) / Kc[k], numeric(m))   # partner means of e
    eb <- vapply(which(keep), function(j) (Bxh[, 1, j] - G %*% Bxx[, 1, j]) / cnt[j], numeric(m))
    ea <- matrix(ea, m); eb <- matrix(eb, m)
    SSr <- SSe - colSums(t(ea^2) * Kc) - colSums(t(eb^2) * cnt[keep])           # two-way residual SS
    s2e <- pmax(SSr, 0) / max(P - n1 - n0k + 1, 1)
    v <- pmax(v - s2e * Mi[1, 1] / P, 0.5 * v)
  }
  out <- list(U = G[, 1], v = v, df = df, a = phiA - rowMeans(phiA), b = phiB - rowMeans(phiB), vind = vind, dfind = dfind, B = G, M = M)
  if (!is.null(cl)) { out$ca <- factor(cl[i1]); out$cb <- factor(cl[i0[keep]]) }
  out
}
.uc_design <- function(meta, formula, tested_term) {
  tl <- attr(stats::terms(formula), "term.labels"); g <- meta[[tested_term]]
  if (length(unique(g)) != 2L || !is.null(.cluster_of(meta))) return(NULL)
  X <- stats::model.matrix(formula, meta); asg <- attr(X, "assign")
  zc <- which(asg != 0 & asg != which(tl == tested_term))
  list(g1 = g == sort(unique(g))[2], Z = if (length(zc)) X[, zc, drop = FALSE] else NULL)
}
.euc_memo <- new.env()
.euc_fit <- function(counts, meta, formula, tested_term) {
  key <- list(counts, meta, formula, tested_term)
  if (!is.null(.euc_memo$key) && identical(.euc_memo$key, key)) return(.euc_memo$val)
  d <- .uc_design(meta, formula, tested_term)
  if (is.null(d)) { val <- .eu_fit_ad(counts, meta, formula, tested_term); val$fallback <- TRUE
  } else {
    depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts); g1 <- d$g1; Z <- d$Z
    rho <- round(.rho_design(depth, g1), 2)
    set.seed(7L); rows <- if (nrow(counts) > 300L) sort(sample.int(nrow(counts), 300L)) else seq_len(nrow(counts))
    ct <- counts[rows, , drop = FALSE]; ct <- ct[rowSums(ct > 0) >= 3, , drop = FALSE]
    med <- function(gm) { u <- .u_stat_cov(ct, depth, g1, Z, gm, "log", rho = rho); stats::median(u$U / sqrt(u$v), na.rm = TRUE) }
    flo <- med(-2); fhi <- med(2)
    gm <- if (is.finite(flo) && is.finite(fhi) && sign(flo) != sign(fhi)) stats::uniroot(med, c(-2, 2), f.lower = flo, f.upper = fhi, tol = 2e-3)$root else 0
    ud <- .u_stat_cov(counts, depth, g1, Z, gm, "det", rho = rho); ul <- .u_stat_cov(counts, depth, g1, Z, gm, "log", rho = rho)
    ok <- rowSums(counts > 0) >= 3 & ud$v > 0 & ul$v > 0; n1 <- sum(g1); n0 <- ncol(ud$b)
    cv <- (rowSums(ud$a * ul$a) / (n1 - 1)) / n1 + (rowSums(ud$b * ul$b) / (n0 - 1)) / n0
    tz <- function(t, df) stats::qnorm(stats::pt(t, df))
    z1 <- ifelse(ok, tz(ud$U / sqrt(ud$v), ud$df), NA); z2 <- ifelse(ok, tz(ul$U / sqrt(ul$v), ul$df), NA)
    rho_c <- ifelse(ok, cv / sqrt(ud$v * ul$v), NA)
    val <- list(feature = rownames(counts), p_det = 2 * stats::pnorm(-abs(z1)), p_log = 2 * stats::pnorm(-abs(z2)),
                p_max = .pmax2(pmax(abs(z1), abs(z2)), rho_c), est = ul$U / log(2), gam = gm, rho = rho, fallback = FALSE)
    attr(val, "rho_c") <- rho_c
  }
  .euc_memo$key <- key; .euc_memo$val <- val; val
}
register_candidate("erdl_uc", function(counts, meta, formula, tested_term) {
  v <- .euc_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_max, estimate = v$est)
}, notes = "it15: erdl_u_ad with covariates handled by pairwise-difference regression (binary exposure)")
register_candidate("erd_uc", function(counts, meta, formula, tested_term) {
  v <- .euc_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_det)
}, notes = "it15: detection-only version of erdl_uc")

# --- it16: unbiased pairwise variance (vcorr) ---
.euc2_memo <- new.env()
.euc_fit2 <- function(counts, meta, formula, tested_term) {
  key <- list(counts, meta, formula, tested_term)
  if (!is.null(.euc2_memo$key) && identical(.euc2_memo$key, key)) return(.euc2_memo$val)
  d <- .uc_design(meta, formula, tested_term)
  if (is.null(d)) { val <- .eu_fit_ad(counts, meta, formula, tested_term); val$fallback <- TRUE
  } else {
    depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts); g1 <- d$g1; Z <- d$Z
    rho <- round(.rho_design(depth, g1), 2)
    set.seed(7L); rows <- if (nrow(counts) > 300L) sort(sample.int(nrow(counts), 300L)) else seq_len(nrow(counts))
    ct <- counts[rows, , drop = FALSE]; ct <- ct[rowSums(ct > 0) >= 3, , drop = FALSE]
    med <- function(gm) { u <- .u_stat_cov(ct, depth, g1, Z, gm, "log", rho = rho, vcorr = TRUE); stats::median(u$U / sqrt(u$v), na.rm = TRUE) }
    flo <- med(-2); fhi <- med(2)
    gm <- if (is.finite(flo) && is.finite(fhi) && sign(flo) != sign(fhi)) stats::uniroot(med, c(-2, 2), f.lower = flo, f.upper = fhi, tol = 2e-3)$root else 0
    ud <- .u_stat_cov(counts, depth, g1, Z, gm, "det", rho = rho, vcorr = TRUE); ul <- .u_stat_cov(counts, depth, g1, Z, gm, "log", rho = rho, vcorr = TRUE)
    ok <- rowSums(counts > 0) >= 3 & ud$v > 0 & ul$v > 0; n1 <- sum(g1); n0 <- ncol(ud$b)
    cv <- (rowSums(ud$a * ul$a) / (n1 - 1)) / n1 + (rowSums(ud$b * ul$b) / (n0 - 1)) / n0
    tz <- function(t, df) stats::qnorm(stats::pt(t, df))
    z1 <- ifelse(ok, tz(ud$U / sqrt(ud$v), ud$df), NA); z2 <- ifelse(ok, tz(ul$U / sqrt(ul$v), ul$df), NA)
    rho_c <- ifelse(ok, pmin(pmax(cv / sqrt(ud$v * ul$v), -0.999), 0.999), NA)
    val <- list(feature = rownames(counts), p_det = 2 * stats::pnorm(-abs(z1)), p_log = 2 * stats::pnorm(-abs(z2)),
                p_max = .pmax2(pmax(abs(z1), abs(z2)), rho_c), est = ul$U / log(2), gam = gm, rho = rho, fallback = FALSE)
  }
  .euc2_memo$key <- key; .euc2_memo$val <- val; val
}
register_candidate("erdl_uc2", function(counts, meta, formula, tested_term) {
  v <- .euc_fit2(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_max, estimate = v$est)
}, notes = "it16: erdl_uc with the pair-residual double count removed from the U-statistic variance")
register_candidate("erd_uc2", function(counts, meta, formula, tested_term) {
  v <- .euc_fit2(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_det)
}, notes = "it16: detection-only erdl_uc2")

# --- it17: balanced designs -> raw values + log-depth covariate; confounded -> pairwise thinning ---
.RAW_D <- 0.25
.depth_d <- function(depth, g1) { l <- log(depth)
  s <- sqrt(((sum(g1) - 1) * stats::var(l[g1]) + (sum(!g1) - 1) * stats::var(l[!g1])) / (length(l) - 2))
  abs(mean(l[g1]) - mean(l[!g1])) / max(s, 1e-8) }
.euc3_memo <- new.env()
.euc_fit3 <- function(counts, meta, formula, tested_term) {
  key <- list(counts, meta, formula, tested_term)
  if (!is.null(.euc3_memo$key) && identical(.euc3_memo$key, key)) return(.euc3_memo$val)
  d <- .uc_design(meta, formula, tested_term)
  depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts)
  if (is.null(d) || .depth_d(depth, d$g1) > .RAW_D) { val <- .euc_fit(counts, meta, formula, tested_term); val$mode <- "thinned"
  } else {
    g1 <- d$g1; Z <- d$Z
    set.seed(7L); rows <- if (nrow(counts) > 300L) sort(sample.int(nrow(counts), 300L)) else seq_len(nrow(counts))
    ct <- counts[rows, , drop = FALSE]; ct <- ct[rowSums(ct > 0) >= 3, , drop = FALSE]
    med <- function(gm) { u <- .u_stat_cov(ct, depth, g1, Z, gm, "log", raw = TRUE); stats::median(u$U / sqrt(u$v), na.rm = TRUE) }
    flo <- med(-2); fhi <- med(2)
    gm <- if (is.finite(flo) && is.finite(fhi) && sign(flo) != sign(fhi)) stats::uniroot(med, c(-2, 2), f.lower = flo, f.upper = fhi, tol = 2e-3)$root else 0
    ud <- .u_stat_cov(counts, depth, g1, Z, gm, "det", raw = TRUE); ul <- .u_stat_cov(counts, depth, g1, Z, gm, "log", raw = TRUE)
    ok <- rowSums(counts > 0) >= 3 & ud$v > 0 & ul$v > 0; n1 <- sum(g1); n0 <- ncol(ud$b)
    cv <- (rowSums(ud$a * ul$a) / (n1 - 1)) / n1 + (rowSums(ud$b * ul$b) / (n0 - 1)) / n0
    tz <- function(t, df) stats::qnorm(stats::pt(t, df))
    z1 <- ifelse(ok, tz(ud$U / sqrt(ud$v), ud$df), NA); z2 <- ifelse(ok, tz(ul$U / sqrt(ul$v), ul$df), NA)
    rc <- ifelse(ok, pmin(pmax(cv / sqrt(ud$v * ul$v), -0.999), 0.999), NA)
    val <- list(feature = rownames(counts), p_det = 2 * stats::pnorm(-abs(z1)), p_log = 2 * stats::pnorm(-abs(z2)),
                p_max = .pmax2(pmax(abs(z1), abs(z2)), rc), est = ul$U / log(2), gam = gm, mode = "raw", fallback = FALSE)
  }
  .euc3_memo$key <- key; .euc3_memo$val <- val; val
}
register_candidate("erdl_ur", function(counts, meta, formula, tested_term) {
  v <- .euc_fit3(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_max, estimate = v$est)
}, notes = "it17: erdl_uc, but depth-balanced designs use raw values with a log-depth covariate (offset by gamma) instead of thinning")
register_candidate("erd_ur", function(counts, meta, formula, tested_term) {
  v <- .euc_fit3(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_det)
}, notes = "it17: detection-only erdl_ur")

# --- it18: repeated measures inside the pairwise test; trimmed compositional centre ---------------
# R23 (5 visits per subject, subject-level exposure) was erdlu_uc's biggest deficit (relative power
# 0.25 vs ADAPT 0.98): subjects forced the global-minimum LM fallback. Now the pairwise test runs with
# a subject-clustered variance. R11 (all DA in one direction): the compositional factor was the root
# of median z over ALL taxa, which one-directional DA pulls; a second pass re-solves it on the taxa
# with |z| < 2 at the first-pass root (a trimmed centre, like the empirical null's).
.uk_design <- function(meta, formula, tested_term) {
  tl <- attr(stats::terms(formula), "term.labels"); g <- meta[[tested_term]]
  if (length(unique(g)) != 2L) return(NULL)
  X <- stats::model.matrix(formula, meta); asg <- attr(X, "assign")
  zc <- which(asg != 0 & asg != which(tl == tested_term))
  list(g1 = g == sort(unique(g))[2], Z = if (length(zc)) X[, zc, drop = FALSE] else NULL, cl = .cluster_of(meta))
}
.euk_memo <- new.env()
.euk_fit <- function(counts, meta, formula, tested_term) {
  key <- list(counts, meta, formula, tested_term)
  if (!is.null(.euk_memo$key) && identical(.euk_memo$key, key)) return(.euk_memo$val)
  d <- .uk_design(meta, formula, tested_term)
  if (is.null(d)) { val <- .eu_fit_ad(counts, meta, formula, tested_term); val$fallback <- TRUE
  } else {
    depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts); g1 <- d$g1; Z <- d$Z; cl <- d$cl
    rho <- round(.rho_design(depth, g1), 2)
    set.seed(7L); rows <- if (nrow(counts) > 300L) sort(sample.int(nrow(counts), 300L)) else seq_len(nrow(counts))
    ct <- counts[rows, , drop = FALSE]; ct <- ct[rowSums(ct > 0) >= 3, , drop = FALSE]
    zf <- function(gm) { u <- .u_stat_cov(ct, depth, g1, Z, gm, "log", rho = rho, cl = cl); u$U / sqrt(u$v) }
    root <- function(sel) { med <- function(gm) stats::median(zf(gm)[sel], na.rm = TRUE); flo <- med(-2); fhi <- med(2)
      if (is.finite(flo) && is.finite(fhi) && sign(flo) != sign(fhi)) stats::uniroot(med, c(-2, 2), f.lower = flo, f.upper = fhi, tol = 2e-3)$root else 0 }
    g0 <- root(rep(TRUE, nrow(ct))); z0 <- zf(g0); sel <- is.finite(z0) & abs(z0) < 2
    gm <- if (sum(sel) >= 30) root(sel) else g0
    ud <- .u_stat_cov(counts, depth, g1, Z, gm, "det", rho = rho, cl = cl); ul <- .u_stat_cov(counts, depth, g1, Z, gm, "log", rho = rho, cl = cl)
    ok <- rowSums(counts > 0) >= 3 & ud$v > 0 & ul$v > 0; n1 <- sum(g1); n0 <- ncol(ud$b)
    if (is.null(cl)) { cv <- (rowSums(ud$a * ul$a) / (n1 - 1)) / n1 + (rowSums(ud$b * ul$b) / (n0 - 1)) / n0
    } else { Ma <- stats::model.matrix(~ ud$ca - 1); Mb <- stats::model.matrix(~ ud$cb - 1); Ga <- ncol(Ma); Gb <- ncol(Mb)
      cv <- rowSums((ud$a %*% Ma) * (ul$a %*% Ma)) / n1^2 * Ga / (Ga - 1) + rowSums((ud$b %*% Mb) * (ul$b %*% Mb)) / n0^2 * Gb / (Gb - 1) }
    tz <- function(t, df) stats::qnorm(stats::pt(t, df))
    z1 <- ifelse(ok, tz(ud$U / sqrt(ud$v), ud$df), NA); z2 <- ifelse(ok, tz(ul$U / sqrt(ul$v), ul$df), NA)
    rc <- ifelse(ok, pmin(pmax(cv / sqrt(ud$v * ul$v), -0.999), 0.999), NA)
    val <- list(feature = rownames(counts), p_det = 2 * stats::pnorm(-abs(z1)), p_log = 2 * stats::pnorm(-abs(z2)),
                p_max = .pmax2(pmax(abs(z1), abs(z2)), rc), est = ul$U / log(2), gam = gm, gam0 = g0, rho = rho, fallback = FALSE)
  }
  .euk_memo$key <- key; .euk_memo$val <- val; val
}
register_candidate("erdl_uk", function(counts, meta, formula, tested_term) {
  v <- .euk_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_max, estimate = v$est)
}, notes = "it18: erdl_uc + subject-clustered pairwise variance (repeated measures) + trimmed compositional centre")
register_candidate("erd_uk", function(counts, meta, formula, tested_term) {
  v <- .euk_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_det)
}, notes = "it18: detection-only erdl_uk")

# --- it19: the pairwise-rarefied probabilistic index (PI) kernel -----------------------------------
# h_ij = P(Y_i > Y_j) - P(Y_i < Y_j), Y_i and Y_j INDEPENDENT random rarefactions of the two samples to
# the pair's common depth (the exposed one at D/c). Under H0 both are draws from the same distribution
# given the composition, so E h_ij = 0 exactly -- a Mann-Whitney comparison done at a common depth.
# (Rank tests on ERD failed in it5/it6 because samples of different precision were ranked together;
# here each comparison is between two equally thinned counts.) With one sample zero it reduces to the
# detection kernel; for abundant taxa it is a rank comparison of abundance, robust to the skew that
# weakens a mean of log counts. Exact from the hypergeometric pmfs when mu_i + mu_j <= 40, normal
# approximation with continuity correction beyond; both-zero pairs contribute 0.
.pi_kernel <- function(yi, Ni, Di, Yj, Nj, Dj, kmax = 90L) {
  m <- length(yi); nk <- ncol(Yj); H <- matrix(0, m, nk)
  Yi <- matrix(yi, m, nk); NI <- Ni; DI <- matrix(Di, m, nk, byrow = TRUE); NJ <- matrix(Nj, m, nk, byrow = TRUE); DJ <- matrix(Dj, m, nk, byrow = TRUE)
  pz <- function(y, N, D) exp(lchoose(N - y, D) - lchoose(N, D))              # P(thinned count = 0)
  a0 <- Yi == 0; b0 <- Yj == 0
  one <- !a0 & b0; H[one] <- 1 - pz(Yi[one], NI, DI[one])
  two <- a0 & !b0; H[two] <- -(1 - pz(Yj[two], NJ[two], DJ[two]))
  both <- which(!a0 & !b0)
  if (length(both)) {
    yb <- Yi[both]; yjb <- Yj[both]; di <- DI[both]; dj <- DJ[both]; nj <- NJ[both]
    mi <- di * yb / NI; mj <- dj * yjb / nj
    vi <- di * (yb / NI) * (1 - yb / NI) * (NI - di) / max(NI - 1, 1); vj <- dj * (yjb / nj) * (1 - yjb / nj) * (nj - dj) / pmax(nj - 1, 1)
    big <- (mi + mj) > 40
    if (any(big)) { d <- mi[big] - mj[big]; s <- sqrt(pmax(vi[big] + vj[big], 1e-8))
      H[both[big]] <- (1 - stats::pnorm((0.5 - d) / s)) - stats::pnorm((-0.5 - d) / s) }
    sm <- which(!big)
    if (length(sm)) {
      K <- 0:kmax
      Pi <- outer(seq_along(sm), K, function(r, k) stats::dhyper(k, yb[sm][r], NI - yb[sm][r], di[sm][r]))
      Pj <- outer(seq_along(sm), K, function(r, k) stats::dhyper(k, yjb[sm][r], nj[sm][r] - yjb[sm][r], dj[sm][r]))
      Fi <- t(apply(Pi, 1, cumsum)); Fj <- t(apply(Pj, 1, cumsum))
      Fi1 <- cbind(0, Fi[, -ncol(Fi), drop = FALSE]); Fj1 <- cbind(0, Fj[, -ncol(Fj), drop = FALSE])   # F(k-1)
      H[both[sm]] <- rowSums(Pi * Fj1) - rowSums(Pj * Fi1)
    }
  }
  H
}

.epi_memo <- new.env()
.epi_fit <- function(counts, meta, formula, tested_term) {
  key <- list(counts, meta, formula, tested_term)
  if (!is.null(.epi_memo$key) && identical(.epi_memo$key, key)) return(.epi_memo$val)
  d <- .uk_design(meta, formula, tested_term)
  if (is.null(d)) { v <- .euk_fit(counts, meta, formula, tested_term); val <- list(feature = v$feature, p = v$p_max, fallback = TRUE)
  } else {
    depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts); g1 <- d$g1; Z <- d$Z; cl <- d$cl
    rho <- round(.rho_design(depth, g1), 2)
    gm <- .euk_fit(counts, meta, formula, tested_term)$gam                       # compositional factor from the log U-statistic
    u <- .u_stat_cov(counts, depth, g1, Z, gm, "pi", rho = rho, cl = cl)
    ok <- rowSums(counts > 0) >= 3 & u$v > 0
    val <- list(feature = rownames(counts), p = ifelse(ok, 2 * stats::pt(-abs(u$U / sqrt(u$v)), u$df), NA), est = u$U, gam = gm, fallback = FALSE)
  }
  .epi_memo$key <- key; .epi_memo$val <- val; val
}
register_candidate("pi_u", function(counts, meta, formula, tested_term) {
  v <- .epi_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p)
}, notes = "it19: pairwise-rarefied probabilistic index (Mann-Whitney at the pair's common depth), design rho, covariates, clusters")

# --- it20: permutation-calibrated tails in depth-balanced designs -------------------------------------
# Where the design is depth-balanced (rho = 1 by the it14 rule), has no covariates or subjects, and
# n <= 120, labels are exchangeable under H0, so the null distribution of the standardised pairwise
# statistic can be read off label permutations. All ordered-pair kernels K_ij = f_i(D_ij) - f_j(D_ij),
# D_ij = min(N_i, N_j), are computed once (up to 400 taxa); each permutation then only re-indexes K.
# The permuted z's, pooled over taxa, form an empirical null; each observed z (with its compositional
# centring, from erdl_uc) is mapped through it -- p = (1 + #{|z*| >= |z|}) / (1 + #z*), a Gaussian tail
# beyond the permuted range -- and back to z, then combined by the usual max test. This corrects the
# tail shape in either direction; the observed statistic itself is unchanged.
.pair_kernels <- function(counts, depth, h, rho = 1) {
  Y <- as.matrix(counts); storage.mode(Y) <- "double"; m <- nrow(Y); n <- ncol(Y)
  fcol <- function(Ysub, Nvec, Dvec) if (h == "det") .erd_f_var(Ysub, Nvec, Dvec) else .erl_f_var(Ysub, Nvec, Dvec)
  K <- array(0, c(m, n, n))
  for (i in seq_len(n - 1L)) { jj <- (i + 1L):n
    D <- pmax(floor(rho * pmin(depth[i], depth[jj])), 1)
    H <- fcol(Y[, rep(i, length(jj)), drop = FALSE], rep(depth[i], length(jj)), D) - fcol(Y[, jj, drop = FALSE], depth[jj], D)
    K[, i, jj] <- H; K[, jj, i] <- -H }
  K
}
.perm_z <- function(K, g1, B = 100L, seed = 11L) {
  m <- dim(K)[1]; n <- dim(K)[2]; n1 <- sum(g1); out <- numeric(0)
  set.seed(seed)
  for (b in seq_len(B)) { gb <- sample(g1); C <- which(gb); Tt <- which(!gb)
    A <- K[, C, Tt, drop = FALSE]                                            # m x n1 x n0
    a <- apply(A, c(1, 2), mean); bb <- apply(A, c(1, 3), mean); U <- rowMeans(a)
    va <- apply(a, 1, stats::var) / length(C); vb <- apply(bb, 1, stats::var) / length(Tt); v <- va + vb
    df <- v^2 / (va^2 / (length(C) - 1) + vb^2 / (length(Tt) - 1))
    z <- stats::qnorm(stats::pt(U / sqrt(v), df)); out <- c(out, z[is.finite(z)]) }
  out
}
.calib_z <- function(z, zstar) {
  az <- sort(abs(zstar)); N <- length(az); s_t <- stats::quantile(az, 0.99) / stats::qnorm(0.995)
  p <- vapply(abs(z), function(x) if (!is.finite(x)) NA_real_ else if (x > az[N]) 2 * stats::pnorm(-x / s_t) else
    (1 + N - findInterval(x, az, left.open = TRUE)) / (1 + N), numeric(1))
  sign(z) * stats::qnorm(1 - pmin(pmax(p, 1e-300), 1) / 2)
}
.eup_memo <- new.env()
.eup_fit <- function(counts, meta, formula, tested_term, B = 100L, sub = 400L) {
  key <- list(counts, meta, formula, tested_term)
  if (!is.null(.eup_memo$key) && identical(.eup_memo$key, key)) return(.eup_memo$val)
  v <- .euc_fit(counts, meta, formula, tested_term)
  d <- .uc_design(meta, formula, tested_term); depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts)
  use <- !is.null(d) && is.null(d$Z) && ncol(counts) <= 120L && isTRUE(v$rho == 1)
  if (!use) { v$perm <- FALSE; .eup_memo$key <- key; .eup_memo$val <- v; return(v) }
  set.seed(7L); rows <- which(rowSums(counts > 0) >= 3); if (length(rows) > sub) rows <- sort(sample(rows, sub))
  zs_d <- .perm_z(.pair_kernels(counts[rows, , drop = FALSE], depth, "det"), d$g1, B)
  zs_l <- .perm_z(.pair_kernels(counts[rows, , drop = FALSE], depth, "log"), d$g1, B)
  zd <- stats::qnorm(pmax(v$p_det, 1e-300) / 2, lower.tail = FALSE); zl <- stats::qnorm(pmax(v$p_log, 1e-300) / 2, lower.tail = FALSE)
  zd <- .calib_z(zd, zs_d); zl <- .calib_z(zl, zs_l)
  # correlation for the max test: recover it from the uncalibrated pair (same as erdl_uc used)
  rc <- attr(v, "rho_c"); if (is.null(rc)) rc <- rep(0.7, length(zd))
  pm <- .pmax2(pmax(abs(zd), abs(zl)), rc)
  val <- v; val$p_max <- pm; val$p_det <- 2 * stats::pnorm(-abs(zd)); val$p_log <- 2 * stats::pnorm(-abs(zl)); val$perm <- TRUE
  val$tail <- c(det = sd(zs_d), log = sd(zs_l), det99 = unname(stats::quantile(abs(zs_d), .999)), log99 = unname(stats::quantile(abs(zs_l), .999)))
  .eup_memo$key <- key; .eup_memo$val <- val; val
}
register_candidate("erdl_up", function(counts, meta, formula, tested_term) {
  v <- .eup_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_max, estimate = v$est)
}, notes = "it20: erdl_uc with permutation-calibrated component tails in depth-balanced designs")

# --- it21: depth-weighted pairs ---
.eucw_memo <- new.env()
.euc_fitw <- function(counts, meta, formula, tested_term, wk = 0.5) {
  key <- list(counts, meta, formula, tested_term, wk)
  if (!is.null(.eucw_memo$key) && identical(.eucw_memo$key, key)) return(.eucw_memo$val)
  d <- .uc_design(meta, formula, tested_term)
  if (is.null(d)) { val <- .euc_fit(counts, meta, formula, tested_term)
  } else {
    depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts); g1 <- d$g1; Z <- d$Z
    rho <- round(.rho_design(depth, g1), 2)
    set.seed(7L); rows <- if (nrow(counts) > 300L) sort(sample.int(nrow(counts), 300L)) else seq_len(nrow(counts))
    ct <- counts[rows, , drop = FALSE]; ct <- ct[rowSums(ct > 0) >= 3, , drop = FALSE]
    med <- function(gm) { u <- .u_stat_cov(ct, depth, g1, Z, gm, "log", rho = rho, wk = wk); stats::median(u$U / sqrt(u$v), na.rm = TRUE) }
    flo <- med(-2); fhi <- med(2)
    gm <- if (is.finite(flo) && is.finite(fhi) && sign(flo) != sign(fhi)) stats::uniroot(med, c(-2, 2), f.lower = flo, f.upper = fhi, tol = 2e-3)$root else 0
    ud <- .u_stat_cov(counts, depth, g1, Z, gm, "det", rho = rho, wk = wk); ul <- .u_stat_cov(counts, depth, g1, Z, gm, "log", rho = rho, wk = wk)
    ok <- rowSums(counts > 0) >= 3 & ud$v > 0 & ul$v > 0; n1 <- sum(g1); n0 <- ncol(ud$b)
    cv <- (rowSums(ud$a * ul$a) / (n1 - 1)) / n1 + (rowSums(ud$b * ul$b) / (n0 - 1)) / n0
    tz <- function(t, df) stats::qnorm(stats::pt(t, df))
    z1 <- ifelse(ok, tz(ud$U / sqrt(ud$v), ud$df), NA); z2 <- ifelse(ok, tz(ul$U / sqrt(ul$v), ul$df), NA)
    rho_c <- ifelse(ok, cv / sqrt(ud$v * ul$v), NA)
    val <- list(feature = rownames(counts), p_det = 2 * stats::pnorm(-abs(z1)), p_log = 2 * stats::pnorm(-abs(z2)),
                p_max = .pmax2(pmax(abs(z1), abs(z2)), rho_c), est = ul$U / log(2), gam = gm, rho = rho, fallback = FALSE)
    attr(val, "rho_c") <- rho_c
  }
  .eucw_memo$key <- key; .eucw_memo$val <- val; val
}
for (.wk in c(0.5, 1)) local({ wk <- .wk; tag <- if (wk == 0.5) "w05" else "w10"
  register_candidate(paste0("erdl_uc_", tag), function(counts, meta, formula, tested_term) {
    v <- .euc_fitw(counts, meta, formula, tested_term, wk = wk); data.frame(feature = v$feature, p = v$p_max, estimate = v$est)
  }, notes = sprintf("it21: erdl_uc with pair weights (common depth)^%.1f", wk))
})

# --- it22: censored (Tobit-score) kernel ---------------------------------------------------------------
# ADAPT's likely edge on unconfounded data: it treats a zero as "below detection" (left-censored log
# abundance) instead of as a value, so the detection and magnitude information enter ONE location
# parameter instead of two tests joined by a max. Kept inside the exact pairwise frame: at the pair's
# common depth D both thinned counts go through the same transform
#     g_D(y) = log y (y >= 1),   g_D(0) = E[l | l < t],  l ~ N(log D + mu, s^2),  t = log 0.5,
# i.e. a zero is imputed at the conditional mean of the latent log count below the detection limit --
# exactly the Tobit score for a location shift of l. (mu, s) per feature come from a label-free Tobit
# fit of log relative abundance over all samples at their own depths (zero <=> a < log(0.5 / N_i)),
# so g_D is a fixed function of the pooled data and h_ij = E g_D(Y_i) - E g_D(Y_j) keeps mean zero
# under H0 whatever the model's fit. Expanding, h = dL0 + w_D * dDet with w_D = -E[l | l < t] > 0:
# the log-count kernel (log y on positives) plus the detection kernel weighted by how far below
# detection the feature's zeros lie at that depth -- rare features lean on detection, common ones on
# magnitude, continuously, with no multiplicity cost.
.tobit_fit <- function(counts, depth, t = log(0.5), iters = 150L, smin = 0.3, tol = 1e-5) {
  Y <- as.matrix(counts); pos <- Y > 0; m <- nrow(Y)
  A <- log(Y) - rep(log(depth), each = m); A[!pos] <- 0
  cth <- matrix(t - log(depth), m, ncol(Y), byrow = TRUE)
  np <- pmax(rowSums(pos), 1); mu <- rowSums(A * pos) / np
  s <- pmax(sqrt(rowSums(((A - mu) * pos)^2) / pmax(np - 1, 1)), 1)
  for (it in seq_len(iters)) {
    al <- (cth - mu) / s; lam <- exp(stats::dnorm(al, log = TRUE) - stats::pnorm(al, log.p = TRUE))
    Ea <- ifelse(pos, A, mu - s * lam); Va <- ifelse(pos, 0, s^2 * pmax(1 - al * lam - lam^2, 0))
    mu2 <- rowMeans(Ea); s2 <- pmax(sqrt(rowMeans((Ea - mu2)^2 + Va)), smin)
    d <- max(abs(mu2 - mu), abs(s2 - s)); mu <- mu2; s <- s2; if (d < tol) break }
  list(mu = mu, s = s, t = t, iters = it)
}
.cen_kernel <- function(Yi, Ni, Di, Yj, Nj, Dj, tob) {
  Hd <- .erd_f_var(Yi, Ni, Di) - .erd_f_var(Yj, Nj, Dj)
  Hl <- .erl_f_var(Yi, Ni, Di, kmax = 130L, big_mu = 50, off = 0) - .erl_f_var(Yj, Nj, Dj, kmax = 130L, big_mu = 50, off = 0)
  mD <- outer(tob$mu, log(Dj), "+"); al <- (tob$t - mD) / tob$s
  lam <- exp(stats::dnorm(al, log = TRUE) - stats::pnorm(al, log.p = TRUE))
  Hl - (mD - tob$s * lam) * Hd                                                # dL0 + w_D dDet
}
.cov_uu <- function(u1, u2, n1) { n0 <- ncol(u1$b); (rowSums(u1$a * u2$a) / (n1 - 1)) / n1 + (rowSums(u1$b * u2$b) / (n0 - 1)) / n0 }
.ecen_memo <- new.env()
.ecen_fit <- function(counts, meta, formula, tested_term, t = log(0.5)) {
  key <- list(counts, meta, formula, tested_term, t)
  if (!is.null(.ecen_memo$key) && identical(.ecen_memo$key, key)) return(.ecen_memo$val)
  d <- .uc_design(meta, formula, tested_term)
  if (is.null(d)) { v <- .euc_fit(counts, meta, formula, tested_term); val <- v; val$p_cen <- val$p_dc <- val$p_lc <- v$p_max
  } else {
    depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts); g1 <- d$g1; Z <- d$Z
    v0 <- .euc_fit(counts, meta, formula, tested_term); gm <- v0$gam; rho <- v0$rho      # same centre and rho as the lead
    tob <- .tobit_fit(counts, depth, t)
    ud <- .u_stat_cov(counts, depth, g1, Z, gm, "det", rho = rho); ul <- .u_stat_cov(counts, depth, g1, Z, gm, "log", rho = rho)
    uc <- .u_stat_cov(counts, depth, g1, Z, gm, "cen", rho = rho, tob = tob)
    ok <- rowSums(counts > 0) >= 3 & ud$v > 0 & ul$v > 0 & uc$v > 0; n1 <- sum(g1)
    tz <- function(u) ifelse(ok, stats::qnorm(stats::pt(u$U / sqrt(u$v), u$df)), NA)
    zd <- tz(ud); zl <- tz(ul); zc <- tz(uc)
    r_dc <- ifelse(ok, .cov_uu(ud, uc, n1) / sqrt(ud$v * uc$v), NA); r_lc <- ifelse(ok, .cov_uu(ul, uc, n1) / sqrt(ul$v * uc$v), NA)
    val <- list(feature = rownames(counts), p_cen = 2 * stats::pnorm(-abs(zc)), p_max = v0$p_max,
                p_dc = .pmax2(pmax(abs(zd), abs(zc)), r_dc), p_lc = .pmax2(pmax(abs(zl), abs(zc)), r_lc),
                est = v0$est, z = cbind(det = zd, log = zl, cen = zc), gam = gm, rho = rho, tob = tob, fallback = FALSE)
  }
  .ecen_memo$key <- key; .ecen_memo$val <- val; val
}
register_candidate("cen_uc", function(counts, meta, formula, tested_term) {
  v <- .ecen_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_cen, estimate = v$est)
}, notes = "it22: pairwise U-statistic on the censored (Tobit-score) kernel: log y on positives, zeros at E[latent | below detection]")
register_candidate("cend_uc", function(counts, meta, formula, tested_term) {
  v <- .ecen_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_dc, estimate = v$est)
}, notes = "it22: max(detection, censored) pairwise test")

# --- it23: depth as a confounder of detection, not only of read sampling ---------------------------------
# mid R18/R19 (depth x4 / x9 with the exposure): the lead's FDR is 0.107 / 0.183 on MIDASim alone. The
# false positives are rare null taxa (control prevalence 0.02-0.18) present in 56-100% of the deep case
# samples, with 3-24x higher mean relative abundance there: MIDASim ties presence to library size through
# its copula thresholds, so a deep sample is NOT a thinned-up shallow one -- the premise of thinning.
# (Real pipelines do the same to a lesser degree: singleton removal, per-sample abundance thresholds.)
# ADAPT survives because its censoring at 1/N_i happens to match that generator. A remedy that assumes no
# model for how depth moves detection: treat depth as a confounder of the pair comparison and adjust for
# it by design -- (a) the pair's log-depth difference as a covariate of the pairwise regression (theta is
# the effect at equal depth), (b) overlap weights exp(-(dlogN)^2 / 2h^2) that concentrate on depth-matched
# pairs, or both. Both are functions of depths and labels only, so the thinning-exact null is untouched.
.eucd_memo <- new.env()
.eucd_fit <- function(counts, meta, formula, tested_term, mode = "cov", lo = 0.25, hmult = 1, rho_fix = NULL) {
  key <- list(counts, meta, formula, tested_term, mode, lo, hmult, rho_fix)
  if (!is.null(.eucd_memo$key) && identical(.eucd_memo$key, key)) return(.eucd_memo$val)
  d <- .uc_design(meta, formula, tested_term)
  depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts)
  dd <- if (is.null(d)) 0 else .depth_d(depth, d$g1)
  if (is.null(d) || dd <= lo) { val <- .euc_fit(counts, meta, formula, tested_term); val$dd <- dd
  } else {
    g1 <- d$g1; Z <- d$Z; ld <- log(depth)
    rho <- if (!is.null(rho_fix)) rho_fix else round(.rho_design(depth, g1), 2)
    if (mode %in% c("cov", "wcov")) Z <- cbind(Z, ldepth = ld)
    wfun <- NULL
    if (mode %in% c("w", "wcov")) {
      s <- sqrt(((sum(g1) - 1) * stats::var(ld[g1]) + (sum(!g1) - 1) * stats::var(ld[!g1])) / (length(ld) - 2)); h <- hmult * s
      wfun <- function(i, jj) exp(-(ld[i] - ld[jj])^2 / (2 * h^2)) }
    set.seed(7L); rows <- if (nrow(counts) > 300L) sort(sample.int(nrow(counts), 300L)) else seq_len(nrow(counts))
    ct <- counts[rows, , drop = FALSE]; ct <- ct[rowSums(ct > 0) >= 3, , drop = FALSE]
    med <- function(gm) { u <- .u_stat_cov(ct, depth, g1, Z, gm, "log", rho = rho, wfun = wfun); stats::median(u$U / sqrt(u$v), na.rm = TRUE) }
    flo <- med(-2); fhi <- med(2)
    gm <- if (is.finite(flo) && is.finite(fhi) && sign(flo) != sign(fhi)) stats::uniroot(med, c(-2, 2), f.lower = flo, f.upper = fhi, tol = 2e-3)$root else 0
    ud <- .u_stat_cov(counts, depth, g1, Z, gm, "det", rho = rho, wfun = wfun); ul <- .u_stat_cov(counts, depth, g1, Z, gm, "log", rho = rho, wfun = wfun)
    ok <- rowSums(counts > 0) >= 3 & ud$v > 0 & ul$v > 0; n1 <- sum(g1)
    tz <- function(u) ifelse(ok, stats::qnorm(stats::pt(u$U / sqrt(u$v), u$df)), NA)
    z1 <- tz(ud); z2 <- tz(ul); rho_c <- ifelse(ok, .cov_uu(ud, ul, n1) / sqrt(ud$v * ul$v), NA)
    val <- list(feature = rownames(counts), p_det = 2 * stats::pnorm(-abs(z1)), p_log = 2 * stats::pnorm(-abs(z2)),
                p_max = .pmax2(pmax(abs(z1), abs(z2)), rho_c), est = ul$U / log(2), gam = gm, rho = rho, dd = dd, fallback = FALSE)
  }
  .eucd_memo$key <- key; .eucd_memo$val <- val; val
}

# --- it24: variance moderation across taxa (empirical Bayes, limma-style) --------------------------------
# Where the lead is furthest behind with depth unconfounded, its per-taxon variance rests on few units:
# R01 (10 samples per group) and R23 (5 visits x 20 subjects, exposure between subjects, so the honest
# variance is cluster-robust with ~G-1 df; the lead even falls back to the global-minimum LM there, TPR
# 0.044 vs 0.17 for the methods that ignore subjects, whose null FPR is only 0.057-0.068). Borrow strength
# across taxa: treat each variance estimate as s0^2 * chisq_df / df around a prior fitted over all taxa
# (log-scale moments, limma's fitFDist), and use the posterior variance (d0 s0^2 + df v) / (d0 + df) with
# d0 + df degrees of freedom. Independent designs: prior = smooth trend in taxon prevalence. Clustered:
# the prior is the design effect -- v_cluster / v_independent -- pooled over taxa, so each taxon keeps its
# own precise independent-sample variance and only the clustering inflation is shared.
.trigamma_inv <- function(x) {
  y <- 0.5 + 1 / x
  for (i in 1:50) { tri <- trigamma(y); dif <- tri * (1 - tri / x) / psigamma(y, 2); y <- y + dif; if (max(-dif / y) < 1e-8) break }
  y }
.squeeze <- function(v, df, x = NULL, off = NULL, dfmax = Inf) {
  ok <- is.finite(v) & v > 0 & is.finite(df) & df > 0; if (is.null(off)) off <- rep(1, length(v))
  ok <- ok & is.finite(off) & off > 0
  e <- log(v / off) - digamma(df / 2) + log(df / 2); tr <- rep(NA_real_, length(v))
  if (!is.null(x) && sum(ok) >= 30) { sp <- stats::smooth.spline(x[ok], e[ok], df = 4); tr[ok] <- stats::predict(sp, x[ok])$y; np <- 4
  } else { tr[ok] <- mean(e[ok]); np <- 1 }
  n <- sum(ok); evar <- sum((e[ok] - tr[ok])^2) / max(n - np, 1) - mean(trigamma(df[ok] / 2))
  d0 <- if (evar > 0) 2 * .trigamma_inv(evar) else Inf
  s0 <- off * exp(tr + if (is.finite(d0)) digamma(d0 / 2) - log(d0 / 2) else 0)
  vp <- if (is.finite(d0)) (d0 * s0 + df * v) / (d0 + df) else s0
  dfp <- pmin(df + d0, dfmax)
  vp[!ok] <- v[!ok]; dfp[!ok] <- df[!ok]
  list(v = vp, df = dfp, d0 = d0)
}
.eum_memo <- new.env()
.eum_fit <- function(counts, meta, formula, tested_term, moderate = TRUE) {
  key <- list(counts, meta, formula, tested_term, moderate)
  if (!is.null(.eum_memo$key) && identical(.eum_memo$key, key)) return(.eum_memo$val)
  d <- .uk_design(meta, formula, tested_term)
  if (is.null(d)) { val <- .euc_fit(counts, meta, formula, tested_term); val$fallback <- TRUE
  } else {
    depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts); g1 <- d$g1; Z <- d$Z; cl <- d$cl
    rho <- round(.rho_design(depth, g1), 2)
    set.seed(7L); rows <- if (nrow(counts) > 300L) sort(sample.int(nrow(counts), 300L)) else seq_len(nrow(counts))
    ct <- counts[rows, , drop = FALSE]; ct <- ct[rowSums(ct > 0) >= 3, , drop = FALSE]
    med <- function(gm) { u <- .u_stat_cov(ct, depth, g1, Z, gm, "log", rho = rho, cl = cl); stats::median(u$U / sqrt(u$v), na.rm = TRUE) }
    flo <- med(-2); fhi <- med(2)
    gm <- if (is.finite(flo) && is.finite(fhi) && sign(flo) != sign(fhi)) stats::uniroot(med, c(-2, 2), f.lower = flo, f.upper = fhi, tol = 2e-3)$root else 0
    ud <- .u_stat_cov(counts, depth, g1, Z, gm, "det", rho = rho, cl = cl); ul <- .u_stat_cov(counts, depth, g1, Z, gm, "log", rho = rho, cl = cl)
    ok <- rowSums(counts > 0) >= 3 & ud$v > 0 & ul$v > 0; n1 <- sum(g1); n0 <- ncol(ud$b)
    if (is.null(cl)) { cv <- .cov_uu(ud, ul, n1)
    } else { Ma <- stats::model.matrix(~ ud$ca - 1); Mb <- stats::model.matrix(~ ud$cb - 1); Ga <- ncol(Ma); Gb <- ncol(Mb)
      cv <- rowSums((ud$a %*% Ma) * (ul$a %*% Ma)) / n1^2 * Ga / (Ga - 1) + rowSums((ud$b %*% Mb) * (ul$b %*% Mb)) / n0^2 * Gb / (Gb - 1) }
    rc <- ifelse(ok, pmin(pmax(cv / sqrt(ud$v * ul$v), -0.999), 0.999), NA)
    mod <- function(u) {
      if (!moderate) return(list(v = u$v, df = u$df, d0 = NA))
      if (is.null(cl)) .squeeze(ifelse(ok, u$v, NA), u$df, x = stats::qlogis(pmin(pmax(rowMeans(counts > 0), 0.01), 0.99)))
      else .squeeze(ifelse(ok, u$v, NA), u$df, off = u$vind, dfmax = u$dfind) }
    md <- mod(ud); ml <- mod(ul)
    tz <- function(u, m) ifelse(ok, stats::qnorm(stats::pt(u$U / sqrt(m$v), m$df)), NA)
    z1 <- tz(ud, md); z2 <- tz(ul, ml)
    val <- list(feature = rownames(counts), p_det = 2 * stats::pnorm(-abs(z1)), p_log = 2 * stats::pnorm(-abs(z2)),
                p_max = .pmax2(pmax(abs(z1), abs(z2)), rc), est = ul$U / log(2), gam = gm, rho = rho, d0 = c(det = md$d0, log = ml$d0),
                df_med = c(raw = stats::median(ul$df[ok]), mod = stats::median(ml$df[ok])), fallback = FALSE)
  }
  .eum_memo$key <- key; .eum_memo$val <- val; val
}
register_candidate("erdl_um", function(counts, meta, formula, tested_term) {
  v <- .eum_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_max, estimate = v$est)
}, notes = "it24: pairwise max test with empirical-Bayes variance moderation across taxa (trend in prevalence; clustered designs: pooled design effect)")
register_candidate("erdl_uk0", function(counts, meta, formula, tested_term) {
  v <- .eum_fit(counts, meta, formula, tested_term, moderate = FALSE); data.frame(feature = v$feature, p = v$p_max, estimate = v$est)
}, notes = "it24 control: same pairwise path (clusters handled pairwise), no moderation")

# --- it25: depth slope pooled across taxa -------------------------------------------------------------
# it23's per-taxon log-depth covariate fixes mid R18/R19 (null FPR 0.10 -> 0.05) but spends the whole
# depth-exposure overlap on every taxon: TP 4 -> 0 on house R19 and mid R18 twinsuk. The depth effect
# that breaks thinning on mid is systematic -- it depends on how rare a taxon is, not on which taxon --
# so estimate it once across taxa: per-taxon slopes b_j of the pair kernel on the pair's log-depth
# difference (from the covariate regression), smoothed against taxon prevalence (robust loess), and
# subtract m(x_j) * (log N_i - log N_j) from every pair. On a sampling-model simulator m ~ 0 and the test
# is the lead; on MIDASim m carries the rare-taxon inflation. The slope curve uses depths and pooled
# counts, not labels beyond the pairing, and costs no per-taxon degrees of freedom. Only when depth is
# confounded with the exposure (the rho trigger, d > 0.25).
.dslope_curve <- function(b, x) {
  ok <- is.finite(b) & is.finite(x); if (sum(ok) < 30) return(rep(0, length(b)))
  lo <- stats::loess(b ~ x, data = data.frame(b = b[ok], x = x[ok]), span = 0.6, degree = 1, family = "symmetric")
  m <- rep(0, length(b)); m[is.finite(x)] <- stats::predict(lo, data.frame(x = x[is.finite(x)])); m[!is.finite(m)] <- 0; m }
.eucs_memo <- new.env()
.eucs_fit <- function(counts, meta, formula, tested_term, rho_fix = NULL, lo = 0.25) {
  key <- list(counts, meta, formula, tested_term, rho_fix, lo)
  if (!is.null(.eucs_memo$key) && identical(.eucs_memo$key, key)) return(.eucs_memo$val)
  d <- .uc_design(meta, formula, tested_term)
  depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts)
  dd <- if (is.null(d)) 0 else .depth_d(depth, d$g1)
  if (is.null(d) || dd <= lo) { val <- .euc_fit(counts, meta, formula, tested_term); val$dd <- dd
  } else {
    g1 <- d$g1; Z <- d$Z; ld <- log(depth); v0 <- .euc_fit(counts, meta, formula, tested_term)
    rho <- if (!is.null(rho_fix)) rho_fix else v0$rho
    x <- stats::qlogis(pmin(pmax(rowMeans(counts > 0), 0.01), 0.99)); Zd <- cbind(Z, ldepth = ld)
    sl <- function(h) { u <- .u_stat_cov(counts, depth, g1, Zd, v0$gam, h, rho = rho); .dslope_curve(u$B[, ncol(u$B)], x) }
    md <- sl("det"); ml <- sl("log")
    set.seed(7L); rows <- if (nrow(counts) > 300L) sort(sample.int(nrow(counts), 300L)) else seq_len(nrow(counts))
    rows <- rows[rowSums(counts[rows, , drop = FALSE] > 0) >= 3]; ct <- counts[rows, , drop = FALSE]
    med <- function(gm) { u <- .u_stat_cov(ct, depth, g1, Z, gm, "log", rho = rho, dslope = ml[rows]); stats::median(u$U / sqrt(u$v), na.rm = TRUE) }
    flo <- med(-2); fhi <- med(2)
    gm <- if (is.finite(flo) && is.finite(fhi) && sign(flo) != sign(fhi)) stats::uniroot(med, c(-2, 2), f.lower = flo, f.upper = fhi, tol = 2e-3)$root else 0
    ud <- .u_stat_cov(counts, depth, g1, Z, gm, "det", rho = rho, dslope = md); ul <- .u_stat_cov(counts, depth, g1, Z, gm, "log", rho = rho, dslope = ml)
    ok <- rowSums(counts > 0) >= 3 & ud$v > 0 & ul$v > 0; n1 <- sum(g1)
    tz <- function(u) ifelse(ok, stats::qnorm(stats::pt(u$U / sqrt(u$v), u$df)), NA)
    z1 <- tz(ud); z2 <- tz(ul); rc <- ifelse(ok, .cov_uu(ud, ul, n1) / sqrt(ud$v * ul$v), NA)
    val <- list(feature = rownames(counts), p_det = 2 * stats::pnorm(-abs(z1)), p_log = 2 * stats::pnorm(-abs(z2)),
                p_max = .pmax2(pmax(abs(z1), abs(z2)), rc), est = ul$U / log(2), gam = gm, rho = rho, dd = dd,
                slope = cbind(det = md, log = ml), fallback = FALSE)
  }
  .eucs_memo$key <- key; .eucs_memo$val <- val; val
}
register_candidate("erdl_us", function(counts, meta, formula, tested_term) {
  v <- .eucs_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_max, estimate = v$est)
}, notes = "it25: erdl_uc with a depth slope pooled across taxa (smooth in prevalence) removed from every pair when depth is confounded")
register_candidate("erdl_us1", function(counts, meta, formula, tested_term) {
  v <- .eucs_fit(counts, meta, formula, tested_term, rho_fix = 1); data.frame(feature = v$feature, p = v$p_max, estimate = v$est)
}, notes = "it25: erdl_us without thinning below the pair's common depth (rho = 1)")

# --- it26: model-based transform, design-based test -----------------------------------------------------
# ADAPT's edge on unconfounded data is efficiency: every read of every sample enters a latent-abundance
# model, where the pairwise thinning keeps only the shallower member's depth. When depth is balanced
# across groups (the rho = 1 branch, d <= 0.25), ANY label-free per-sample transform gives a valid
# two-sample test by exchangeability, so the model can be used for the transform without being trusted
# for inference. Transform: the posterior mode of the latent log relative abundance a_ij under a
# Poisson-lognormal fitted to each taxon without labels, y ~ Pois(N_i e^a), a ~ N(mu_j, s_j^2)
# (Laplace EM), which shrinks low counts toward the taxon mean and places each zero according to how much
# a zero at THAT depth says -- a deep zero far below, a shallow zero near the mean. Compositional factor:
# exposed samples enter with N_i e^gamma (a -> a + gamma exactly on the latent scale), gamma chosen so the
# median z is 0 as before. Tested by the pairwise U-statistic machinery (covariates, clusters, moderation)
# on s_i - s_j, combined with the thinned detection test by the max. Depth-confounded designs keep the lead.
.pln_mode <- function(Y, Nm, mu, s2, iters = 25L) {
  a <- log((Y + 0.5) / Nm); a <- pmin(a, 0)
  for (it in seq_len(iters)) { ea <- Nm * exp(a); g <- Y - ea - (a - mu) / s2; h <- ea + 1 / s2
    st <- g / h; st <- pmax(pmin(st, 3), -3); a <- a + st; if (max(abs(st)) < 1e-7) break }
  list(a = a, v = 1 / (Nm * exp(a) + 1 / s2))
}
.pln_fit <- function(counts, depth, iters = 60L, smin = 0.2) {
  Y <- as.matrix(counts); m <- nrow(Y); Nm <- matrix(depth, m, ncol(Y), byrow = TRUE)
  pos <- Y > 0; mu <- log(pmax(rowSums(Y), 0.5) / sum(depth)); s2 <- rep(1, m)
  for (it in seq_len(iters)) { md <- .pln_mode(Y, Nm, mu, s2); mu2 <- rowMeans(md$a); s22 <- pmax(rowMeans((md$a - mu2)^2 + md$v), smin^2)
    dd <- max(abs(mu2 - mu), abs(sqrt(s22) - sqrt(s2))); mu <- mu2; s2 <- s22; if (dd < 1e-4) break }
  list(mu = mu, s2 = s2)
}
.pln_scores <- function(counts, depth, g1, gam, fit) {
  Y <- as.matrix(counts); Neff <- depth * exp(gam * g1); Nm <- matrix(Neff, nrow(Y), ncol(Y), byrow = TRUE)
  .pln_mode(Y, Nm, fit$mu, fit$s2)$a
}
.epln_memo <- new.env()
.epln_fit <- function(counts, meta, formula, tested_term, lo = 0.25, moderate = FALSE) {
  key <- list(counts, meta, formula, tested_term, lo, moderate)
  if (!is.null(.epln_memo$key) && identical(.epln_memo$key, key)) return(.epln_memo$val)
  d <- .uk_design(meta, formula, tested_term)
  depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts)
  dd <- if (is.null(d)) Inf else .depth_d(depth, d$g1)
  if (is.null(d) || dd > lo) { val <- if (moderate) .eum_fit(counts, meta, formula, tested_term) else .euc_fit(counts, meta, formula, tested_term)
    val$p_pln <- val$p_max; val$p_plnd <- val$p_max; val$branch <- "lead"
  } else {
    g1 <- d$g1; Z <- d$Z; cl <- d$cl; fit <- .pln_fit(counts, depth)
    set.seed(7L); rows <- if (nrow(counts) > 300L) sort(sample.int(nrow(counts), 300L)) else seq_len(nrow(counts))
    rows <- rows[rowSums(counts[rows, , drop = FALSE] > 0) >= 3]; fr <- list(mu = fit$mu[rows], s2 = fit$s2[rows])
    med <- function(gm) { S <- .pln_scores(counts[rows, , drop = FALSE], depth, g1, gm, fr)
      u <- .u_stat_cov(counts[rows, , drop = FALSE], depth, g1, Z, 0, "pre", S = S, cl = cl); stats::median(u$U / sqrt(u$v), na.rm = TRUE) }
    flo <- med(-2); fhi <- med(2)
    gm <- if (is.finite(flo) && is.finite(fhi) && sign(flo) != sign(fhi)) stats::uniroot(med, c(-2, 2), f.lower = flo, f.upper = fhi, tol = 2e-3)$root else 0
    S <- .pln_scores(counts, depth, g1, gm, fit)
    up <- .u_stat_cov(counts, depth, g1, Z, 0, "pre", S = S, cl = cl)
    ud <- .u_stat_cov(counts, depth, g1, Z, gm, "det", rho = 1, cl = cl)
    ok <- rowSums(counts > 0) >= 3 & ud$v > 0 & up$v > 0; n1 <- sum(g1)
    mod <- function(u) { if (!moderate) return(list(v = u$v, df = u$df))
      if (is.null(cl)) .squeeze(ifelse(ok, u$v, NA), u$df, x = stats::qlogis(pmin(pmax(rowMeans(counts > 0), 0.01), 0.99)))
      else .squeeze(ifelse(ok, u$v, NA), u$df, off = u$vind, dfmax = u$dfind) }
    if (!is.null(cl) && !moderate) { }  # clustered variance already in u$v
    mp <- mod(up); md <- mod(ud)
    tz <- function(u, m) ifelse(ok, stats::qnorm(stats::pt(u$U / sqrt(m$v), m$df)), NA)
    zp <- tz(up, mp); zd <- tz(ud, md)
    cv <- if (is.null(cl)) .cov_uu(ud, up, n1) else {
      Ma <- stats::model.matrix(~ ud$ca - 1); Mb <- stats::model.matrix(~ ud$cb - 1); Ga <- ncol(Ma); Gb <- ncol(Mb); n0 <- ncol(ud$b)
      rowSums((ud$a %*% Ma) * (up$a %*% Ma)) / n1^2 * Ga / (Ga - 1) + rowSums((ud$b %*% Mb) * (up$b %*% Mb)) / n0^2 * Gb / (Gb - 1) }
    rc <- ifelse(ok, pmin(pmax(cv / sqrt(ud$v * up$v), -0.999), 0.999), NA)
    val <- list(feature = rownames(counts), p_pln = 2 * stats::pnorm(-abs(zp)), p_det = 2 * stats::pnorm(-abs(zd)),
                p_plnd = .pmax2(pmax(abs(zp), abs(zd)), rc), est = up$U / log(2), gam = gm, dd = dd, branch = "pln", fallback = FALSE)
    val$p_max <- val$p_plnd
  }
  .epln_memo$key <- key; .epln_memo$val <- val; val
}
register_candidate("pln_u", function(counts, meta, formula, tested_term) {
  v <- .epln_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_pln, estimate = v$est)
}, notes = "it26: balanced depth -> pairwise test on Poisson-lognormal posterior-mode latent abundance (label-free fit); confounded -> lead")
register_candidate("plnd_u", function(counts, meta, formula, tested_term) {
  v <- .epln_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_plnd, estimate = v$est)
}, notes = "it26: max(thinned detection, PLN latent) when depth is balanced; lead otherwise")
register_candidate("plnd_um", function(counts, meta, formula, tested_term) {
  v <- .epln_fit(counts, meta, formula, tested_term, moderate = TRUE); data.frame(feature = v$feature, p = v$p_plnd, estimate = v$est)
}, notes = "it26 + it24: plnd_u with variance moderation (and pairwise clustered handling)")

# --- it27: empirical-Bayes depth adjustment ------------------------------------------------------------
# it23 (log-depth covariate per taxon) fixes mid R18/R19 but costs the whole depth-exposure overlap for
# every taxon; it25 (one slope curve pooled over taxa) is cheap but does not fix mid -- the depth effect
# there is taxon-specific. In between: the per-taxon depth slope b_j of the pair kernel (identified from
# the depth spread WITHIN the exposure groups) shrunk toward 0 with a prior variance tau^2 estimated
# across taxa of similar prevalence (method of moments, five bins): B_j = tau^2 / (tau^2 + se_j^2).
# Estimate theta_EB = theta_cov + (1 - B_j) b_j dbar -- the lead when B = 0 (OLS: mean h = theta +
# b dbar), the covariate-adjusted effect when B = 1. B depends on se_j, not on b_j, so theta_EB is
# linear in the pair kernels and its variance comes from the same projections. Where thinning holds
# (house, msq) tau^2 ~ 0 and this is the lead; on MIDASim it adjusts the taxa whose detection follows depth.
.eucb_memo <- new.env()
.eucb_fit <- function(counts, meta, formula, tested_term, lo = 0.25, nbin = 5L) {
  key <- list(counts, meta, formula, tested_term, lo, nbin)
  if (!is.null(.eucb_memo$key) && identical(.eucb_memo$key, key)) return(.eucb_memo$val)
  d <- .uc_design(meta, formula, tested_term)
  depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts)
  dd <- if (is.null(d)) 0 else .depth_d(depth, d$g1)
  v0 <- .euc_fit(counts, meta, formula, tested_term)
  if (is.null(d) || dd <= lo) { val <- v0; val$dd <- dd
  } else {
    g1 <- d$g1; Zd <- cbind(d$Z, ldepth = log(depth)); pb <- ncol(Zd) + 1L; rho <- v0$rho; gm <- v0$gam; n1 <- sum(g1)
    x <- stats::qlogis(pmin(pmax(rowMeans(counts > 0), 0.01), 0.99)); ok0 <- rowSums(counts > 0) >= 3
    eb <- function(h) {
      u1 <- .u_stat_cov(counts, depth, g1, Zd, gm, h, rho = rho, coef = 1L); u2 <- .u_stat_cov(counts, depth, g1, Zd, gm, h, rho = rho, coef = pb)
      b <- u1$B[, pb]; se2 <- u2$v; dbar <- u1$M[1, pb] / u1$M[1, 1]
      bins <- cut(x, unique(stats::quantile(x[ok0], seq(0, 1, length.out = nbin + 1))), include.lowest = TRUE)
      tau2 <- ave(ifelse(ok0 & is.finite(b), b^2 - se2, NA), bins, FUN = function(z) max(mean(z, na.rm = TRUE), 0))
      B <- ifelse(is.finite(tau2) & tau2 > 0, tau2 / (tau2 + se2), 0); w <- (1 - B) * dbar
      a <- u1$a + w * u2$a; bb <- u1$b + w * u2$b; n0 <- ncol(bb)
      va <- rowSums(a^2) / (n1 - 1) / n1; vb <- rowSums(bb^2) / (n0 - 1) / n0; v <- va + vb
      list(U = u1$U + w * b, v = v, df = v^2 / (va^2 / (n1 - 1) + vb^2 / (n0 - 1)), a = a, b = bb, B = B, tau2 = tapply(tau2, bins, `[`, 1), dbar = dbar) }
    ud <- eb("det"); ul <- eb("log")
    ok <- ok0 & ud$v > 0 & ul$v > 0
    tz <- function(u) ifelse(ok, stats::qnorm(stats::pt(u$U / sqrt(u$v), u$df)), NA)
    z1 <- tz(ud); z2 <- tz(ul); rc <- ifelse(ok, pmin(pmax(.cov_uu(ud, ul, n1) / sqrt(ud$v * ul$v), -0.999), 0.999), NA)
    val <- list(feature = rownames(counts), p_det = 2 * stats::pnorm(-abs(z1)), p_log = 2 * stats::pnorm(-abs(z2)),
                p_max = .pmax2(pmax(abs(z1), abs(z2)), rc), est = ul$U / log(2), gam = gm, rho = rho, dd = dd,
                Bmed = c(det = stats::median(ud$B[ok]), log = stats::median(ul$B[ok])), tau2 = rbind(det = ud$tau2, log = ul$tau2), fallback = FALSE)
  }
  .eucb_memo$key <- key; .eucb_memo$val <- val; val
}
register_candidate("erdl_ub", function(counts, meta, formula, tested_term) {
  v <- .eucb_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_max, estimate = v$est)
}, notes = "it27: erdl_uc with an empirical-Bayes per-taxon depth-slope adjustment when depth is confounded with the exposure")

# --- it28: the combination -- moderation (it24, all designs) + EB depth adjustment (it27, confounded) ---
.eb_depth <- function(counts, depth, g1, Z, gm, h, rho, nbin = 5L) {
  Zd <- cbind(Z, ldepth = log(depth)); pb <- ncol(Zd) + 1L; n1 <- sum(g1)
  x <- stats::qlogis(pmin(pmax(rowMeans(counts > 0), 0.01), 0.99)); ok0 <- rowSums(counts > 0) >= 3
  u1 <- .u_stat_cov(counts, depth, g1, Zd, gm, h, rho = rho, coef = 1L); u2 <- .u_stat_cov(counts, depth, g1, Zd, gm, h, rho = rho, coef = pb)
  b <- u1$B[, pb]; se2 <- u2$v; dbar <- u1$M[1, pb] / u1$M[1, 1]
  bins <- cut(x, unique(stats::quantile(x[ok0], seq(0, 1, length.out = nbin + 1))), include.lowest = TRUE)
  tau2 <- ave(ifelse(ok0 & is.finite(b), b^2 - se2, NA), bins, FUN = function(z) max(mean(z, na.rm = TRUE), 0))
  B <- ifelse(is.finite(tau2) & tau2 > 0, tau2 / (tau2 + se2), 0); w <- (1 - B) * dbar
  a <- u1$a + w * u2$a; bb <- u1$b + w * u2$b; n0 <- ncol(bb)
  va <- rowSums(a^2) / (n1 - 1) / n1; vb <- rowSums(bb^2) / (n0 - 1) / n0; v <- va + vb
  list(U = u1$U + w * b, v = v, df = v^2 / (va^2 / (n1 - 1) + vb^2 / (n0 - 1)), a = a, b = bb, B = B)
}
.eumb_memo <- new.env()
.eumb_fit <- function(counts, meta, formula, tested_term, lo = 0.25) {
  key <- list(counts, meta, formula, tested_term, lo)
  if (!is.null(.eumb_memo$key) && identical(.eumb_memo$key, key)) return(.eumb_memo$val)
  d <- .uk_design(meta, formula, tested_term)
  depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts)
  dd <- if (is.null(d)) 0 else .depth_d(depth, d$g1)
  if (is.null(d) || !is.null(d$cl) || dd <= lo) { val <- .eum_fit(counts, meta, formula, tested_term); val$dd <- dd; val$eb <- FALSE
  } else {
    v0 <- .euc_fit(counts, meta, formula, tested_term); g1 <- d$g1; n1 <- sum(g1)
    ud <- .eb_depth(counts, depth, g1, d$Z, v0$gam, "det", v0$rho); ul <- .eb_depth(counts, depth, g1, d$Z, v0$gam, "log", v0$rho)
    ok <- rowSums(counts > 0) >= 3 & ud$v > 0 & ul$v > 0
    xx <- stats::qlogis(pmin(pmax(rowMeans(counts > 0), 0.01), 0.99))
    md <- .squeeze(ifelse(ok, ud$v, NA), ud$df, x = xx); ml <- .squeeze(ifelse(ok, ul$v, NA), ul$df, x = xx)
    tz <- function(u, m) ifelse(ok, stats::qnorm(stats::pt(u$U / sqrt(m$v), m$df)), NA)
    z1 <- tz(ud, md); z2 <- tz(ul, ml); rc <- ifelse(ok, pmin(pmax(.cov_uu(ud, ul, n1) / sqrt(ud$v * ul$v), -0.999), 0.999), NA)
    val <- list(feature = rownames(counts), p_det = 2 * stats::pnorm(-abs(z1)), p_log = 2 * stats::pnorm(-abs(z2)),
                p_max = .pmax2(pmax(abs(z1), abs(z2)), rc), est = ul$U / log(2), gam = v0$gam, rho = v0$rho, dd = dd, eb = TRUE, fallback = FALSE)
  }
  .eumb_memo$key <- key; .eumb_memo$val <- val; val
}
register_candidate("erdl_umb", function(counts, meta, formula, tested_term) {
  v <- .eumb_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_max, estimate = v$est)
}, notes = "it28: erdl_um (moderation, pairwise clusters) + it27 EB depth adjustment when depth is confounded")

# --- it29: expected rarefied square-root count (ERS) -------------------------------------------------------
# srv9 with comparators on the dev suite: the lead's largest power deficits are on msq (GUniFrac
# SimulateMSeq: LDM 51/5 + 70/5 TP vs 40/1 + 44/1 at R00; R17 56/1 + 64/2 vs 32/2 + 43/1) and mid R00
# (LDM 20/3 + 35/6 vs 12/1 + 20/2). SimulateMSeq multiplies a real template sample's counts by the effect
# and re-draws them multinomially at a new depth: zeros are inherited (no prevalence signal), the signal is
# purely in the magnitude of low-to-moderate counts, where sampling noise matters. The log scale suits
# lognormal biological spread at high counts; for Poisson-like counts the variance-stabilising scale is the
# square root (LDM's arcsin-root). E[sqrt(Y_D)]: exact hypergeometric sum to D p = 50, fourth-order
# expansion above.
.ers_f_var <- function(counts, depth, Dvec, kmax = 70L, big_mu = 20) {
  Y <- as.matrix(counts); N <- matrix(depth, nrow(Y), ncol(Y), byrow = TRUE); Dm <- matrix(Dvec, nrow(Y), ncol(Y), byrow = TRUE)
  mu <- Dm * Y / N; out <- matrix(0, nrow(Y), ncol(Y))
  big <- mu > big_mu
  if (any(big)) { pp <- (Y / N)[big]; Nb <- N[big]; Db <- Dm[big]; fpc <- (Nb - Db) / pmax(Nb - 1, 1)
    v <- Db * pp * (1 - pp) * fpc; k3 <- v * (1 - 2 * pp) * (Nb - 2 * Db) / pmax(Nb - 2, 1); m1 <- mu[big]
    out[big] <- sqrt(m1) - v / (8 * m1^1.5) + k3 / (16 * m1^2.5) - 15 * (3 * v^2 + v * (1 - 6 * pp * (1 - pp))) / (384 * m1^3.5) }
  sm <- which(!big & Y > 0)
  if (length(sm)) { y <- Y[sm]; nn <- N[sm] - Y[sm]; dd <- Dm[sm]; acc <- numeric(length(sm))
    for (k in 1:kmax) { live <- k <= y & k <= dd; if (!any(live)) break
      acc[live] <- acc[live] + stats::dhyper(k, y[live], nn[live], dd[live]) * sqrt(k) }
    out[sm] <- acc }
  out
}

.pmax3 <- function(m, r12, r13, r23, n = 48L) {                   # P(max |Z_k| >= m), trivariate normal, vectorised over taxa
  gl <- .gauss_legendre(n); out <- rep(NA_real_, length(m))
  for (i in seq_along(m)) {
    if (!all(is.finite(c(m[i], r12[i], r13[i], r23[i])))) next
    a <- max(min(r12[i], 0.995), -0.995); b <- max(min(r13[i], 0.995), -0.995); c <- max(min(r23[i], 0.995), -0.995)
    R <- matrix(c(1, a, b, a, 1, c, b, c, 1), 3); ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) < 1e-4) { e2 <- eigen(R, symmetric = TRUE); R <- e2$vectors %*% diag(pmax(e2$values, 1e-4)) %*% t(e2$vectors); R <- stats::cov2cor(R); a <- R[1, 2]; b <- R[1, 3]; c <- R[2, 3] }
    s2 <- sqrt(1 - a^2); be1 <- (b - c * a) / (1 - a^2); be2 <- (c - b * a) / (1 - a^2); s3 <- sqrt(max(1 - (be1 * b + be2 * c), 1e-8))
    z1 <- m[i] * gl$nodes; z2 <- m[i] * gl$nodes; w <- m[i] * gl$weights
    Z1 <- matrix(z1, n, n); Z2 <- matrix(z2, n, n, byrow = TRUE)
    f2 <- stats::dnorm((Z2 - a * Z1) / s2) / s2; mu3 <- be1 * Z1 + be2 * Z2
    p3 <- stats::pnorm((m[i] - mu3) / s3) - stats::pnorm((-m[i] - mu3) / s3)
    inner <- (f2 * p3) %*% w
    out[i] <- 1 - sum(w * stats::dnorm(z1) * inner) }
  pmin(pmax(out, 0), 1)
}
# the lead's statistic with a third kernel: max(|z_det|, |z_sqrt|, |z_log|) against the trivariate normal
.eu3_memo <- new.env()
.eu3_fit <- function(counts, meta, formula, tested_term) {
  key <- list(counts, meta, formula, tested_term)
  if (!is.null(.eu3_memo$key) && identical(.eu3_memo$key, key)) return(.eu3_memo$val)
  d <- .uc_design(meta, formula, tested_term); v0 <- .euc_fit(counts, meta, formula, tested_term)
  if (is.null(d)) { val <- v0; val$p_sqrt <- v0$p_max; val$p3 <- v0$p_max; val$p_ls <- v0$p_max; val$p_ds <- v0$p_max
  } else {
    depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts); g1 <- d$g1; Z <- d$Z; n1 <- sum(g1)
    ud <- .u_stat_cov(counts, depth, g1, Z, v0$gam, "det", rho = v0$rho); ul <- .u_stat_cov(counts, depth, g1, Z, v0$gam, "log", rho = v0$rho)
    us <- .u_stat_cov(counts, depth, g1, Z, v0$gam, "sqrt", rho = v0$rho)
    ok <- rowSums(counts > 0) >= 3 & ud$v > 0 & ul$v > 0 & us$v > 0
    tz <- function(u) ifelse(ok, stats::qnorm(stats::pt(u$U / sqrt(u$v), u$df)), NA)
    zd <- tz(ud); zl <- tz(ul); zs <- tz(us)
    cr <- function(u1, u2) ifelse(ok, .cov_uu(u1, u2, n1) / sqrt(u1$v * u2$v), NA)
    rdl <- cr(ud, ul); rds <- cr(ud, us); rls <- cr(ul, us)
    val <- list(feature = rownames(counts), p_det = 2 * stats::pnorm(-abs(zd)), p_log = 2 * stats::pnorm(-abs(zl)), p_sqrt = 2 * stats::pnorm(-abs(zs)),
                p_max = v0$p_max, p_ls = .pmax2(pmax(abs(zl), abs(zs)), rls), p_ds = .pmax2(pmax(abs(zd), abs(zs)), rds),
                p3 = .pmax3(pmax(abs(zd), abs(zl), abs(zs)), rdl, rds, rls), est = v0$est, gam = v0$gam, rho = v0$rho, fallback = FALSE)
  }
  .eu3_memo$key <- key; .eu3_memo$val <- val; val
}
register_candidate("erdsl_uc", function(counts, meta, formula, tested_term) {
  v <- .eu3_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p3, estimate = v$est)
}, notes = "it29: max(|z|) over detection, square-root and log pairwise kernels, trivariate normal")
register_candidate("ersl_uc", function(counts, meta, formula, tested_term) {
  v <- .eu3_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_ls, estimate = v$est)
}, notes = "it29: max over square-root and log kernels")
register_candidate("ers_uc", function(counts, meta, formula, tested_term) {
  v <- .eu3_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_sqrt, estimate = v$est)
}, notes = "it29: square-root kernel alone")

# --- it30: positive-part kernel (the two-part idea inside the exact pairwise frame) -----------------------
# msq (SimulateMSeq) keeps each template sample's zeros and multiplies its non-zero counts: the signal is
# entirely in magnitude-given-presence, and the log kernel dilutes it with the zeros (log1p(0) = 0 in both
# groups). Compare magnitudes only where BOTH members of a pair are present at their common depth:
#   h_ij = E[1{Y_i>0} 1{Y_j>0} (log Y_i - log Y_j)] = f_j L_i - f_i L_j,   L = E[log Y_D; Y_D > 0], f = ERD,
# over independent thinnings of the two samples. Symmetric weight x antisymmetric difference, so under H0
# (the two thinned counts identically distributed) its mean is exactly zero, like every other kernel here;
# a structural zero (y = 0) removes the pair. Detection carries the presence part, as before.
.pos_kernel <- function(Yi, Ni, Di, Yj, Nj, Dj) {
  fi <- .erd_f_var(Yi, Ni, Di); fj <- .erd_f_var(Yj, Nj, Dj)
  Li <- .erl_f_var(Yi, Ni, Di, kmax = 70L, big_mu = 20, off = 0); Lj <- .erl_f_var(Yj, Nj, Dj, kmax = 70L, big_mu = 20, off = 0)
  fj * Li - fi * Lj
}
.eup_memo <- new.env()
.eupos_fit <- function(counts, meta, formula, tested_term) {
  key <- list(counts, meta, formula, tested_term)
  if (!is.null(.eup_memo$key) && identical(.eup_memo$key, key)) return(.eup_memo$val)
  d <- .uc_design(meta, formula, tested_term); v0 <- .euc_fit(counts, meta, formula, tested_term)
  if (is.null(d)) { val <- v0; val$p_pos <- NA; val$p_dp <- v0$p_max; val$p_dlp <- v0$p_max
  } else {
    depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts); g1 <- d$g1; Z <- d$Z; n1 <- sum(g1)
    ud <- .u_stat_cov(counts, depth, g1, Z, v0$gam, "det", rho = v0$rho); ul <- .u_stat_cov(counts, depth, g1, Z, v0$gam, "log", rho = v0$rho)
    up <- .u_stat_cov(counts, depth, g1, Z, v0$gam, "pos", rho = v0$rho)
    ok <- rowSums(counts > 0) >= 3 & ud$v > 0 & ul$v > 0
    okp <- ok & up$v > 0
    tz <- function(u, o) ifelse(o, stats::qnorm(stats::pt(u$U / sqrt(u$v), u$df)), NA)
    zd <- tz(ud, ok); zl <- tz(ul, ok); zp <- tz(up, okp)
    cr <- function(u1, u2) ifelse(okp, .cov_uu(u1, u2, n1) / sqrt(u1$v * u2$v), NA)
    rdl <- ifelse(ok, .cov_uu(ud, ul, n1) / sqrt(ud$v * ul$v), NA); rdp <- cr(ud, up); rlp <- cr(ul, up)
    pdp <- .pmax2(pmax(abs(zd), abs(zp)), rdp); pdlp <- .pmax3(pmax(abs(zd), abs(zl), abs(zp)), rdl, rdp, rlp)
    fb <- !okp & ok; pdp[fb] <- v0$p_max[fb]; pdlp[fb] <- v0$p_max[fb]           # no pair with both present: detection/log only
    val <- list(feature = rownames(counts), p_det = 2 * stats::pnorm(-abs(zd)), p_log = 2 * stats::pnorm(-abs(zl)), p_pos = 2 * stats::pnorm(-abs(zp)),
                p_max = v0$p_max, p_dp = pdp, p_dlp = pdlp, est = v0$est, gam = v0$gam, rho = v0$rho, fallback = FALSE)
  }
  .eup_memo$key <- key; .eup_memo$val <- val; val
}
register_candidate("erdp_uc", function(counts, meta, formula, tested_term) {
  v <- .eupos_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_dp, estimate = v$est)
}, notes = "it30: max(detection, positive-part log) pairwise kernels -- a two-part test inside the exact pairwise frame")
register_candidate("erdlp_uc", function(counts, meta, formula, tested_term) {
  v <- .eupos_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_dlp, estimate = v$est)
}, notes = "it30: max(detection, log, positive-part log), trivariate normal")

# --- it31: per-taxon, label-free choice of scale (log vs square root) --------------------------------------
# srv9/local: on msq the square-root scale reaches LDM's power (t on arcsin-root 15/13 TP vs log 11/10,
# LDM 16/16) while on house the log scale wins (R08 26 vs 19); a max over both costs ~2x in p and loses on
# both. Which scale is efficient is a property of the taxon's distribution, not of the labels: for a
# multiplicative change in abundance, the efficiency of a transform g is (E[y g'(y)])^2 / Var(g(y)). It is
# estimated per taxon on the POOLED samples (counts scaled to the median depth), so the choice is invariant
# to relabelling and the test given the choice keeps its null distribution. The chosen kernel then enters
# the lead's max test with detection.
.scale_choice <- function(counts, depth) {
  Y <- as.matrix(counts) * rep(stats::median(depth) / depth, each = nrow(counts))
  effl <- rowMeans(Y / (1 + Y))^2 / apply(log1p(Y), 1, stats::var)
  effs <- rowMeans(0.5 * sqrt(Y))^2 / apply(sqrt(Y), 1, stats::var)
  ifelse(is.finite(effs) & is.finite(effl) & effs > effl, "sqrt", "log")
}
.euch_memo <- new.env()
.euch_fit <- function(counts, meta, formula, tested_term) {
  key <- list(counts, meta, formula, tested_term)
  if (!is.null(.euch_memo$key) && identical(.euch_memo$key, key)) return(.euch_memo$val)
  d <- .uc_design(meta, formula, tested_term); v0 <- .euc_fit(counts, meta, formula, tested_term)
  if (is.null(d)) { val <- v0; val$p_ch <- v0$p_max
  } else {
    depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts); g1 <- d$g1; Z <- d$Z; n1 <- sum(g1)
    ch <- .scale_choice(counts, depth)
    ud <- .u_stat_cov(counts, depth, g1, Z, v0$gam, "det", rho = v0$rho); ul <- .u_stat_cov(counts, depth, g1, Z, v0$gam, "log", rho = v0$rho)
    us <- .u_stat_cov(counts, depth, g1, Z, v0$gam, "sqrt", rho = v0$rho)
    sel <- ch == "sqrt"; um <- ul; for (nm in c("U", "v", "df")) um[[nm]][sel] <- us[[nm]][sel]; um$a[sel, ] <- us$a[sel, ]; um$b[sel, ] <- us$b[sel, ]
    ok <- rowSums(counts > 0) >= 3 & ud$v > 0 & um$v > 0
    tz <- function(u) ifelse(ok, stats::qnorm(stats::pt(u$U / sqrt(u$v), u$df)), NA)
    zd <- tz(ud); zm <- tz(um); rc <- ifelse(ok, pmin(pmax(.cov_uu(ud, um, n1) / sqrt(ud$v * um$v), -0.999), 0.999), NA)
    val <- list(feature = rownames(counts), p_det = 2 * stats::pnorm(-abs(zd)), p_mag = 2 * stats::pnorm(-abs(zm)),
                p_ch = .pmax2(pmax(abs(zd), abs(zm)), rc), p_max = v0$p_max, choice = ch, est = v0$est, gam = v0$gam, rho = v0$rho, fallback = FALSE)
  }
  .euch_memo$key <- key; .euch_memo$val <- val; val
}
register_candidate("erdch_uc", function(counts, meta, formula, tested_term) {
  v <- .euch_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_ch, estimate = v$est)
}, notes = "it31: max(detection, magnitude kernel) with the magnitude scale (log / sqrt) chosen per taxon, label-free, by estimated efficiency")

# E[Y_D^lambda], 0 < lambda < 1 (lambda = 0.25: between the log and the square root)
.erp_f_var <- function(counts, depth, Dvec, lambda = 0.25, kmax = 70L, big_mu = 20) {
  Y <- as.matrix(counts); N <- matrix(depth, nrow(Y), ncol(Y), byrow = TRUE); Dm <- matrix(Dvec, nrow(Y), ncol(Y), byrow = TRUE)
  mu <- Dm * Y / N; out <- matrix(0, nrow(Y), ncol(Y)); L <- lambda
  big <- mu > big_mu
  if (any(big)) { pp <- (Y / N)[big]; Nb <- N[big]; Db <- Dm[big]; fpc <- (Nb - Db) / pmax(Nb - 1, 1)
    v <- Db * pp * (1 - pp) * fpc; k3 <- v * (1 - 2 * pp) * (Nb - 2 * Db) / pmax(Nb - 2, 1); m1 <- mu[big]
    out[big] <- m1^L + L * (L - 1) * m1^(L - 2) * v / 2 + L * (L - 1) * (L - 2) * m1^(L - 3) * k3 / 6 +
      L * (L - 1) * (L - 2) * (L - 3) * m1^(L - 4) * (3 * v^2 + v * (1 - 6 * pp * (1 - pp))) / 24 }
  sm <- which(!big & Y > 0)
  if (length(sm)) { y <- Y[sm]; nn <- N[sm] - Y[sm]; dd <- Dm[sm]; acc <- numeric(length(sm))
    for (k in 1:kmax) { live <- k <= y & k <= dd; if (!any(live)) break
      acc[live] <- acc[live] + stats::dhyper(k, y[live], nn[live], dd[live]) * k^L }
    out[sm] <- acc }
  out
}

.euq4_memo <- new.env()
.euq4_fit <- function(counts, meta, formula, tested_term) {
  key <- list(counts, meta, formula, tested_term)
  if (!is.null(.euq4_memo$key) && identical(.euq4_memo$key, key)) return(.euq4_memo$val)
  d <- .uc_design(meta, formula, tested_term); v0 <- .euc_fit(counts, meta, formula, tested_term)
  if (is.null(d)) { val <- v0; val$p_dq <- v0$p_max
  } else {
    depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts); g1 <- d$g1; n1 <- sum(g1)
    ud <- .u_stat_cov(counts, depth, g1, d$Z, v0$gam, "det", rho = v0$rho); uq <- .u_stat_cov(counts, depth, g1, d$Z, v0$gam, "q4", rho = v0$rho)
    ok <- rowSums(counts > 0) >= 3 & ud$v > 0 & uq$v > 0
    tz <- function(u) ifelse(ok, stats::qnorm(stats::pt(u$U / sqrt(u$v), u$df)), NA)
    zd <- tz(ud); zq <- tz(uq); rc <- ifelse(ok, pmin(pmax(.cov_uu(ud, uq, n1) / sqrt(ud$v * uq$v), -0.999), 0.999), NA)
    val <- list(feature = rownames(counts), p_q4 = 2 * stats::pnorm(-abs(zq)), p_dq = .pmax2(pmax(abs(zd), abs(zq)), rc), est = v0$est, gam = v0$gam, rho = v0$rho, fallback = FALSE)
  }
  .euq4_memo$key <- key; .euq4_memo$val <- val; val
}
register_candidate("erdq_uc", function(counts, meta, formula, tested_term) {
  v <- .euq4_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_dq, estimate = v$est)
}, notes = "it31b: max(detection, fourth-root kernel)")

# --- it32: label-free per-taxon kernel choice with the Poisson sensitivity identity -----------------------
# it31's criterion treated a fold change as scaling the observed count. Under Poisson sampling of a latent
# mean mu, d E[g(Y)] / d log mu = E[mu (g(Y+1) - g(Y))] = E[Y (g(Y) - g(Y-1))], estimable from the pooled
# counts without labels. Efficiency of kernel g for a multiplicative change: S_g^2 / Var g(Y), with
#   S_det = P(Y = 1),  S_log = E[Y log((1+Y)/Y)],  S_sqrt = E[Y (sqrt(Y) - sqrt(Y-1))].
# Pick the most efficient kernel per taxon (single kernel, no max penalty), or the more efficient of
# log/sqrt inside the lead's max with detection. The choice uses pooled counts only.
.kernel_eff <- function(counts, depth) {
  Y <- round(as.matrix(counts) * rep(stats::median(depth) / depth, each = nrow(counts)))
  vr <- function(M) apply(M, 1, stats::var)
  Sd <- rowMeans(Y == 1); Sl <- rowMeans(ifelse(Y > 0, Y * log((1 + Y) / pmax(Y, 1)), 0)); Ss <- rowMeans(ifelse(Y > 0, Y * (sqrt(Y) - sqrt(pmax(Y - 1, 0))), 0))
  cbind(det = Sd^2 / vr(Y > 0), log = Sl^2 / vr(log1p(Y)), sqrt = Ss^2 / vr(sqrt(Y)))
}
.euk3_memo <- new.env()
.euk3_fit <- function(counts, meta, formula, tested_term) {
  key <- list(counts, meta, formula, tested_term)
  if (!is.null(.euk3_memo$key) && identical(.euk3_memo$key, key)) return(.euk3_memo$val)
  d <- .uc_design(meta, formula, tested_term); v0 <- .euc_fit(counts, meta, formula, tested_term)
  if (is.null(d)) { val <- v0; val$p_one <- v0$p_max; val$p_dmag <- v0$p_max
  } else {
    depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts); g1 <- d$g1; n1 <- sum(g1)
    E <- .kernel_eff(counts, depth); E[!is.finite(E)] <- 0
    U <- lapply(c(det = "det", log = "log", sqrt = "sqrt"), function(h) .u_stat_cov(counts, depth, g1, d$Z, v0$gam, h, rho = v0$rho))
    ok <- rowSums(counts > 0) >= 3 & U$det$v > 0 & U$log$v > 0 & U$sqrt$v > 0
    z <- lapply(U, function(u) ifelse(ok, stats::qnorm(stats::pt(u$U / sqrt(u$v), u$df)), NA))
    best <- colnames(E)[max.col(E, ties.method = "first")]
    zone <- ifelse(best == "det", z$det, ifelse(best == "log", z$log, z$sqrt))
    mag <- ifelse(E[, "sqrt"] > E[, "log"], "sqrt", "log"); zm <- ifelse(mag == "sqrt", z$sqrt, z$log)
    rdl <- .cov_uu(U$det, U$log, n1) / sqrt(U$det$v * U$log$v); rds <- .cov_uu(U$det, U$sqrt, n1) / sqrt(U$det$v * U$sqrt$v)
    rc <- ifelse(ok, pmin(pmax(ifelse(mag == "sqrt", rds, rdl), -0.999), 0.999), NA)
    val <- list(feature = rownames(counts), p_one = 2 * stats::pnorm(-abs(zone)), p_dmag = .pmax2(pmax(abs(z$det), abs(zm)), rc),
                choice = best, mag = mag, p_max = v0$p_max, est = v0$est, gam = v0$gam, rho = v0$rho, fallback = FALSE)
  }
  .euk3_memo$key <- key; .euk3_memo$val <- val; val
}
register_candidate("eone_uc", function(counts, meta, formula, tested_term) {
  v <- .euk3_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_one, estimate = v$est)
}, notes = "it32: one pairwise kernel per taxon (detection / log / sqrt), chosen label-free by Poisson-identity efficiency")
register_candidate("erdmag_uc", function(counts, meta, formula, tested_term) {
  v <- .euk3_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_dmag, estimate = v$est)
}, notes = "it32: max(detection, log-or-sqrt chosen label-free per taxon)")

# --- it33: kernel selection across taxa (leave-one-out) ---------------------------------------------------
# it31/it32: which scale is efficient depends on the TYPE of alternative (prevalence + abundance on house,
# abundance-only on msq), which pooled counts cannot reveal; max over kernels pays ~2x in p. But the
# other taxa can: the signal of the whole table shows on which scale the exposure acts. For taxon j,
# count over the OTHER taxa how many reach p < 0.001 with each option -- the lead's max(det, log) and the
# square-root kernel -- and test j with the option that wins (ties -> the lead). Taxon j's own statistic
# never enters its choice, so its p-value is used as is; each option is valid on its own.
.eusel_memo <- new.env()
.eusel_fit <- function(counts, meta, formula, tested_term, thr = 1e-3, margin = 2L) {
  key <- list(counts, meta, formula, tested_term, thr, margin)
  if (!is.null(.eusel_memo$key) && identical(.eusel_memo$key, key)) return(.eusel_memo$val)
  v <- .eu3_fit(counts, meta, formula, tested_term)
  if (isTRUE(v$fallback) || is.null(v$p_sqrt) || all(!is.finite(v$p_sqrt))) { val <- v; val$p_sel <- v$p_max; val$sel_sqrt <- 0
  } else {
    P <- cbind(lead = v$p_max, sqrt = v$p_sqrt); hit <- !is.na(P) & P < thr
    tot <- colSums(hit); loo <- sweep(-hit, 2, tot, "+")                          # counts over the other taxa
    use_sqrt <- loo[, "sqrt"] >= loo[, "lead"] + margin                             # leave-one-out with a margin: near ties keep the lead
    val <- v; val$p_sel <- ifelse(use_sqrt, v$p_sqrt, v$p_max); val$sel_sqrt <- mean(use_sqrt); val$hits <- tot
  }
  .eusel_memo$key <- key; .eusel_memo$val <- val; val
}
register_candidate("esel_uc", function(counts, meta, formula, tested_term) {
  v <- .eusel_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_sel, estimate = v$est)
}, notes = "it33: per taxon, the lead's max(det, log) or the sqrt kernel -- whichever gives more p < 0.001 among the OTHER taxa")

# it33 + it24: the lead option is erdl_um (moderated; pairwise clustered designs), the alternative the sqrt kernel
.euselm_memo <- new.env()
.euselm_fit <- function(counts, meta, formula, tested_term, thr = 1e-3, margin = 2L) {
  key <- list(counts, meta, formula, tested_term, thr, margin)
  if (!is.null(.euselm_memo$key) && identical(.euselm_memo$key, key)) return(.euselm_memo$val)
  vm <- .eum_fit(counts, meta, formula, tested_term); v3 <- .eu3_fit(counts, meta, formula, tested_term)
  val <- vm; val$p_sel <- vm$p_max; val$sel_sqrt <- 0
  if (!isTRUE(v3$fallback) && !is.null(v3$p_sqrt) && any(is.finite(v3$p_sqrt)) && is.null(.cluster_of(meta))) {
    P <- cbind(lead = vm$p_max, sqrt = v3$p_sqrt); hit <- !is.na(P) & P < thr
    loo <- sweep(-hit, 2, colSums(hit), "+"); use <- loo[, "sqrt"] >= loo[, "lead"] + margin
    val$p_sel <- ifelse(use, v3$p_sqrt, vm$p_max); val$sel_sqrt <- mean(use) }
  .euselm_memo$key <- key; .euselm_memo$val <- val; val
}
register_candidate("esel_um", function(counts, meta, formula, tested_term) {
  v <- .euselm_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_sel, estimate = v$est)
}, notes = "it33 + it24: erdl_um or the sqrt kernel, chosen per taxon by leave-one-out hit counts over the other taxa (margin 2)")

# --- it34: per-sample compositional size factors inside the thinning ---------------------------------------
# Every comparator that beats the lead somewhere normalises each SAMPLE by a stable reference -- ADAPT by
# its reference taxa, LinDA by the CLR, ZicoSeq by its reference set -- while the lead normalises by total
# reads (thinning to a common depth) plus one group-level factor gamma. When a few dominant taxa fluctuate
# from sample to sample, every other taxon's proportion fluctuates inversely, and that noise enters every
# pair comparison. Per-sample log size factor sf_i = median over common taxa (prevalence >= 0.8) of
# log((y_ij + 0.5) / N_i) - its taxon mean, centred; label-free. In the pair (i, j) both samples are thinned
# to a common EFFECTIVE depth: D_i = De / b_i, D_j = De / b_j with b = exp(gamma 1{exposed} + sf) and
# De = rho min(N_i b_i, N_j b_j) -- the group-level gamma generalised to every sample.
.size_factors <- function(counts, depth, prev_min = 0.8, nmin = 20L) {
  Y <- as.matrix(counts); pv <- rowMeans(Y > 0); R <- which(pv >= prev_min)
  if (length(R) < nmin) R <- order(pv, decreasing = TRUE)[seq_len(min(nmin, nrow(Y)))]
  L <- log((Y[R, , drop = FALSE] + 0.5) / rep(depth, each = length(R))); L <- L - rowMeans(L)
  sf <- apply(L, 2, stats::median); sf - mean(sf)
}
.eusf_memo <- new.env()
.eusf_fit <- function(counts, meta, formula, tested_term) {
  key <- list(counts, meta, formula, tested_term)
  if (!is.null(.eusf_memo$key) && identical(.eusf_memo$key, key)) return(.eusf_memo$val)
  d <- .uc_design(meta, formula, tested_term)
  if (is.null(d)) { val <- .euc_fit(counts, meta, formula, tested_term)
  } else {
    depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts); g1 <- d$g1; Z <- d$Z; n1 <- sum(g1)
    rho <- round(.rho_design(depth, g1), 2); sf <- .size_factors(counts, depth)
    set.seed(7L); rows <- if (nrow(counts) > 300L) sort(sample.int(nrow(counts), 300L)) else seq_len(nrow(counts))
    ct <- counts[rows, , drop = FALSE]; ct <- ct[rowSums(ct > 0) >= 3, , drop = FALSE]
    med <- function(gm) { u <- .u_stat_cov(ct, depth, g1, Z, gm, "log", rho = rho, sf = sf); stats::median(u$U / sqrt(u$v), na.rm = TRUE) }
    flo <- med(-2); fhi <- med(2)
    gm <- if (is.finite(flo) && is.finite(fhi) && sign(flo) != sign(fhi)) stats::uniroot(med, c(-2, 2), f.lower = flo, f.upper = fhi, tol = 2e-3)$root else 0
    ud <- .u_stat_cov(counts, depth, g1, Z, gm, "det", rho = rho, sf = sf); ul <- .u_stat_cov(counts, depth, g1, Z, gm, "log", rho = rho, sf = sf)
    ok <- rowSums(counts > 0) >= 3 & ud$v > 0 & ul$v > 0
    tz <- function(u) ifelse(ok, stats::qnorm(stats::pt(u$U / sqrt(u$v), u$df)), NA)
    z1 <- tz(ud); z2 <- tz(ul); rc <- ifelse(ok, .cov_uu(ud, ul, n1) / sqrt(ud$v * ul$v), NA)
    val <- list(feature = rownames(counts), p_det = 2 * stats::pnorm(-abs(z1)), p_log = 2 * stats::pnorm(-abs(z2)),
                p_max = .pmax2(pmax(abs(z1), abs(z2)), rc), est = ul$U / log(2), gam = gm, rho = rho, sf_sd = stats::sd(sf), fallback = FALSE)
  }
  .eusf_memo$key <- key; .eusf_memo$val <- val; val
}
register_candidate("erdl_usf", function(counts, meta, formula, tested_term) {
  v <- .eusf_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_max, estimate = v$est)
}, notes = "it34: erdl_uc with per-sample compositional size factors in the common-depth thinning")

# --- it35: the combination -- size factors (it34) + moderation and pairwise clusters (it24) + scale chosen
# across taxa (it33). Each piece is label-free or leave-one-out, so the composite keeps the lead's null.
.efull_memo <- new.env()
.efull_fit <- function(counts, meta, formula, tested_term, thr = 1e-3, margin = 2L, select = TRUE, sf_balanced_only = FALSE, winsor = 0, opts3 = FALSE) {
  key <- list(counts, meta, formula, tested_term, thr, margin, select, sf_balanced_only, winsor, opts3)
  if (!is.null(.efull_memo$key) && identical(.efull_memo$key, key)) return(.efull_memo$val)
  d <- .uk_design(meta, formula, tested_term)
  if (is.null(d)) { val <- .euc_fit(counts, meta, formula, tested_term); val$p_full <- val$p_max; val$p_ms <- val$p_max
  } else {
    depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts); g1 <- d$g1; Z <- d$Z; cl <- d$cl; n1 <- sum(g1)
    rho <- round(.rho_design(depth, g1), 2); sf <- .size_factors(counts, depth)
    if (sf_balanced_only && rho < 1) sf <- NULL                                    # size factors only when depth is balanced
    if (winsor > 0 && rho == 1) counts <- .winsorize_counts(counts, depth, winsor) # it36: cap each taxon's top proportions (label-free)
    set.seed(7L); rows <- if (nrow(counts) > 300L) sort(sample.int(nrow(counts), 300L)) else seq_len(nrow(counts))
    ct <- counts[rows, , drop = FALSE]; ct <- ct[rowSums(ct > 0) >= 3, , drop = FALSE]
    med <- function(gm) { u <- .u_stat_cov(ct, depth, g1, Z, gm, "log", rho = rho, cl = cl, sf = sf); stats::median(u$U / sqrt(u$v), na.rm = TRUE) }
    flo <- med(-2); fhi <- med(2)
    gm <- if (is.finite(flo) && is.finite(fhi) && sign(flo) != sign(fhi)) stats::uniroot(med, c(-2, 2), f.lower = flo, f.upper = fhi, tol = 2e-3)$root else 0
    U <- lapply(c(det = "det", log = "log", sqrt = "sqrt"), function(h) .u_stat_cov(counts, depth, g1, Z, gm, h, rho = rho, cl = cl, sf = sf))
    ok <- rowSums(counts > 0) >= 3 & U$det$v > 0 & U$log$v > 0 & U$sqrt$v > 0
    xx <- stats::qlogis(pmin(pmax(rowMeans(counts > 0), 0.01), 0.99))
    mod <- function(u) if (is.null(cl)) .squeeze(ifelse(ok, u$v, NA), u$df, x = xx) else .squeeze(ifelse(ok, u$v, NA), u$df, off = u$vind, dfmax = u$dfind)
    M <- lapply(U, mod)
    z <- mapply(function(u, m) ifelse(ok, stats::qnorm(stats::pt(u$U / sqrt(m$v), m$df)), NA), U, M, SIMPLIFY = FALSE)
    cv <- if (is.null(cl)) .cov_uu(U$det, U$log, n1) else {
      Ma <- stats::model.matrix(~ U$det$ca - 1); Mb <- stats::model.matrix(~ U$det$cb - 1); Ga <- ncol(Ma); Gb <- ncol(Mb); n0 <- ncol(U$det$b)
      rowSums((U$det$a %*% Ma) * (U$log$a %*% Ma)) / n1^2 * Ga / (Ga - 1) + rowSums((U$det$b %*% Mb) * (U$log$b %*% Mb)) / n0^2 * Gb / (Gb - 1) }
    rc <- ifelse(ok, pmin(pmax(cv / sqrt(U$det$v * U$log$v), -0.999), 0.999), NA)
    pl <- .pmax2(pmax(abs(z$det), abs(z$log)), rc); ps <- 2 * stats::pnorm(-abs(z$sqrt))
    use <- rep(FALSE, length(pl)); pfull <- pl
    if (select && !opts3) { P <- cbind(lead = pl, sqrt = ps); hit <- !is.na(P) & P < thr; loo <- sweep(-hit, 2, colSums(hit), "+"); use <- loo[, "sqrt"] >= loo[, "lead"] + margin; pfull <- ifelse(use, ps, pl) }
    if (select && opts3) {                                                         # it37: the lead, log alone or sqrt alone
      plg <- 2 * stats::pnorm(-abs(z$log)); P <- cbind(lead = pl, log = plg, sqrt = ps); hit <- !is.na(P) & P < thr; loo <- sweep(-hit, 2, colSums(hit), "+")
      alt <- ifelse(loo[, "sqrt"] >= loo[, "log"], "sqrt", "log"); altn <- pmax(loo[, "sqrt"], loo[, "log"])
      use <- altn >= loo[, "lead"] + margin; pfull <- ifelse(use, ifelse(alt == "sqrt", ps, plg), pl) }
    val <- list(feature = rownames(counts), p_full = pfull, p_ms = pl, p_det = 2 * stats::pnorm(-abs(z$det)), p_log = 2 * stats::pnorm(-abs(z$log)),
                p_sqrt = ps, sel_sqrt = mean(use), est = U$log$U / log(2), gam = gm, rho = rho, fallback = FALSE)
  }
  .efull_memo$key <- key; .efull_memo$val <- val; val
}
register_candidate("efull", function(counts, meta, formula, tested_term) {
  v <- .efull_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_full, estimate = v$est)
}, notes = "it35: size factors + moderation + pairwise clusters + scale chosen across taxa")
register_candidate("emsf", function(counts, meta, formula, tested_term) {
  v <- .efull_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_ms, estimate = v$est)
}, notes = "it35 without the scale selection: size factors + moderation + pairwise clusters")

register_candidate("efull_b", function(counts, meta, formula, tested_term) {
  v <- .efull_fit(counts, meta, formula, tested_term, sf_balanced_only = TRUE); data.frame(feature = v$feature, p = v$p_full, estimate = v$est)
}, notes = "it35b: efull with size factors only when depth is balanced across groups (rho = 1)")

# --- it36: winsorised proportions (ZicoSeq's default: top 3% per taxon) -----------------------------------
# p09: msq is where the lead loses most (calibrated TPR 0.075 vs ZicoSeq 0.105, LDM 0.073 at higher FDR;
# ~1.5-2x per regime). ZicoSeq as the benchmark runs it: square-root link on reference-normalised,
# posterior-sampled proportions WINSORISED at each taxon's top 3%. A few samples with extreme proportions
# dominate a mean-difference statistic on the sqrt (and partly the log) scale. Cap each count at its
# taxon's 97th-percentile proportion times the sample's depth (positives stay >= 1); label-free, only
# when depth is balanced (the cap scales with depth, so thinning stays approximately consistent).
.winsorize_counts <- function(counts, depth, pct = 0.03) {
  Y <- as.matrix(counts); R <- Y / rep(depth, each = nrow(Y))
  q <- apply(R, 1, stats::quantile, probs = 1 - pct, names = FALSE)
  cap <- ceiling(outer(q, depth)); Yw <- ifelse(Y > 0, pmax(pmin(Y, cap), 1), 0); dimnames(Yw) <- dimnames(Y); Yw
}
register_candidate("efull_bw", function(counts, meta, formula, tested_term) {
  v <- .efull_fit(counts, meta, formula, tested_term, sf_balanced_only = TRUE, winsor = 0.03); data.frame(feature = v$feature, p = v$p_full, estimate = v$est)
}, notes = "it36: efull_b on counts winsorised at each taxon's top 3% of proportions (balanced depth only)")

# --- it37: three options for the leave-one-out selection ------------------------------------------------
# On msq the detection half of max(det, log) only costs (zeros are inherited, no prevalence signal): local
# tongue cells efull_b 8 / 5 vs log alone 9 / 9. Let the other taxa also choose "log alone".
register_candidate("efull_b3", function(counts, meta, formula, tested_term) {
  v <- .efull_fit(counts, meta, formula, tested_term, sf_balanced_only = TRUE, opts3 = TRUE); data.frame(feature = v$feature, p = v$p_full, estimate = v$est)
}, notes = "it37: efull_b with the selection over max(det, log) / log alone / sqrt alone")
register_candidate("efull_b3w", function(counts, meta, formula, tested_term) {
  v <- .efull_fit(counts, meta, formula, tested_term, sf_balanced_only = TRUE, opts3 = TRUE, winsor = 0.03); data.frame(feature = v$feature, p = v$p_full, estimate = v$est)
}, notes = "it36 + it37: winsorised counts and three-option selection")
register_candidate("efull_bc", function(counts, meta, formula, tested_term) {
  v <- .efull_fit(counts, meta, formula, tested_term, sf_balanced_only = TRUE); data.frame(feature = v$feature, p = v$p_full, estimate = v$est)
}, notes = "efull_b (same as efull_b; registered for the dev suite comparison)")
