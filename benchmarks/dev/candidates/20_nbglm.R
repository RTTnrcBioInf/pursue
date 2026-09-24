# -----------------------------------------------------------------------------
# 20_nbglm.R -- iteration 2: NB GLM on ALL cells, depth offset, per-model dispersion, centred
#
#   X_ij ~ NB( mu_ij = N_i * exp(x_i' b_j),  theta_j )        (MASS::glm.nb, theta by ML per fit)
#
# Why: every candidate so far used only the detection pattern (logistic, sharedlink) or only the
# detected cells (0.2's abundance arm, which recovers ~0% of abundance-only signal). None used zeros
# and magnitudes in one model. The ZINB mu-test (claude/results-diagnosis-and-0.3-direction.md)
# failed on DISPERSION -- theta fixed at a value shrunk to a cross-feature trend, so the variance
# was understated and the LRT anti-conservative (FPR 0.10-0.13 on a global null) -- not on its
# offset. For counts, log-depth at slope 1 is the sampling physics (expected reads scale with
# depth), and it holds for MIDASim too. So: estimate theta per feature and re-estimate it under
# the null, which is what an LRT needs to be calibrated; centre the tested coefficient at the
# empirical-null centre, as 0.2 does, for the compositional shift that total-read depth carries.
# Also covers features detected in every sample, which no detection test can see.
#
# nbglm: LRT of b_tested = delta (empirical-null centre); nbglm_raw: against 0.
# -----------------------------------------------------------------------------
.nb_memo <- new.env()

.nb_one <- function(y, X, off, tc, fix = NULL) {
  # fix = NULL: full model. Otherwise the tested columns are held at `fix` via the offset.
  if (is.null(fix)) { Xf <- X; o <- off } else { Xf <- X[, -tc, drop = FALSE]; o <- off + drop(X[, tc, drop = FALSE] %*% fix) }
  f <- tryCatch(suppressWarnings(MASS::glm.nb(y ~ 0 + Xf + offset(o), control = stats::glm.control(maxit = 50))),
                error = function(e) NULL)
  if (is.null(f) || !is.finite(f$twologlik)) return(NULL)
  f
}

.nb_fit <- function(counts, meta, formula, tested_term) {
  key <- list(counts, meta, formula, tested_term)
  if (!is.null(.nb_memo$key) && identical(.nb_memo$key, key)) return(.nb_memo$val)
  X <- stats::model.matrix(formula, meta); asg <- attr(X, "assign")
  tc <- which(asg == which(attr(stats::terms(formula), "term.labels") == tested_term)); k <- length(tc)
  depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts)
  off <- log(depth); m <- nrow(counts)
  full <- vector("list", m); b <- matrix(NA_real_, m, k); se <- matrix(NA_real_, m, k); seh <- matrix(NA_real_, m, k)
  for (j in seq_len(m)) {
    y <- as.numeric(counts[j, ]); if (sum(y > 0) < 3) next
    f <- .nb_one(y, X, off, tc); if (is.null(f)) next
    full[[j]] <- f$twologlik
    cf <- stats::coef(f); V <- tryCatch(stats::vcov(f), error = function(e) NULL)
    b[j, ] <- cf[tc]; if (!is.null(V)) se[j, ] <- sqrt(pmax(diag(V)[tc], 0))
    # sandwich (HC3): a variance that does not trust the NB variance function -- the fix for
    # count-model LRTs being anti-conservative on realistic, misspecified counts
    Vh <- tryCatch(sandwich::vcovHC(f, type = "HC3"), error = function(e) NULL)
    if (!is.null(Vh)) seh[j, ] <- sqrt(pmax(diag(Vh)[tc], 0))
  }
  delta <- vapply(seq_len(k), function(c) {
    okc <- is.finite(b[, c]) & is.finite(se[, c]) & se[, c] < 5          # drop separation blow-ups
    ce <- tryCatch(get("enull_center", envir = asNamespace("PURSUE"))(ifelse(okc, b[, c], NA), ifelse(okc, se[, c], NA)),
                   error = function(e) NULL)
    if (is.null(ce)) stats::median(b[okc, c]) else ce$delta }, numeric(1))
  lrt <- function(val) vapply(seq_len(m), function(j) {
    if (is.null(full[[j]])) return(NA_real_)
    f0 <- .nb_one(as.numeric(counts[j, ]), X, off, tc, fix = val); if (is.null(f0)) return(NA_real_)
    stats::pchisq(max(full[[j]] - f0$twologlik, 0), k, lower.tail = FALSE) }, numeric(1))
  # HC3 Wald, centred at its own empirical-null centre (k = 1: z-test; k > 1 would need the full
  # sandwich covariance of the tested block -- not needed on the dev suite, flagged NA)
  delta_h <- if (k == 1) { okc <- is.finite(b[, 1]) & is.finite(seh[, 1]) & seh[, 1] < 5
    ce <- tryCatch(get("enull_center", envir = asNamespace("PURSUE"))(ifelse(okc, b[, 1], NA), ifelse(okc, seh[, 1], NA)), error = function(e) NULL)
    if (is.null(ce)) stats::median(b[okc, 1]) else ce$delta } else NA_real_
  p_hc <- if (k == 1) 2 * stats::pnorm(-abs((b[, 1] - delta_h) / seh[, 1])) else rep(NA_real_, m)
  val <- list(feature = rownames(counts), p = lrt(delta), p_raw = lrt(rep(0, k)), p_hc = p_hc,
              est = b[, 1] - delta[1], delta = delta, b = b[, 1], se = se[, 1], seh = seh[, 1])
  .nb_memo$key <- key; .nb_memo$val <- val; val
}

register_candidate("nbglm", function(counts, meta, formula, tested_term) {
  v <- .nb_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p, estimate = v$est)
}, notes = "it2: NB GLM on all cells, log-depth offset, theta by ML in each fit, effect tested against the empirical-null centre")

register_candidate("nbglm_raw", function(counts, meta, formula, tested_term) {
  v <- .nb_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_raw)
}, notes = "it2: as nbglm, tested against 0 -- isolates the centring")

register_candidate("nbglm_hc3", function(counts, meta, formula, tested_term) {
  v <- .nb_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_hc)
}, notes = "it2: NB GLM point estimates, sandwich HC3 Wald test against the empirical-null centre -- robust to a wrong variance function")
