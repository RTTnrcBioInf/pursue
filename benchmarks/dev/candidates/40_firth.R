# -----------------------------------------------------------------------------
# 40_firth.R -- iteration 3: Jeffreys-penalised (Firth) detection models, penalised LRT
#
# it1 diagnosis (diag_fp.R, 12 cells): the shared link's null z-scores are calibrated in the BULK
# (|z| median 0.67, 90% 1.63 vs N(0,1) 0.67 / 1.64) but its false positives are rare features
# (prevalence ~0.15 vs 0.29 for other nulls) with large realised detection shifts -- the tail of the
# chi-square approximation for sparse binary data, where likelihood-ratio and Wald tests are known to
# be anti-conservative. logistic_depth shows the same pattern, milder. The standard remedy for sparse
# logistic regression is Firth's penalty, 0.5 log det I(beta) (Jeffreys prior), with a penalised LRT
# whose penalty uses the full model's information at both fits (Heinze & Schemper 2002).
#
# One penalised fitter serves any link P = plogis(g(eta)) with g monotone: plain logistic has
# g = identity; the shared learned link has the spline g from 10_sharedlink.R. Weight
# w = P(1-P) g'^2 and dw/deta = P(1-P)(1-2P) g'^3 + 2 P(1-P) g' g''.
# -----------------------------------------------------------------------------
.lk_logit <- list(fP = function(x) pmin(pmax(stats::plogis(x), 1e-12), 1 - 1e-12),
                  fgp = function(x) rep(1, length(x)), fgpp = function(x) rep(0, length(x)))

.pen_feature <- function(d, off, X, L, start, fixed = integer(0), fixed_val = numeric(0)) {
  free <- setdiff(seq_len(ncol(X)), fixed)
  off2 <- off + if (length(fixed)) drop(X[, fixed, drop = FALSE] %*% fixed_val) else 0
  Xf <- X[, free, drop = FALSE]
  parts <- function(b) {
    e <- off2 + drop(Xf %*% b); P <- L$fP(e); g1 <- L$fgp(e); w <- P * (1 - P) * g1^2
    I <- crossprod(X * sqrt(w)); ch <- tryCatch(chol(I), error = function(e) NULL)
    list(e = e, P = P, g1 = g1, w = w, I = I, ch = ch)
  }
  obj <- function(b) { z <- parts(b); if (is.null(z$ch)) return(1e10)
    -(sum(d * log(z$P) + (1 - d) * log(1 - z$P)) + sum(log(diag(z$ch)))) }       # 0.5 log det I = sum log diag chol
  grad <- function(b) { z <- parts(b); if (is.null(z$ch)) return(rep(0, length(b)))
    Ii <- chol2inv(z$ch); qd <- rowSums((X %*% Ii) * X)
    dw <- z$P * (1 - z$P) * ((1 - 2 * z$P) * z$g1^3 + 2 * z$g1 * L$fgpp(z$e))
    -drop(crossprod(Xf, (d - z$P) * z$g1 + 0.5 * dw * qd)) }
  f <- tryCatch(stats::nlminb(start[free], obj, grad), error = function(e) NULL)
  if (is.null(f) || !is.finite(f$objective) || f$objective >= 1e10) return(NULL)
  b <- numeric(ncol(X)); b[free] <- f$par; if (length(fixed)) b[fixed] <- fixed_val
  z <- parts(f$par); V <- if (is.null(z$ch)) NULL else chol2inv(z$ch)
  list(b = b, pl = -f$objective, V = V)
}

# plain logistic + centred log-depth, Firth-penalised, penalised LRT against 0
register_candidate("firth_depth", function(counts, meta, formula, tested_term) {
  X <- stats::model.matrix(formula, meta); a <- attr(X, "assign")
  tc <- which(a == which(attr(stats::terms(formula), "term.labels") == tested_term))
  ld <- log(meta$depth); X <- cbind(X, ld = ld - mean(ld)); D <- counts > 0; n <- ncol(D)
  p <- vapply(seq_len(nrow(D)), function(j) {
    d <- as.numeric(D[j, ]); if (sum(d) < 3 || sum(d) > n - 3) return(NA_real_)
    st <- c(stats::qlogis(mean(d)), rep(0, ncol(X) - 1L))
    f1 <- .pen_feature(d, 0, X, .lk_logit, st); if (is.null(f1)) return(NA_real_)
    f0 <- .pen_feature(d, 0, X, .lk_logit, f1$b, fixed = tc, fixed_val = rep(0, length(tc))); if (is.null(f0)) return(NA_real_)
    stats::pchisq(max(2 * (f1$pl - f0$pl), 0), length(tc), lower.tail = FALSE) }, numeric(1))
  data.frame(feature = rownames(counts), p = p)
}, notes = "it3: logistic_depth with Firth's penalty and penalised LRT (sparse-feature tail)")

# the shared learned link, Firth-penalised per feature; robust empirical-null centre; penalised LRT
.slf_memo <- new.env()
.slf_fit <- function(counts, meta, formula, tested_term) {
  key <- list(counts, meta, formula, tested_term)
  if (!is.null(.slf_memo$key) && identical(.slf_memo$key, key)) return(.slf_memo$val)
  v <- .sl_fit(counts, meta, formula, tested_term)            # link learned as in it1 (memoised)
  X <- v$X; tc <- v$tc; k <- length(tc); L <- v$L; m <- nrow(counts)
  lN <- log(if (!is.null(meta$depth)) meta$depth else colSums(counts)); D <- counts > 0
  B <- v$B; pl <- rep(NA_real_, m); se <- matrix(NA_real_, m, k)
  for (j in v$ok) { f <- .pen_feature(as.numeric(D[j, ]), lN, X, L, B[j, ]); if (is.null(f)) next
    B[j, ] <- f$b; pl[j] <- f$pl; if (!is.null(f$V)) se[j, ] <- sqrt(pmax(diag(f$V)[tc], 0)) }
  delta <- vapply(seq_len(k), function(c) enull_robust(B[, tc[c]], se[, c])$delta, numeric(1))
  lrt <- function(val) vapply(seq_len(m), function(j) {
    if (!is.finite(pl[j])) return(NA_real_)
    r <- .pen_feature(as.numeric(D[j, ]), lN, X, L, B[j, ], fixed = tc, fixed_val = val); if (is.null(r)) return(NA_real_)
    stats::pchisq(max(2 * (pl[j] - r$pl), 0), k, lower.tail = FALSE) }, numeric(1))
  val <- list(feature = rownames(counts), p = lrt(delta), p_raw = lrt(rep(0, k)), est = B[, tc[1]] - delta[1], b = B[, tc[1]], se = se[, 1])
  .slf_memo$key <- key; .slf_memo$val <- val; val
}
register_candidate("sharedlink_firth", function(counts, meta, formula, tested_term) {
  v <- .slf_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p, estimate = v$est)
}, notes = "it3: shared learned link + Firth penalty per feature + robust centre, penalised LRT")
register_candidate("sharedlink_firth_raw", function(counts, meta, formula, tested_term) {
  v <- .slf_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_raw)
}, notes = "it3: sharedlink_firth tested against 0 (no centring)")

# --- Firth only where the data are sparse ------------------------------------------------------------
# it3: Firth fixed the shared link's tail (FDR 0.133 -> 0.047) but cost 15-25% of TP where nothing was
# wrong. diag_fp (12 cells): 26 of the shared link's 36 false positives have a minor-class count
# (min(#detected, #not detected)) below 20, against 15% of its true positives. The chi-square
# approximation fails where the minor class is small, so penalise there and nowhere else. The rule
# uses only the pooled detection counts, never the labels.
.SPARSE_K <- 20L
register_candidate("sharedlink_hyb", function(counts, meta, formula, tested_term) {
  v <- .sl_fit(counts, meta, formula, tested_term); w <- .slf_fit(counts, meta, formula, tested_term)
  D <- counts > 0; minor <- pmin(rowSums(D), ncol(D) - rowSums(D))
  p <- ifelse(minor < .SPARSE_K, w$p, v$p_rc)
  data.frame(feature = v$feature, p = p, estimate = ifelse(minor < .SPARSE_K, w$est, v$b - v$delta_rc[1]))
}, notes = "it6: shared link, robust centre; Firth-penalised LRT only for features with minor-class count < 20")
