# -----------------------------------------------------------------------------
# elementary.R -- baseline methods every benchmark must carry (Hawinkel 2019; Wirbel
# 2024; Pelto 2025). Uniform wrapper signature:
#   fn(counts [features x samples], meta, formula, tested_term, args = list())
#     -> data.frame(feature, arm, p, q, estimate, se, ci_lo, ci_hi, status)
# `arm` is "single" for one-part methods. Estimates are log2 fold changes where the
# method has one on that scale; NA otherwise.
# -----------------------------------------------------------------------------

.empty_result <- function(features, arm = "single", status = "failed") {
  data.frame(feature = features, arm = arm, p = NA_real_, q = NA_real_, estimate = NA_real_, se = NA_real_,
             ci_lo = NA_real_, ci_hi = NA_real_, status = status, stringsAsFactors = FALSE)
}
# A method's own multiple-testing procedure is PART OF THE METHOD (ALDEx2's expected-BH across
# Monte-Carlo instances, ZicoSeq's permutation FDR, LOCOM's, corncob's...). Overwriting it with
# plain BH benchmarks a method the authors never published. So: keep `q` when the wrapper
# supplies one, and fall back to BH only when the method offers nothing.
.finish <- function(df) {
  # `.st` lets a wrapper flag how it got its q (e.g. "bh_fallback_no_method_q") without
  # hard-coding the column order; folded into status so the smoke test reports it.
  if (".st" %in% names(df)) { i <- !is.na(df$.st) & is.na(df$status); df$status[i] <- df$.st[i]; df$.st <- NULL }
  if (!"q" %in% names(df) || all(is.na(df$q))) df$q <- stats::p.adjust(df$p, "BH")
  df$status[is.na(df$status)] <- "ok"; rownames(df) <- NULL; df
}
# Pull a q-vector out of a result object, trying several column names; NULL -> caller uses BH.
.q_named <- function(x, feats) { v <- .as_pvec(x, feats); if (is.null(v) || all(is.na(v))) NULL else v }
.tested_is_binary <- function(meta, tested_term) {
  v <- meta[[tested_term]]; !is.null(v) && (is.factor(v) || is.character(v) || length(unique(v)) == 2L) && length(unique(v)) == 2L
}
.logtss <- function(counts) { N <- colSums(counts); log(sweep(counts + 0.5, 2L, N, "/")) }

method_wilcoxon_tss <- function(counts, meta, formula, tested_term, args = list()) {
  feats <- rownames(counts)
  if (!.tested_is_binary(meta, tested_term) || length(all.vars(formula)) > 1L)
    return(.empty_result(feats, status = "not_applicable_needs_binary_no_covariates"))
  g <- factor(meta[[tested_term]]); R <- sweep(counts, 2L, colSums(counts), "/")
  p <- apply(R, 1L, function(r) tryCatch(stats::wilcox.test(r[g == levels(g)[2]], r[g == levels(g)[1]], exact = FALSE)$p.value, error = function(e) NA_real_))
  est <- apply(R, 1L, function(r) log2((mean(r[g == levels(g)[2]]) + 1e-12) / (mean(r[g == levels(g)[1]]) + 1e-12)))
  .finish(data.frame(feature = feats, arm = "single", p = p, q = NA, estimate = est, se = NA, ci_lo = NA, ci_hi = NA, status = NA, stringsAsFactors = FALSE))
}

method_lm_logtss <- function(counts, meta, formula, tested_term, args = list()) {
  feats <- rownames(counts); Y <- .logtss(counts)
  X <- stats::model.matrix(formula, meta); a <- attr(X, "assign"); labs <- attr(stats::terms(formula), "term.labels")
  tc <- which(a == which(labs == tested_term)); k <- length(tc)
  fit <- stats::lm.fit(X, t(Y)); df <- ncol(counts) - ncol(X)
  XtXi <- chol2inv(qr.R(fit$qr))[order(fit$qr$pivot), order(fit$qr$pivot)]
  s2 <- colSums(fit$residuals^2) / df
  B <- t(fit$coefficients)[, tc, drop = FALSE]
  if (k == 1L) { se <- sqrt(s2 * XtXi[tc, tc]); t <- B[, 1] / se; p <- 2 * stats::pt(-abs(t), df)
    est <- B[, 1] / log(2); se2 <- se / log(2); ci <- stats::qt(0.975, df) * se2
    return(.finish(data.frame(feature = feats, arm = "single", p = p, q = NA, estimate = est, se = se2, ci_lo = est - ci, ci_hi = est + ci, status = NA, stringsAsFactors = FALSE))) }
  V <- XtXi[tc, tc]; Fs <- sapply(seq_along(feats), function(j) as.numeric(t(B[j, ]) %*% solve(V, B[j, ])) / (k * s2[j]))
  p <- stats::pf(Fs, k, df, lower.tail = FALSE)
  .finish(data.frame(feature = feats, arm = "single", p = p, q = NA, estimate = NA, se = NA, ci_lo = NA, ci_hi = NA, status = NA, stringsAsFactors = FALSE))
}

method_limma_logtss <- function(counts, meta, formula, tested_term, args = list()) {
  feats <- rownames(counts); Y <- .logtss(counts)
  X <- stats::model.matrix(formula, meta); a <- attr(X, "assign"); labs <- attr(stats::terms(formula), "term.labels")
  tc <- which(a == which(labs == tested_term))
  fit <- limma::eBayes(limma::lmFit(Y, X), trend = TRUE, robust = TRUE)
  tt <- limma::topTable(fit, coef = tc, number = Inf, sort.by = "none")
  est <- if (length(tc) == 1L) tt$logFC / log(2) else NA
  se2 <- if (length(tc) == 1L) (fit$stdev.unscaled[, tc] * sqrt(fit$s2.post)) / log(2) else NA
  ci <- if (length(tc) == 1L) stats::qt(0.975, fit$df.total) * se2 else NA
  .finish(data.frame(feature = feats, arm = "single", p = tt$P.Value, q = NA, estimate = est, se = se2,
                     ci_lo = est - ci, ci_hi = est + ci, status = NA, stringsAsFactors = FALSE))
}

method_logistic_presence <- function(counts, meta, formula, tested_term, args = list()) {
  feats <- rownames(counts); D <- counts > 0
  X <- stats::model.matrix(formula, meta); a <- attr(X, "assign"); labs <- attr(stats::terms(formula), "term.labels")
  tc <- which(a == which(labs == tested_term)); k <- length(tc)
  res <- lapply(seq_along(feats), function(j) {
    d <- as.integer(D[j, ]); if (sum(d) < 3 || sum(d) > length(d) - 3) return(c(NA, NA, NA))
    f1 <- tryCatch(suppressWarnings(stats::glm.fit(X, d, family = stats::binomial())), error = function(e) NULL)
    f0 <- tryCatch(suppressWarnings(stats::glm.fit(X[, -tc, drop = FALSE], d, family = stats::binomial())), error = function(e) NULL)
    if (is.null(f1) || is.null(f0)) return(c(NA, NA, NA))
    lrt <- f0$deviance - f1$deviance
    c(stats::pchisq(max(lrt, 0), k, lower.tail = FALSE), if (k == 1L) f1$coefficients[tc] else NA, NA)
  })
  res <- do.call(rbind, res)
  .finish(data.frame(feature = feats, arm = "single", p = res[, 1], q = NA, estimate = res[, 2], se = NA, ci_lo = NA, ci_hi = NA,
                     status = ifelse(is.na(res[, 1]), "no_detection_variation", NA), stringsAsFactors = FALSE))
}

method_pursue <- function(counts, meta, formula, tested_term, args = list()) {
  feats <- rownames(counts)
  fit <- do.call(PURSUE::pursue, c(list(otu = t(counts), meta = meta, formula = formula, tested_term = tested_term,
                                        min_prevalence = 0, verbose = FALSE), args))
  r <- fit$results; r <- r[match(feats, r$feature), ]
  one_col <- "abund_lfc2" %in% names(r)
  ab <- data.frame(feature = feats, arm = "abundance", p = r$abund_p, q = r$abund_q,
                   estimate = if (one_col) r$abund_lfc2 else NA, se = if (one_col) r$abund_se2 else NA,
                   ci_lo = if (one_col) r$abund_ci2_lo else NA, ci_hi = if (one_col) r$abund_ci2_hi else NA,
                   status = r$abund_status, stringsAsFactors = FALSE)
  pr <- data.frame(feature = feats, arm = "presence", p = r$pres_p, q = r$pres_q,
                   estimate = if ("pres_logor" %in% names(r)) r$pres_logor else NA, se = if ("pres_se" %in% names(r)) r$pres_se else NA,
                   ci_lo = if ("pres_ci_lo" %in% names(r)) r$pres_ci_lo else NA, ci_hi = if ("pres_ci_hi" %in% names(r)) r$pres_ci_hi else NA,
                   status = r$pres_status, stringsAsFactors = FALSE)
  cb <- data.frame(feature = feats, arm = "combined", p = r$any_p, q = r$any_q, estimate = NA, se = NA, ci_lo = NA, ci_hi = NA, status = "ok", stringsAsFactors = FALSE)
  out <- rbind(.finish(ab), .finish(pr), .finish(cb))
  attr(out, "centre") <- fit$centre
  out
}
