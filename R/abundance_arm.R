# -----------------------------------------------------------------------------
# abundance_arm.R -- log relative abundance over detected cells
#
# Per feature: y = log(X / N) on rows with X > 0, regressed on the model design plus
# a centred log-depth covariate (the depth-truncation correction), optional
# winsorization, empirical-Bayes variance moderation across features (limma), and
# either moderated t / F statistics or HC3 sandwich Wald statistics.
# Effect estimates are then centred across features (see centering.R).
# -----------------------------------------------------------------------------

#' Fit the abundance arm for every feature
#' @keywords internal
fit_abundance_arm <- function(otu, depth, design, depth_adjust = TRUE,
                              min_nonzero = 8L, winsorize = TRUE, winsor_quantile = 0.03,
                              winsor_min_n = 20L, robust_se = FALSE, n_cores = 1L) {
  m <- ncol(otu); n <- nrow(otu)
  X_full <- if (depth_adjust) cbind(design$X, log_depth = design$log_depth) else design$X
  tcols <- design$tested_cols; k <- length(tcols); p <- ncol(X_full)
  cn <- design$coef_names[tcols]

  fit_one <- function(j) {
    xj <- otu[, j]; det <- xj > 0; n_obs <- sum(det)
    out <- list(n_obs = n_obs, status = "ok", s2 = NA_real_, df = NA_real_,
                beta = rep(NA_real_, k), V = matrix(NA_real_, k, k), se_hc3 = rep(NA_real_, k),
                V_hc3 = matrix(NA_real_, k, k), amean = NA_real_, n_wins = 0L)
    if (n_obs < min_nonzero) { out$status <- "too_few_detected"; return(out) }
    y <- log(xj[det] / depth[det]); out$amean <- mean(y)
    if (winsorize && n_obs >= winsor_min_n) {
      qs <- stats::quantile(y, c(winsor_quantile, 1 - winsor_quantile), names = FALSE)
      out$n_wins <- sum(y < qs[1] | y > qs[2]); y <- pmin(pmax(y, qs[1]), qs[2])
    }
    Xd <- X_full[det, , drop = FALSE]
    qrX <- qr(Xd)
    if (qrX$rank < p) { out$status <- "design_rank_deficient"; return(out) }
    df <- n_obs - p
    if (df < 2L) { out$status <- "df_too_small"; return(out) }
    coef <- qr.coef(qrX, y); res <- qr.resid(qrX, y)
    XtXi <- chol2inv(qr.R(qrX))[order(qrX$pivot), order(qrX$pivot), drop = FALSE]
    out$s2 <- sum(res^2) / df; out$df <- df
    out$beta <- unname(coef[tcols]); out$V <- XtXi[tcols, tcols, drop = FALSE]
    if (robust_se) {
      hc <- tryCatch(sandwich_hc3(Xd, res, XtXi), error = function(e) NULL)
      if (!is.null(hc)) { out$V_hc3 <- hc[tcols, tcols, drop = FALSE]; out$se_hc3 <- sqrt(diag(out$V_hc3)) }
    }
    out
  }
  fits <- if (n_cores > 1L) parallel::mclapply(seq_len(m), fit_one, mc.cores = n_cores) else lapply(seq_len(m), fit_one)

  status <- vapply(fits, function(f) f$status, character(1))
  n_obs  <- vapply(fits, function(f) f$n_obs, numeric(1))
  ok     <- status == "ok"
  s2 <- vapply(fits, function(f) f$s2, numeric(1)); dfv <- vapply(fits, function(f) f$df, numeric(1))
  amean <- vapply(fits, function(f) f$amean, numeric(1))
  beta <- t(vapply(fits, function(f) f$beta, numeric(k))); if (k == 1L) beta <- matrix(beta, ncol = 1L)
  colnames(beta) <- cn

  # empirical-Bayes variance moderation (limma), with abundance trend and robust prior
  d0 <- NA_real_; var_post <- s2; df_mod <- dfv
  if (sum(ok) >= 3L && !robust_se) {
    sq <- limma::squeezeVar(s2[ok], df = dfv[ok], covariate = amean[ok], robust = TRUE)
    var_post[ok] <- sq$var.post; d0 <- if (length(sq$df.prior) == 1L) sq$df.prior else stats::median(sq$df.prior)
    df_mod[ok] <- pmin(dfv[ok] + sq$df.prior, sum(dfv[ok]))      # same cap as limma::eBayes
  }
  # per-column SE (moderated, or HC3)
  se <- matrix(NA_real_, m, k, dimnames = list(NULL, cn))
  for (j in which(ok)) {
    if (robust_se) se[j, ] <- fits[[j]]$se_hc3
    else se[j, ] <- sqrt(var_post[j] * diag(fits[[j]]$V))
  }
  df_use <- if (robust_se) dfv else df_mod
  list(fits = fits, status = status, n_obs = n_obs, s2 = s2, df = dfv, df_mod = df_use,
       var_post = var_post, d0 = d0, beta = beta, se = se, amean = amean, k = k,
       coef_names = cn, robust_se = robust_se, ok = ok)
}

#' HC3 sandwich covariance from a QR fit
#' @keywords internal
sandwich_hc3 <- function(Xd, res, XtXi) {
  h <- rowSums((Xd %*% XtXi) * Xd)                    # leverages
  u <- res / (1 - h)
  meat <- crossprod(Xd * u)
  XtXi %*% meat %*% XtXi
}

#' Test statistics for the abundance arm after (optional) centring
#'
#' @param arm output of fit_abundance_arm()
#' @param center output of center_effects() or NULL
#' @return data.frame with per-feature omnibus p, per-column estimates/SE/p
#' @keywords internal
abundance_tests <- function(arm, center = NULL) {
  m <- length(arm$status); k <- arm$k; cn <- arm$coef_names
  est <- arm$beta; se <- arm$se
  if (!is.null(center)) {
    est <- sweep(est, 2L, center$delta, "-")
    se  <- sqrt(sweep(se^2, 2L, center$se_delta^2, "+"))
  }
  tstat <- est / se
  p_col <- 2 * stats::pt(-abs(tstat), df = arm$df_mod)
  p_omni <- p_col[, 1L]
  if (k > 1L) {
    p_omni <- rep(NA_real_, m)
    for (j in which(arm$ok)) {
      V <- if (arm$robust_se) arm$fits[[j]]$V_hc3 else arm$var_post[j] * arm$fits[[j]]$V
      if (!is.null(center)) V <- V + diag(center$se_delta^2, k)
      Fj <- tryCatch(as.numeric(t(est[j, ]) %*% solve(V, est[j, ])) / k, error = function(e) NA_real_)
      p_omni[j] <- stats::pf(Fj, k, arm$df_mod[j], lower.tail = FALSE)
    }
  }
  out <- data.frame(abund_n_obs = arm$n_obs, abund_status = arm$status, abund_p = p_omni,
                    stringsAsFactors = FALSE)
  for (i in seq_len(k)) {
    out[[paste0("abund_lfc2_", cn[i])]] <- est[, i] / log(2)
    out[[paste0("abund_se2_",  cn[i])]] <- se[, i] / log(2)
    out[[paste0("abund_t_",    cn[i])]] <- tstat[, i]
    if (k > 1L) out[[paste0("abund_p_", cn[i])]] <- p_col[, i]
  }
  if (k == 1L) {
    names(out) <- sub(paste0("_", cn[1L], "$"), "", names(out))
    ci <- stats::qt(0.975, df = arm$df_mod) * out$abund_se2
    out$abund_ci2_lo <- out$abund_lfc2 - ci; out$abund_ci2_hi <- out$abund_lfc2 + ci
  }
  out
}
