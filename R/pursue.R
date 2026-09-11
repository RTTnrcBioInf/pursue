# -----------------------------------------------------------------------------
# pursue.R -- user entry point
# -----------------------------------------------------------------------------

#' Two-part differential abundance analysis
#'
#' Fits two separately reported components for every feature of a count table:
#'
#' * **Abundance**: a linear model on log relative abundance over detected cells
#'   (X > 0) with a centred log-depth covariate, mild winsorization, empirical-Bayes
#'   variance moderation across features (limma), and empirical-null mixture centring
#'   of the tested effect across features so that the reported log fold change is
#'   relative to the community's typical (null) feature. Optional HC3 sandwich SEs.
#' * **Presence**: a per-feature zero-inflated negative binomial with a depth offset,
#'   cross-feature dispersion shrinkage, a ridge-penalised presence logit, and a
#'   likelihood-ratio test of the tested term on structural presence.
#'
#' Both components use the same formula. The abundance estimand is the log2 fold change in
#' relative abundance among samples where the feature is present, relative to the
#' empirical-null centre; it equals the absolute-abundance fold change only under the
#' assumption that the null features form the modal cluster of effects (diagnostics
#' report `pi0` and the centre). The presence estimand is the change in log odds of
#' structural presence, net of detection.
#'
#' @param otu numeric count matrix or data.frame, **rows = samples, columns = features**.
#' @param meta data.frame of sample metadata, rows aligned with `otu`.
#' @param formula RHS-only model formula over columns of `meta`, e.g. `~ group + age`.
#' @param tested_term the term of `formula` to test (e.g. `"group"`). Continuous terms
#'   and factors with several levels are both supported (multi-df tests).
#' @param min_prevalence features detected in fewer than this fraction of samples are
#'   not tested (default 0.10). Set 0 to test everything with enough detections.
#' @param min_nonzero minimum detected cells for the abundance arm (default 8).
#' @param depth_adjust add centred log-depth as an abundance-arm covariate (default TRUE).
#' @param winsorize,winsor_quantile,winsor_min_n winsorize log abundances at the given
#'   quantile when at least `winsor_min_n` cells are detected (default 0.03 / 20).
#' @param robust_se use HC3 sandwich SEs in the abundance arm instead of moderated
#'   variances (recommended for strongly unbalanced heteroscedastic designs).
#' @param center centre abundance effects across features with the empirical-null
#'   mixture (default TRUE).
#' @param d0_prior prior weight (pseudo-observations) for dispersion shrinkage in the
#'   presence arm (default 10). @param tau_ridge SD of the ridge on presence coefficients.
#' @param combine also report a Cauchy-combined "any effect" p-value (default TRUE).
#' @param q_alpha FDR threshold used for the `*_sig` flags (default 0.05).
#' @param n_cores cores for the per-feature loops (forking; default 1).
#' @param verbose print progress.
#' @return a list of class `pursue`: `results` (one row per feature), `centre`
#'   (per-column empirical-null centre and its SE on the log2 scale, and pi0),
#'   `diagnostics`, `design`, `call`.
#' @examples
#' set.seed(1)
#' n <- 40; m <- 60
#' meta <- data.frame(group = factor(rep(c("a", "b"), each = n / 2)), age = rnorm(n))
#' depth <- round(exp(rnorm(n, log(2e4), 0.3)))
#' lam <- exp(matrix(rnorm(n * m, -6, 1.5), n, m))
#' lam[meta$group == "b", 1:6] <- lam[meta$group == "b", 1:6] * 4
#' otu <- matrix(rpois(n * m, depth * lam), n, m, dimnames = list(NULL, paste0("f", 1:m)))
#' res <- pursue(otu, meta, ~ group + age, "group", verbose = FALSE)
#' head(res$results[order(res$results$abund_p), c("feature", "abund_lfc2", "abund_q", "pres_p")])
#' @export
pursue <- function(otu, meta, formula, tested_term,
                   min_prevalence = 0.10, min_nonzero = 8L, depth_adjust = TRUE,
                   winsorize = TRUE, winsor_quantile = 0.03, winsor_min_n = 20L,
                   robust_se = FALSE, center = TRUE,
                   d0_prior = 10, tau_ridge = 5,
                   combine = TRUE, q_alpha = 0.05, n_cores = 1L, verbose = TRUE) {
  t0 <- Sys.time()
  otu <- as_count_matrix(otu); meta <- align_meta(otu, meta)
  if (nrow(otu) < 4L) stop("Need at least 4 samples.")
  depth <- rowSums(otu)
  if (any(depth <= 0)) stop("Every sample must have a positive library size.")
  design <- build_design(formula, meta, tested_term, depth)

  prev <- colMeans(otu > 0)
  tested <- prev >= min_prevalence & colSums(otu > 0) >= 3
  if (!any(tested)) stop("No feature passes the prevalence filter.")
  if (verbose) message(sprintf("PURSUE: %d samples, %d features (%d tested), term '%s' (%d df)",
                               nrow(otu), ncol(otu), sum(tested), tested_term, design$n_tested))
  otu_t <- otu[, tested, drop = FALSE]

  if (verbose) message("  abundance arm ...")
  arm <- fit_abundance_arm(otu_t, depth, design, depth_adjust = depth_adjust, min_nonzero = min_nonzero,
                           winsorize = winsorize, winsor_quantile = winsor_quantile, winsor_min_n = winsor_min_n,
                           robust_se = robust_se, n_cores = n_cores)
  ctr <- if (center) center_effects(arm) else NULL
  ab <- abundance_tests(arm, ctr)

  if (verbose) message("  presence arm ...")
  pr <- fit_presence_arm(otu_t, depth, design, d0_prior = d0_prior, tau_ridge = tau_ridge, n_cores = n_cores)

  res <- data.frame(feature = colnames(otu_t), prevalence = prev[tested], stringsAsFactors = FALSE)
  res <- cbind(res, ab, pr)
  res$abund_q <- stats::p.adjust(res$abund_p, "BH")
  res$pres_q  <- stats::p.adjust(res$pres_p, "BH")
  if (combine) {
    res$any_p <- cauchy_combine(cbind(res$abund_p, res$pres_p))
    res$any_q <- stats::p.adjust(res$any_p, "BH")
  }
  res$abund_sig <- !is.na(res$abund_q) & res$abund_q <= q_alpha
  res$pres_sig  <- !is.na(res$pres_q)  & res$pres_q  <= q_alpha
  # features not tested: placeholder rows
  if (any(!tested)) {
    skip <- res[rep(NA_integer_, sum(!tested)), , drop = FALSE]
    skip$feature <- colnames(otu)[!tested]; skip$prevalence <- prev[!tested]
    skip$abund_status <- "below_prevalence_filter"; skip$pres_status <- "below_prevalence_filter"
    skip$abund_sig <- FALSE; skip$pres_sig <- FALSE
    res <- rbind(res, skip)
  }
  res <- res[match(colnames(otu), res$feature), , drop = FALSE]; rownames(res) <- NULL

  diag <- list(
    n_samples = nrow(otu), n_features = ncol(otu), n_tested = sum(tested),
    depth_summary = stats::quantile(depth, c(0, .25, .5, .75, 1)),
    abundance = list(d0_prior_df = arm$d0, n_ok = sum(arm$ok), status_table = table(arm$status),
                     robust_se = robust_se),
    # reported on the log2 scale, like the abundance effects in `results`
    centre = if (center) list(delta = ctr$delta / log(2), se_delta = ctr$se_delta / log(2), pi0 = ctr$pi0,
                              median_delta = ctr$median_delta / log(2), notes = ctr$notes) else NULL,
    presence = list(status_table = table(pr$pres_status), theta_range = range(pr$pres_theta, na.rm = TRUE)),
    both_arms_significant = sum(res$abund_sig & res$pres_sig, na.rm = TRUE),
    heteroscedasticity_hint = hetero_hint(otu_t, meta, design),
    runtime_s = as.numeric(Sys.time() - t0, units = "secs"))
  if (verbose) {
    message(sprintf("  done in %.1f s: %d abundance calls, %d presence calls at q <= %g%s",
                    diag$runtime_s, sum(res$abund_sig), sum(res$pres_sig), q_alpha,
                    if (center) sprintf("; centre = %s, pi0 = %s",
                                        paste(round(ctr$delta / log(2), 2), collapse = "/"),
                                        paste(round(ctr$pi0, 2), collapse = "/")) else ""))
  }
  structure(list(results = res, centre = diag$centre, diagnostics = diag, design = design,
                 call = match.call()), class = "pursue")
}

#' @export
print.pursue <- function(x, ...) {
  d <- x$diagnostics
  cat("PURSUE two-part differential abundance\n")
  cat(sprintf("  %d samples, %d features (%d tested); tested term: %s\n", d$n_samples, d$n_features, d$n_tested, x$design$tested_term))
  cat(sprintf("  abundance arm: %d features fitted, prior df %.1f%s\n", d$abundance$n_ok,
              ifelse(is.na(d$abundance$d0_prior_df), NA, d$abundance$d0_prior_df),
              if (d$abundance$robust_se) " (HC3 SEs)" else ""))
  if (!is.null(x$centre)) cat(sprintf("  centre (log2): %s; pi0: %s\n",
                                     paste(round(x$centre$delta, 3), collapse = ", "),
                                     paste(round(x$centre$pi0, 2), collapse = ", ")))
  cat(sprintf("  calls at q<=0.05: %d abundance, %d presence, %d both\n",
              sum(x$results$abund_sig), sum(x$results$pres_sig), d$both_arms_significant))
  invisible(x)
}

#' Cauchy (ACAT) combination of p-values across the columns of a matrix
#' @keywords internal
cauchy_combine <- function(P, weights = NULL) {
  P <- as.matrix(P); k <- ncol(P)
  if (is.null(weights)) weights <- rep(1 / k, k)
  apply(P, 1L, function(p) {
    ok <- is.finite(p); if (!any(ok)) return(NA_real_)
    p <- pmin(pmax(p[ok], 1e-15), 1 - 1e-15); w <- weights[ok] / sum(weights[ok])
    T <- sum(w * tan((0.5 - p) * pi))
    if (T > 1e15) return(1 / (T * pi))
    0.5 - atan(T) / pi
  })
}

#' Rough heteroscedasticity screen for a binary tested term
#' @keywords internal
hetero_hint <- function(otu, meta, design) {
  if (design$n_tested != 1L) return(NULL)
  xcol <- design$X[, design$tested_cols[1L]]
  if (length(unique(xcol)) != 2L) return(NULL)
  g <- xcol == max(xcol); n1 <- sum(g); n0 <- sum(!g)
  depth <- rowSums(otu)
  Y <- log((otu + 0.5) / depth)
  v1 <- apply(Y[g, , drop = FALSE], 2L, stats::var); v0 <- apply(Y[!g, , drop = FALSE], 2L, stats::var)
  ratio <- stats::median(v1 / v0, na.rm = TRUE)
  list(group_sizes = c(n0, n1), imbalance = max(n0, n1) / min(n0, n1), median_variance_ratio = ratio,
       suggest_robust_se = (max(n0, n1) / min(n0, n1) > 2) && (ratio > 2 || ratio < 0.5))
}
