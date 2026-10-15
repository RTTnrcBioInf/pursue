# -----------------------------------------------------------------------------
# 91_rm_refs.R -- subject-aware versions of the reference comparators, for repeated-measures cells (dev/tools/rm.R,
# 2026-10-15). The benchmark runs every comparator without subject structure (R23: they ignore subjects); these give
# each one the subject handling its package supports, so the comparison is against the best available use:
#   ref_linda_re       LinDA with a subject random intercept, formula ~ group + (1|subject) (lmerTest, Satterthwaite)
#   ref_ldm_cl         LDM with cluster.id = subject: whole subjects permuted when the exposure is fixed per subject,
#                      permutation within subject when it varies within subject
#   ref_zicoseq_st     ZicoSeq with strata = subject (permutation within subject); within-subject exposure only
#   ref_logistic_cr    logistic regression on presence with a subject-clustered (CR1) sandwich, t on G - 1 df
# All read meta$subject. A method whose package is missing errors in its cells and is absent from the comparison.
# -----------------------------------------------------------------------------
local({
  root <- Sys.getenv("PURSUE_BENCH_ROOT", unset = "benchmarks")
  e <- new.env()
  for (f in c("R/methods/elementary.R", "R/methods/external.R")) sys.source(file.path(root, f), envir = e)
  within_subject <- function(meta, tested_term) {
    if (is.null(meta$subject)) stop("meta has no subject column")
    any(tapply(as.character(meta[[tested_term]]), meta$subject, function(v) length(unique(v)) > 1L)) }
  fin <- function(r) {
    if (all(!is.finite(r$p))) stop("comparator returned no p-values: ", paste(unique(r$status), collapse = ","))
    r <- r[r$arm %in% c("single", "combined") | !duplicated(r$feature), ]
    data.frame(feature = r$feature, p = r$p, q = r$q, estimate = r$estimate) }

  register_candidate("ref_linda_re", function(counts, meta, formula, tested_term) {
    within_subject(meta, tested_term)
    feats <- rownames(counts); cn <- e$.coef_name(formula, meta, tested_term)
    a <- list(feature.dat.type = "count", prev.filter = 0, mean.abund.filter = 0, max.abund.filter = 0, is.winsor = TRUE,
              outlier.pct = 0.03, adaptive = TRUE, p.adj.method = "BH", alpha = 0.05, n.cores = 1, verbose = FALSE)
    out <- e$.quiet(do.call(MicrobiomeStat::linda, c(list(feature.dat = counts, meta.dat = meta,
                                                            formula = paste(e$.formula_chr(formula), "+ (1|subject)")), a)))
    tab <- out$output[[cn]]; tab <- tab[match(feats, rownames(tab)), ]
    fin(data.frame(feature = feats, arm = "single", p = tab$pvalue, q = tab$padj, estimate = tab$log2FoldChange, status = NA))
  }, notes = "LinDA with (1|subject) random intercept")

  register_candidate("ref_ldm_cl", function(counts, meta, formula, tested_term) {
    w <- within_subject(meta, tested_term)
    fin(e$method_ldm(counts, meta, formula, tested_term, args = list(cluster.id = "subject",
          perm.within.type = if (w) "free" else "none", perm.between.type = if (w) "none" else "free")))
  }, notes = "LDM with cluster.id = subject (between-subject: whole subjects permuted; within-subject: permutation within subject)")

  register_candidate("ref_zicoseq_st", function(counts, meta, formula, tested_term) {
    if (!within_subject(meta, tested_term)) stop("not applicable: exposure is fixed per subject (ZicoSeq permutes within strata only)")
    fin(e$method_zicoseq(counts, meta, formula, tested_term, args = list(strata = factor(meta$subject))))
  }, notes = "ZicoSeq with strata = subject (within-subject exposure only)")

  register_candidate("ref_logistic_cr", function(counts, meta, formula, tested_term) {
    within_subject(meta, tested_term)
    feats <- rownames(counts); D <- counts > 0; cl <- factor(meta$subject); G <- nlevels(cl); Zc <- stats::model.matrix(~ cl - 1)
    X <- stats::model.matrix(formula, meta); a <- attr(X, "assign"); labs <- attr(stats::terms(formula), "term.labels")
    tc <- which(a == which(labs == tested_term)); if (length(tc) != 1L) stop("needs a one-df tested term")
    res <- t(vapply(seq_along(feats), function(j) {
      d <- as.integer(D[j, ]); if (sum(d) < 3 || sum(d) > length(d) - 3) return(c(NA_real_, NA_real_))
      f <- tryCatch(suppressWarnings(stats::glm.fit(X, d, family = stats::binomial())), error = function(err) NULL)
      if (is.null(f) || !f$converged) return(c(NA_real_, NA_real_))
      mu <- f$fitted.values; W <- mu * (1 - mu); Bi <- tryCatch(solve(crossprod(X * sqrt(W))), error = function(err) NULL)
      if (is.null(Bi)) return(c(NA_real_, NA_real_))
      U <- crossprod(Zc, X * (d - mu)); V <- Bi %*% crossprod(U) %*% Bi * G / (G - 1)
      b <- f$coefficients[tc]; se <- sqrt(max(V[tc, tc], 0)); if (!is.finite(b) || se <= 0 || abs(b) > 15) return(c(NA_real_, NA_real_))
      c(2 * stats::pt(-abs(b / se), G - 1), b) }, numeric(2)))
    data.frame(feature = feats, p = res[, 1], q = { q <- rep(NA_real_, nrow(res)); ok <- is.finite(res[, 1]); q[ok] <- stats::p.adjust(res[ok, 1], "BH"); q },
               estimate = res[, 2])
  }, notes = "logistic presence with subject-clustered CR1 sandwich, t on G - 1 df")
})
