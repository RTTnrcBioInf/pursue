# -----------------------------------------------------------------------------
# external.R -- wrappers for published methods. Same signature as elementary.R.
#
# These are adapted from the wrappers that ran in the original PURSUE benchmark
# (benchmark scripts/*.R, verified on the cluster) and extended to arbitrary formulas.
# The ones marked [UNVERIFIED] were written against the package documentation and have
# not been executed here (CRAN/Bioconductor are not reachable from the development
# container). hpc/smoke_test_methods.R runs every wrapper on a small dataset; run it
# on the cluster before any benchmark job.
# -----------------------------------------------------------------------------

.coef_name <- function(formula, meta, tested_term) {
  X <- stats::model.matrix(formula, meta); a <- attr(X, "assign"); labs <- attr(stats::terms(formula), "term.labels")
  colnames(X)[a == which(labs == tested_term)]
}
.formula_chr <- function(formula) paste(deparse(formula), collapse = "")
.rhs_chr <- function(formula) sub("^~\\s*", "", .formula_chr(formula))
.nuisance_matrix <- function(formula, meta, tested_term) {
  X <- stats::model.matrix(formula, meta); a <- attr(X, "assign"); labs <- attr(stats::terms(formula), "term.labels")
  Z <- X[, a != which(labs == tested_term) & a != 0, drop = FALSE]; if (ncol(Z) == 0L) NULL else Z
}
.quiet <- function(expr) { r <- NULL; suppressMessages(suppressWarnings(invisible(utils::capture.output(r <- expr, file = nullfile())))); r }

# ---- LinDA (MicrobiomeStat) -- verified pattern ----
method_linda <- function(counts, meta, formula, tested_term, args = list()) {
  feats <- rownames(counts); cn <- .coef_name(formula, meta, tested_term)
  if (length(cn) != 1L) return(.empty_result(feats, status = "not_applicable_multi_df"))
  a <- utils::modifyList(list(feature.dat.type = "count", prev.filter = 0, mean.abund.filter = 0, max.abund.filter = 0,
                              is.winsor = TRUE, outlier.pct = 0.03, adaptive = TRUE, p.adj.method = "BH", alpha = 0.05, n.cores = 1, verbose = FALSE), args)
  out <- .quiet(do.call(MicrobiomeStat::linda, c(list(feature.dat = counts, meta.dat = meta, formula = .formula_chr(formula)), a)))
  tab <- out$output[[cn]]; tab <- tab[match(feats, rownames(tab)), ]
  .finish(data.frame(feature = feats, arm = "single", p = tab$pvalue, q = NA, estimate = tab$log2FoldChange, se = tab$lfcSE,
                     ci_lo = tab$log2FoldChange - 1.96 * tab$lfcSE, ci_hi = tab$log2FoldChange + 1.96 * tab$lfcSE, status = NA, stringsAsFactors = FALSE))
}

# ---- ANCOM-BC2 ----
# Calling convention drifted: ANCOMBC <= 2.2 took (data = matrix, taxa_are_rows, meta_data);
# 2.4 dropped both and wants a TreeSummarizedExperiment / phyloseq in `data`. The cluster smoke
# test of 2026-09-11 (ANCOMBC 2.4.0) failed on the old form, so the call is built from the
# installed version's formals.
.ancombc2_data_args <- function(counts, meta) {
  fm <- names(formals(ANCOMBC::ancombc2))
  if ("taxa_are_rows" %in% fm) return(list(data = counts, taxa_are_rows = TRUE, meta_data = meta))
  if (requireNamespace("TreeSummarizedExperiment", quietly = TRUE))
    return(list(data = TreeSummarizedExperiment::TreeSummarizedExperiment(
      assays = list(counts = as.matrix(counts)), colData = meta), assay_name = "counts"))
  if (requireNamespace("phyloseq", quietly = TRUE))
    return(list(data = phyloseq::phyloseq(phyloseq::otu_table(as.matrix(counts), taxa_are_rows = TRUE),
                                          phyloseq::sample_data(meta))))
  stop("ANCOMBC >= 2.4 needs TreeSummarizedExperiment or phyloseq to carry the count table")
}

method_ancombc2 <- function(counts, meta, formula, tested_term, args = list()) {
  feats <- rownames(counts); cn <- .coef_name(formula, meta, tested_term)
  if (length(cn) != 1L) return(.empty_result(feats, status = "not_applicable_multi_df"))
  a <- utils::modifyList(list(p_adj_method = "BH", pseudo_sens = TRUE, prv_cut = 0, lib_cut = 0, s0_perc = 0.05,
                              struc_zero = FALSE, neg_lb = FALSE, alpha = 0.05, n_cl = 1, verbose = FALSE,
                              global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE), args)
  out <- .quiet(do.call(ANCOMBC::ancombc2,
                        c(.ancombc2_data_args(counts, meta), list(fix_formula = .rhs_chr(formula)), a)))
  r <- as.data.frame(out$res); r <- r[match(feats, r$taxon), , drop = FALSE]
  get <- function(pfx) if (paste0(pfx, cn) %in% names(r)) r[[paste0(pfx, cn)]] else NULL
  p <- as.numeric(get("p_")); qa <- as.numeric(get("q_"))
  lfc <- as.numeric(get("lfc_")) / log(2); se <- as.numeric(get("se_")) / log(2)
  # the authors' recommended call is "differential AND passed the pseudo-count sensitivity check":
  # 2.2 exposed that jointly as diff_robust_*, 2.4 exposes diff_* and passed_ss_* separately.
  rob <- get("diff_robust_")
  if (is.null(rob)) { d <- get("diff_"); ss <- get("passed_ss_")
    rob <- if (is.null(d)) rep(TRUE, nrow(r)) else if (is.null(ss)) as.logical(d) else as.logical(d) & as.logical(ss) }
  rob <- as.logical(rob); rob[is.na(rob)] <- FALSE
  res <- .finish(data.frame(feature = feats, arm = "single", p = p, q = NA, estimate = lfc, se = se,
                            ci_lo = lfc - 1.96 * se, ci_hi = lfc + 1.96 * se,
                            status = ifelse(is.na(p), "not_returned", NA), stringsAsFactors = FALSE))
  # Deviation from the original benchmark wrapper, declared under protocol section 6: p is left raw so
  # that calibration (KS) and pAUC measure the model, and the sensitivity filter gates only the call (q).
  res$q <- ifelse(rob, ifelse(is.na(qa), res$q, qa), 1)
  res
}

# ---- corncob -- verified pattern ----
method_corncob <- function(counts, meta, formula, tested_term, args = list()) {
  feats <- rownames(counts)
  labs <- attr(stats::terms(formula), "term.labels"); nuis <- setdiff(labs, tested_term)
  f_null <- stats::as.formula(paste("~", if (length(nuis)) paste(nuis, collapse = " + ") else "1"))
  a <- utils::modifyList(list(test = "LRT", boot = FALSE, fdr_cutoff = 0.05), args)
  out <- .quiet(do.call(corncob::differentialTest, c(list(formula = formula, phi.formula = formula, formula_null = f_null, phi.formula_null = formula,
                                                          data = counts, sample_data = meta, taxa_are_rows = TRUE), a)))
  p <- out$p[feats]
  .finish(data.frame(feature = feats, arm = "single", p = as.numeric(p), q = NA, estimate = NA, se = NA, ci_lo = NA, ci_hi = NA,
                     status = ifelse(feats %in% c(out$discriminant_taxa_DA, out$discriminant_taxa_DV), "discriminant_excluded", NA), stringsAsFactors = FALSE))
}

# ---- LDM -- verified pattern (formula: y | confounders ~ tested) ----
method_ldm <- function(counts, meta, formula, tested_term, args = list()) {
  feats <- rownames(counts); y <- t(counts); storage.mode(y) <- "numeric"
  labs <- attr(stats::terms(formula), "term.labels"); nuis <- setdiff(labs, tested_term)
  fstr <- if (length(nuis)) paste0("y | ", paste(nuis, collapse = " + "), " ~ ", tested_term) else paste0("y ~ ", tested_term)
  # LDM resolves the response by name from the GLOBAL environment -- not from `data`, and not
  # from the formula's environment. The 2026-09-11 smoke test failed with "object 'y' not
  # found" when it was only in a local env; this is the pattern the original benchmark used.
  had <- exists("y", envir = globalenv(), inherits = FALSE)
  old_y <- if (had) get("y", envir = globalenv()) else NULL
  assign("y", y, envir = globalenv())
  on.exit({ if (had) assign("y", old_y, envir = globalenv()) else suppressWarnings(rm("y", envir = globalenv())) }, add = TRUE)
  a <- utils::modifyList(list(fdr.nominal = 0.05, seed = 1, n.perm.max = 5000, verbose = FALSE), args)
  out <- .quiet(do.call(LDM::ldm, c(list(formula = stats::as.formula(fstr), data = meta), a)))
  pick <- function(m) { m <- as.matrix(m); rn <- rownames(m)
    k <- if (!is.null(rn) && tested_term %in% rn) tested_term else if (!is.null(rn) && "cov1" %in% rn) "cov1" else nrow(m)
    v <- as.numeric(m[k, ]); names(v) <- colnames(m); v }
  p <- pick(out$p.otu.omni)
  q <- tryCatch(pick(out$q.otu.omni), error = function(e) stats::setNames(rep(NA_real_, length(p)), names(p)))
  res <- .finish(data.frame(feature = feats, arm = "single", p = unname(p[feats]), q = NA, estimate = NA, se = NA,
                            ci_lo = NA, ci_hi = NA, status = NA, stringsAsFactors = FALSE))
  # LDM does its own FDR control on permutation p-values; prefer it over BH where available.
  if (!all(is.na(q))) res$q <- unname(q[feats])
  res
}

method_locom <- function(counts, meta, formula, tested_term, args = list()) {
  feats <- rownames(counts)
  if (!.tested_is_binary(meta, tested_term)) return(.empty_result(feats, status = "not_applicable_needs_binary"))
  g <- factor(meta[[tested_term]]); Y <- as.integer(g == levels(g)[2]); C <- .nuisance_matrix(formula, meta, tested_term)
  a <- utils::modifyList(list(fdr.nominal = 0.05, seed = 1, n.perm.max = 20000, n.cores = 1), args)
  call <- c(list(otu.table = t(counts), Y = Y), if (!is.null(C)) list(C = C), a)
  out <- .quiet(do.call(LOCOM::locom, call))
  p <- as.numeric(out$p.otu[1, ]); names(p) <- colnames(out$p.otu); es <- as.numeric(out$effect.size); names(es) <- colnames(out$p.otu)
  .finish(data.frame(feature = feats, arm = "single", p = p[feats], q = NA, estimate = es[feats], se = NA, ci_lo = NA, ci_hi = NA, status = NA, stringsAsFactors = FALSE))
}

# ---- LOCOM2 [UNVERIFIED] -- CRAN package LOCOM2; assumed API mirrors LOCOM (otu.table, Y, C) ----
method_locom2 <- function(counts, meta, formula, tested_term, args = list()) {
  feats <- rownames(counts)
  if (!.tested_is_binary(meta, tested_term)) return(.empty_result(feats, status = "not_applicable_needs_binary"))
  g <- factor(meta[[tested_term]]); Y <- as.integer(g == levels(g)[2]); C <- .nuisance_matrix(formula, meta, tested_term)
  a <- utils::modifyList(list(fdr.nominal = 0.05, seed = 1, n.cores = 1), args)
  fn <- get("locom2", envir = asNamespace("LOCOM2"))
  out <- .quiet(do.call(fn, c(list(otu.table = t(counts), Y = Y), if (!is.null(C)) list(C = C), a)))
  p <- as.numeric(out$p.otu[1, ]); names(p) <- colnames(out$p.otu)
  .finish(data.frame(feature = feats, arm = "single", p = p[feats], q = NA, estimate = NA, se = NA, ci_lo = NA, ci_hi = NA, status = NA, stringsAsFactors = FALSE))
}

# ---- ALDEx2 -- verified pattern (Wilcoxon for binary w/o covariates, GLM otherwise) ----
method_aldex2 <- function(counts, meta, formula, tested_term, args = list()) {
  feats <- rownames(counts); mc <- if (is.null(args$mc.samples)) 128 else args$mc.samples
  if (.tested_is_binary(meta, tested_term) && length(all.vars(formula)) == 1L) {
    g <- factor(meta[[tested_term]]); conds <- as.character(g)
    clr <- .quiet(ALDEx2::aldex.clr(counts, conds, mc.samples = mc, denom = "all", verbose = FALSE))
    tt <- .quiet(ALDEx2::aldex.ttest(clr, paired.test = FALSE, hist.plot = FALSE, verbose = FALSE))
    ef <- tryCatch(.quiet(ALDEx2::aldex.effect(clr, verbose = FALSE)), error = function(e) NULL)
    p <- tt$wi.ep[match(feats, rownames(tt))]; est <- if (!is.null(ef)) ef$diff.btw[match(feats, rownames(ef))] else NA
    return(.finish(data.frame(feature = feats, arm = "single", p = p, q = NA, estimate = est, se = NA, ci_lo = NA, ci_hi = NA, status = NA, stringsAsFactors = FALSE)))
  }
  mm <- stats::model.matrix(formula, meta); cn <- .coef_name(formula, meta, tested_term)
  clr <- .quiet(ALDEx2::aldex.clr(counts, mm, mc.samples = mc, denom = "all", verbose = FALSE))
  gl <- .quiet(ALDEx2::aldex.glm(clr, verbose = FALSE, fdr.method = "BH"))
  pc <- paste0(cn[1], ":pval"); p <- gl[[pc]][match(feats, rownames(gl))]
  .finish(data.frame(feature = feats, arm = "single", p = p, q = NA, estimate = NA, se = NA, ci_lo = NA, ci_hi = NA,
                     status = if (length(cn) > 1L) "multi_df_first_column_only" else NA, stringsAsFactors = FALSE))
}

# ---- MaAsLin 3 -- verified pattern; reported per arm (abundance / prevalence) ----
method_maaslin3 <- function(counts, meta, formula, tested_term, args = list()) {
  feats <- rownames(counts); cn <- .coef_name(formula, meta, tested_term)
  if (length(cn) != 1L) return(.empty_result(feats, status = "not_applicable_multi_df"))
  md <- meta; md$log_depth <- log(pmax(colSums(counts), 1))
  fstr <- paste(.formula_chr(formula), "+ log_depth")
  tmp <- file.path(tempdir(), paste0("maaslin3_", Sys.getpid(), "_", as.integer(stats::runif(1, 1, 1e9)))); dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp, recursive = TRUE, force = TRUE), add = TRUE)
  a <- utils::modifyList(list(normalization = "TSS", transform = "LOG", augment = TRUE, standardize = TRUE, median_comparison_abundance = TRUE,
                              median_comparison_prevalence = FALSE, plot_summary_plot = FALSE, plot_associations = FALSE, save_models = FALSE,
                              max_significance = 1, min_abundance = 0, min_prevalence = 0, max_prevalence = 1.01, min_variance = 0, verbosity = "ERROR"), args)
  .quiet(do.call(maaslin3::maaslin3, c(list(input_data = as.data.frame(t(counts)), input_metadata = md, output = tmp, formula = fstr), a)))
  res <- utils::read.delim(file.path(tmp, "all_results.tsv"), check.names = FALSE)
  res <- res[res$name == cn, , drop = FALSE]
  pick <- function(model) { r <- res[res$model == model, ]; r <- r[match(feats, r$feature), ]
    data.frame(feature = feats, arm = if (model == "abundance") "abundance" else "presence", p = as.numeric(r$pval_individual), q = NA,
               estimate = as.numeric(r$coef), se = as.numeric(r$stderr), ci_lo = NA, ci_hi = NA, status = NA, stringsAsFactors = FALSE) }
  ab <- .finish(pick("abundance")); pr <- .finish(pick("prevalence"))
  rj <- res[match(feats, res$feature), ]; cb <- data.frame(feature = feats, arm = "combined", p = as.numeric(rj$pval_joint), q = NA, estimate = NA, se = NA, ci_lo = NA, ci_hi = NA, status = NA, stringsAsFactors = FALSE)
  rbind(ab, pr, .finish(cb))
}

# ---- ZicoSeq (GUniFrac) -- verified pattern ----
method_zicoseq <- function(counts, meta, formula, tested_term, args = list()) {
  feats <- rownames(counts); labs <- attr(stats::terms(formula), "term.labels"); nuis <- setdiff(labs, tested_term)
  a <- utils::modifyList(list(feature.dat.type = "count", prev.filter = 0, mean.abund.filter = 0, max.abund.filter = 0, min.prop = 0,
                              is.winsor = TRUE, outlier.pct = 0.03, is.post.sample = TRUE, post.sample.no = 25, link.func = list(function(x) sign(x) * sqrt(abs(x))),
                              stats.combine.func = max, perm.no = 999, strata = NULL, ref.pct = 0.5, stage.no = 6, excl.pct = 0.2, is.fwer = FALSE, verbose = FALSE, return.feature.dat = FALSE), args)
  out <- .quiet(do.call(GUniFrac::ZicoSeq, c(list(meta.dat = meta, feature.dat = counts, grp.name = tested_term, adj.name = if (length(nuis)) nuis else NULL), a)))
  p <- out$p.raw[feats]
  .finish(data.frame(feature = feats, arm = "single", p = as.numeric(p), q = NA, estimate = NA, se = NA, ci_lo = NA, ci_hi = NA, status = NA, stringsAsFactors = FALSE))
}

# ---- fastANCOM [UNVERIFIED] -- fastANCOM(Y = samples x features, x = exposure, Z = covariates) ----
method_fastancom <- function(counts, meta, formula, tested_term, args = list()) {
  feats <- rownames(counts)
  x <- meta[[tested_term]]; if (is.factor(x) || is.character(x)) { g <- factor(x); if (nlevels(g) != 2L) return(.empty_result(feats, status = "not_applicable_multi_level")); x <- as.integer(g == levels(g)[2]) }
  Z <- .nuisance_matrix(formula, meta, tested_term)
  out <- .quiet(do.call(fastANCOM::fastANCOM, c(list(Y = t(counts), x = x), if (!is.null(Z)) list(Z = Z), args)))
  r <- out$results$final; r <- r[match(feats, rownames(r)), ]
  .finish(data.frame(feature = feats, arm = "single", p = as.numeric(r$log2FC.pval), q = NA, estimate = as.numeric(r$log2FC), se = as.numeric(r$log2FC.SD),
                     ci_lo = NA, ci_hi = NA, status = NA, stringsAsFactors = FALSE))
}

# ---- ADAPT [UNVERIFIED] -- Bioconductor ADAPT::adapt(input_data = phyloseq, cond.var, base.cond, adj.var) ----
method_adapt <- function(counts, meta, formula, tested_term, args = list()) {
  feats <- rownames(counts)
  if (!.tested_is_binary(meta, tested_term)) return(.empty_result(feats, status = "not_applicable_needs_binary"))
  labs <- attr(stats::terms(formula), "term.labels"); nuis <- setdiff(labs, tested_term)
  ps <- phyloseq::phyloseq(phyloseq::otu_table(counts, taxa_are_rows = TRUE), phyloseq::sample_data(meta))
  g <- factor(meta[[tested_term]])
  out <- .quiet(do.call(ADAPT::adapt, c(list(input_data = ps, cond.var = tested_term, base.cond = levels(g)[1], adj.var = if (length(nuis)) nuis else NULL), args)))
  r <- as.data.frame(out@details); r <- r[match(feats, r$Taxa), ]
  .finish(data.frame(feature = feats, arm = "single", p = as.numeric(r$pval), q = NA, estimate = as.numeric(r$log10foldchange) * log2(10), se = NA, ci_lo = NA, ci_hi = NA, status = NA, stringsAsFactors = FALSE))
}

# ---- radEmu / fastEmu [UNVERIFIED] -- emuFit(formula, Y = samples x features, covariate_data) with robust score tests ----
method_fastemu <- function(counts, meta, formula, tested_term, args = list()) {
  feats <- rownames(counts); cn <- .coef_name(formula, meta, tested_term)
  if (length(cn) != 1L) return(.empty_result(feats, status = "not_applicable_multi_df"))
  fn <- if (requireNamespace("fastEmu", quietly = TRUE)) get("fastEmuFit", envir = asNamespace("fastEmu")) else radEmu::emuFit
  # radEmu renamed `covariate_data` to `data`; 2.3.2 rejects the old name with
  # "both formula and data containing covariates ... must be provided" (smoke test 2026-09-11).
  dat_arg <- if ("data" %in% names(formals(fn))) "data" else "covariate_data"
  call_args <- c(list(formula = formula, Y = t(counts), run_score_tests = TRUE), stats::setNames(list(meta), dat_arg), args)
  out <- .quiet(do.call(fn, call_args))
  co <- out$coef; co <- co[co$covariate == cn, ]; co <- co[match(feats, co$category), ]
  .finish(data.frame(feature = feats, arm = "single", p = as.numeric(co$pval), q = NA, estimate = as.numeric(co$estimate) / log(2), se = as.numeric(co$se) / log(2),
                     ci_lo = as.numeric(co$lower) / log(2), ci_hi = as.numeric(co$upper) / log(2), status = NA, stringsAsFactors = FALSE))
}
