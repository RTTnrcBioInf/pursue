# -----------------------------------------------------------------------------
# 92_refs2.R -- the comparators conf1 / conf2 left out of axes A-C, as dev-suite REFERENCES (2026-11-04): MaAsLin 3,
# LOCOM2, ANCOM-BC2, ALDEx2, corncob through the benchmark's own wrappers (benchmarks/R/methods/external.R), each with
# its own q. MaAsLin 3 reports three arms (abundance, prevalence, combined); as in the benchmark's tables the combined
# arm is the call (90_comparators.R's wrapper keeps the first row per feature, which for MaAsLin 3 would be the
# abundance arm, so these get their own wrapper). Which package version runs is whatever the R library path finds
# first; hpc/run_srv35.sh records it.
# -----------------------------------------------------------------------------
local({
  root <- Sys.getenv("PURSUE_BENCH_ROOT", unset = "benchmarks")
  e <- new.env()
  for (f in c("R/methods/elementary.R", "R/methods/external.R")) sys.source(file.path(root, f), envir = e)
  wrap <- function(fn) { force(fn); function(counts, meta, formula, tested_term) {
    r <- e[[fn]](counts, meta, formula, tested_term)
    k <- r[r$arm %in% c("single", "combined"), , drop = FALSE]; if (!nrow(k)) k <- r[!duplicated(r$feature), , drop = FALSE]
    if (all(!is.finite(k$p))) stop("comparator returned no p-values: ", paste(unique(r$status), collapse = ","))
    data.frame(feature = k$feature, p = k$p, q = k$q, estimate = k$estimate)
  } }
  # LOCOM2 1.1 (CRAN 2026-10-29) replaced locom2(otu.table, Y, C, ...) with a formula interface, locom2(otu.table | C ~ Y),
  # and dropped otu.table / Y / C as arguments (srv35's smoke: "unused arguments (otu.table = ...)"). The benchmark wrapper
  # targets 1.0; this one calls whichever interface the installed version has. Same choices as the benchmark: the Wald
  # test (fall back to permutation, then asymptotic), its own q, no internal filter (the common prevalence filter is applied
  # first), fdr.nominal 0.05, seed 1, one core.
  e$method_locom2_any <- function(counts, meta, formula, tested_term, args = list()) {
    fn <- get("locom2", envir = asNamespace("LOCOM2"))
    if (!"formula" %in% names(formals(fn))) return(e$method_locom2(counts, meta, formula, tested_term, args))
    feats <- rownames(counts)
    if (!e$.tested_is_binary(meta, tested_term)) return(e$.empty_result(feats, status = "not_applicable_needs_binary"))
    g <- factor(meta[[tested_term]]); ev <- new.env(); ev$otu <- t(counts); ev$Y <- as.integer(g == levels(g)[2])
    C <- e$.nuisance_matrix(formula, meta, tested_term); if (!is.null(C)) ev$C <- C
    fo <- stats::as.formula(if (is.null(C)) "otu ~ Y" else "otu | C ~ Y", env = ev)
    out <- e$.quiet(fn(formula = fo, fdr.nominal = 0.05, filter = FALSE, seed = 1, n.cores = 1, verbose = FALSE))
    p <- NULL; pnm <- NA_character_
    for (nm in c("p.otu.Wald", "p.otu.perm", "p.otu.asymptotic")) { p <- e$.as_pvec(out[[nm]], feats); if (!is.null(p)) { pnm <- nm; break } }
    if (is.null(p)) return(e$.empty_result(feats, status = "no_pvalues"))
    q <- e$.as_pvec(out[[sub("^p", "q", pnm)]], feats)
    res <- e$.finish(data.frame(feature = feats, arm = "single", p = p, q = NA, estimate = NA, se = NA, ci_lo = NA, ci_hi = NA,
                                status = paste("LOCOM2", utils::packageVersion("LOCOM2"), pnm), stringsAsFactors = FALSE))
    if (!is.null(q) && !all(is.na(q))) res$q <- q
    res
  }
  for (m in c("maaslin3", "ancombc2", "aldex2", "corncob"))
    register_candidate(paste0("ref_", m), wrap(paste0("method_", m)), notes = paste("reference comparator:", m, "(benchmark wrapper, own q; combined arm for MaAsLin 3)"))
  register_candidate("ref_locom2", wrap("method_locom2_any"), notes = "reference comparator: LOCOM2 (Wald test, own q; 1.0 or 1.1 interface)")
})
