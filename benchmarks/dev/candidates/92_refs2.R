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
  for (m in c("maaslin3", "locom2", "ancombc2", "aldex2", "corncob"))
    register_candidate(paste0("ref_", m), wrap(paste0("method_", m)), notes = paste("reference comparator:", m, "(benchmark wrapper, own q; combined arm for MaAsLin 3)"))
})
