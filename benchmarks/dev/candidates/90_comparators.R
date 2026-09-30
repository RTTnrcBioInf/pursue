# -----------------------------------------------------------------------------
# 90_comparators.R -- the benchmark's strongest comparators as dev-suite REFERENCES (not candidates):
# ADAPT, ZicoSeq, LinDA, LDM through the benchmark's own wrappers (benchmarks/R/methods/external.R),
# scored on the dev templates and seeds exactly like the candidates, each with its own q where it
# reports one (as the benchmark scores it). They give the inner loop a per-setting target: where, on
# the tuning templates, a calibrated comparator finds more than the lead. A comparator whose package
# is not installed errors in its cells and is simply absent from the comparison.
# -----------------------------------------------------------------------------
local({
  root <- Sys.getenv("PURSUE_BENCH_ROOT", unset = "benchmarks")
  e <- new.env()
  for (f in c("R/methods/elementary.R", "R/methods/external.R")) sys.source(file.path(root, f), envir = e)
  wrap <- function(fn) { force(fn); function(counts, meta, formula, tested_term) {
    r <- e[[fn]](counts, meta, formula, tested_term)
    if (all(!is.finite(r$p))) stop("comparator returned no p-values: ", paste(unique(r$status), collapse = ","))
    r <- r[r$arm %in% c("single", "combined") | !duplicated(r$feature), ]
    data.frame(feature = r$feature, p = r$p, q = r$q, estimate = r$estimate)
  } }
  for (m in c("adapt", "zicoseq", "linda", "ldm"))
    register_candidate(paste0("ref_", m), wrap(paste0("method_", m)), notes = paste("reference comparator:", m, "(benchmark wrapper, own q)"))
})
