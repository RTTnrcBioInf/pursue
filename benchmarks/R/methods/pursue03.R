# -----------------------------------------------------------------------------
# pursue03.R -- PURSUE 0.3 candidates from the R&D programme, as benchmark methods.
# The implementations live in benchmarks/dev/candidates/ (05_enull.R, 60_erd.R) and are sourced from
# there into a private environment, so the benchmark runs exactly the code the dev suite scored.
# Design and results: claude/rd-notebook.md (it5, it8, it11, srv2, srv3).
#
#   pursue03_erdl  expected-rarefaction test, detection AND log transforms, depth-scaling
#                  compositional centring, max(|z|) over the pair (it11 `erdl_max`)
#   pursue03_erdc  expected rarefied detection with depth-scaling centring (it8 `erd_c`)
#   pursue03_erlc  expected rarefied log count with depth-scaling centring (it11 `erl_c`)
#   pursue03_erdlu pairwise common-depth U-statistics, detection AND log, max(|z|) (it12 `erdl_u`)
#   pursue03_erdu  pairwise common-depth detection U-statistic (it12 `erd_u2`)
#   pursue03_erlu  pairwise common-depth log-count U-statistic (it12 `erl_u`)
# The U-statistic methods need a binary exposure with no other terms and no repeated subjects;
# otherwise they fall back to the depth-centred LM versions (erdl / erdc / erlc).
# Estimates: erdl / erlc report the tested coefficient on the expected-rarefied log-count scale,
# divided by log(2) -- close to log2 FC for common taxa, attenuated for rare ones. erdc's
# coefficient is a difference in detection probability, not a fold change, so it reports none.
# All three use only base R + stats; cluster-robust variance when meta$subject repeats.
# -----------------------------------------------------------------------------
.p03 <- local({
  e <- new.env(parent = globalenv())
  e$register_candidate <- function(name, fn, notes = "") invisible(NULL)
  dir <- file.path(Sys.getenv("PURSUE_BENCH_ROOT", unset = "."), "dev", "candidates")
  for (f in c("05_enull.R", "60_erd.R")) sys.source(file.path(dir, f), envir = e)
  e
})
.p03_result <- function(feats, p, est = NA_real_) {
  .finish(data.frame(feature = feats, arm = "single", p = p, q = NA, estimate = est, se = NA, ci_lo = NA, ci_hi = NA,
                     status = ifelse(is.finite(p), NA, "not_tested"), stringsAsFactors = FALSE))
}
method_pursue03_erdl <- function(counts, meta, formula, tested_term, args = list()) {
  v <- .p03$.erdl_fit(counts, meta, formula, tested_term); i <- match(rownames(counts), v$feature)
  .p03_result(rownames(counts), v$p_max[i], v$est_log[i] / log(2))
}
method_pursue03_erdc <- function(counts, meta, formula, tested_term, args = list()) {
  r <- .p03$.erd_c(counts, meta, formula, tested_term); i <- match(rownames(counts), r$feature)
  .p03_result(rownames(counts), r$p[i])
}
method_pursue03_erlc <- function(counts, meta, formula, tested_term, args = list()) {
  v <- .p03$.erdl_fit(counts, meta, formula, tested_term); i <- match(rownames(counts), v$feature)
  .p03_result(rownames(counts), v$p_log[i], v$est_log[i] / log(2))
}

method_pursue03_erdlu <- function(counts, meta, formula, tested_term, args = list()) {
  v <- .p03$.eu_fit(counts, meta, formula, tested_term); i <- match(rownames(counts), v$feature)
  .p03_result(rownames(counts), v$p_max[i], v$est[i])
}
method_pursue03_erdu <- function(counts, meta, formula, tested_term, args = list()) {
  v <- .p03$.eu_fit(counts, meta, formula, tested_term); i <- match(rownames(counts), v$feature)
  .p03_result(rownames(counts), v$p_det[i])
}
method_pursue03_erlu <- function(counts, meta, formula, tested_term, args = list()) {
  v <- .p03$.eu_fit(counts, meta, formula, tested_term); i <- match(rownames(counts), v$feature)
  .p03_result(rownames(counts), v$p_log[i], v$est[i])
}

# it13: the pairwise tests with both pair members thinned to rho x their common depth (srv4 showed
# rho = 1 is anti-conservative on mid under depth confounding). _r07 = rho 0.7, _r05 = rho 0.5.
for (.r in c(0.7, 0.5)) local({ r <- .r; tag <- sprintf("r%02d", round(10 * r))
  assign(paste0("method_pursue03_erdlu_", tag), function(counts, meta, formula, tested_term, args = list()) {
    v <- .p03$.eu_fit(counts, meta, formula, tested_term, rho = r); i <- match(rownames(counts), v$feature)
    .p03_result(rownames(counts), v$p_max[i], v$est[i]) }, envir = globalenv())
  assign(paste0("method_pursue03_erdu_", tag), function(counts, meta, formula, tested_term, args = list()) {
    v <- .p03$.eu_fit(counts, meta, formula, tested_term, rho = r); i <- match(rownames(counts), v$feature)
    .p03_result(rownames(counts), v$p_det[i]) }, envir = globalenv())
})
