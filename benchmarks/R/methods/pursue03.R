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

# it14: thinning factor chosen from the design (rho 1 when depth is balanced across groups -> 0.5)
method_pursue03_erdlu_ad <- function(counts, meta, formula, tested_term, args = list()) {
  v <- .p03$.eu_fit_ad(counts, meta, formula, tested_term); i <- match(rownames(counts), v$feature)
  .p03_result(rownames(counts), v$p_max[i], v$est[i])
}
method_pursue03_erdu_ad <- function(counts, meta, formula, tested_term, args = list()) {
  v <- .p03$.eu_fit_ad(counts, meta, formula, tested_term); i <- match(rownames(counts), v$feature)
  .p03_result(rownames(counts), v$p_det[i])
}

# it15: covariates inside the pairwise test (pairwise-difference regression); design-chosen rho
method_pursue03_erdlu_uc <- function(counts, meta, formula, tested_term, args = list()) {
  v <- .p03$.euc_fit(counts, meta, formula, tested_term); i <- match(rownames(counts), v$feature)
  .p03_result(rownames(counts), v$p_max[i], v$est[i])
}
method_pursue03_erdu_uc <- function(counts, meta, formula, tested_term, args = list()) {
  v <- .p03$.euc_fit(counts, meta, formula, tested_term); i <- match(rownames(counts), v$feature)
  .p03_result(rownames(counts), v$p_det[i])
}

# it22: censored (Tobit-score) kernel -- zeros imputed at E[latent log count | below detection] at the
# pair's common depth; cendu = max(detection, censored)
method_pursue03_cendu_uc <- function(counts, meta, formula, tested_term, args = list()) {
  v <- .p03$.ecen_fit(counts, meta, formula, tested_term); i <- match(rownames(counts), v$feature)
  .p03_result(rownames(counts), v$p_dc[i], v$est[i])
}
method_pursue03_cenu_uc <- function(counts, meta, formula, tested_term, args = list()) {
  v <- .p03$.ecen_fit(counts, meta, formula, tested_term); i <- match(rownames(counts), v$feature)
  .p03_result(rownames(counts), v$p_cen[i], v$est[i])
}

# it24: erdl_uc + empirical-Bayes variance moderation across taxa; clustered designs through the pairwise
# test with a pooled design effect
method_pursue03_erdlum <- function(counts, meta, formula, tested_term, args = list()) {
  v <- .p03$.eum_fit(counts, meta, formula, tested_term); i <- match(rownames(counts), v$feature)
  .p03_result(rownames(counts), v$p_max[i], v$est[i])
}
# it33: erdl_um or the square-root kernel per taxon, chosen by leave-one-out hit counts over the other taxa
method_pursue03_eselum <- function(counts, meta, formula, tested_term, args = list()) {
  v <- .p03$.euselm_fit(counts, meta, formula, tested_term); i <- match(rownames(counts), v$feature)
  .p03_result(rownames(counts), v$p_sel[i], v$est[i])
}
# it35: size factors + moderation + pairwise clusters + scale chosen across taxa; _b: size factors only
# when depth is balanced across groups
method_pursue03_efull <- function(counts, meta, formula, tested_term, args = list()) {
  v <- .p03$.efull_fit(counts, meta, formula, tested_term); i <- match(rownames(counts), v$feature)
  .p03_result(rownames(counts), v$p_full[i], v$est[i])
}
method_pursue03_efullb <- function(counts, meta, formula, tested_term, args = list()) {
  v <- .p03$.efull_fit(counts, meta, formula, tested_term, sf_balanced_only = TRUE); i <- match(rownames(counts), v$feature)
  .p03_result(rownames(counts), v$p_full[i], v$est[i])
}
# it36 + it37: efull_b on winsorised counts (top 3% per taxon, balanced depth only), three-option selection
method_pursue03_efullb3w <- function(counts, meta, formula, tested_term, args = list()) {
  v <- .p03$.efull_fit(counts, meta, formula, tested_term, sf_balanced_only = TRUE, opts3 = TRUE, winsor = 0.03); i <- match(rownames(counts), v$feature)
  .p03_result(rownames(counts), v$p_full[i], v$est[i])
}
# it38: efull_b3w with exact (Monte-Carlo, studentised) permutation p-values where samples are exchangeable
# under the null (two groups, no covariates, no repeated subjects, balanced depth); efull_b3w elsewhere
method_pursue03_eperm <- function(counts, meta, formula, tested_term, args = list()) {
  v <- .p03$.eperm_fit(counts, meta, formula, tested_term); i <- match(rownames(counts), v$feature)
  .p03_result(rownames(counts), v$p_full[i], v$est[i])
}
# it38r: eperm with the permutation path also at chance-level depth imbalance (rho > 0.5)
method_pursue03_epermr <- function(counts, meta, formula, tested_term, args = list()) {
  v <- .p03$.eperm_fit(counts, meta, formula, tested_term, rho_min = 0.51); i <- match(rownames(counts), v$feature)
  .p03_result(rownames(counts), v$p_full[i], v$est[i])
}
# it39: eperm_r with the compositional centre re-estimated in every permutation (label-free scores + linearised centre)
method_pursue03_epermc <- function(counts, meta, formula, tested_term, args = list()) {
  v <- .p03$.eperm_fit(counts, meta, formula, tested_term, rho_min = 0.51, recentre = TRUE); i <- match(rownames(counts), v$feature)
  .p03_result(rownames(counts), v$p_full[i], v$est[i])
}
# it39 + it40: eperm_c with the scale chosen by BH discoveries among the other taxa
method_pursue03_epermcs <- function(counts, meta, formula, tested_term, args = list()) {
  v <- .p03$.eperm_fit(counts, meta, formula, tested_term, rho_min = 0.51, recentre = TRUE, sel = "bh"); i <- match(rownames(counts), v$feature)
  .p03_result(rownames(counts), v$p_full[i], v$est[i])
}
# it41: eperm_cs with margin 1 (cs1); and with the unthinned lin / sqp scales as further options (x1: margin 1, x2: margin 2)
method_pursue03_epermcs1 <- function(counts, meta, formula, tested_term, args = list()) {
  v <- .p03$.eperm_fit(counts, meta, formula, tested_term, rho_min = 0.51, recentre = TRUE, sel = "bh", margin = 1L); i <- match(rownames(counts), v$feature)
  .p03_result(rownames(counts), v$p_full[i], v$est[i])
}
method_pursue03_epermx1 <- function(counts, meta, formula, tested_term, args = list()) {
  v <- .p03$.eperm_fit(counts, meta, formula, tested_term, rho_min = 0.51, recentre = TRUE, sel = "bh", margin = 1L, extra = c("lin", "sqp")); i <- match(rownames(counts), v$feature)
  .p03_result(rownames(counts), v$p_full[i], v$est[i])
}
method_pursue03_epermx2 <- function(counts, meta, formula, tested_term, args = list()) {
  v <- .p03$.eperm_fit(counts, meta, formula, tested_term, rho_min = 0.51, recentre = TRUE, sel = "bh", margin = 2L, extra = c("lin", "sqp")); i <- match(rownames(counts), v$feature)
  .p03_result(rownames(counts), v$p_full[i], v$est[i])
}
