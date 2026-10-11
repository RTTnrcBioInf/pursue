# --- it49: designed depth confounding -- log alone, thinned harder (Phase 2, target 4: mid R17 / R19) -------------------
# Under designed depth confounding (two groups, rho < 0.51) eperm_c5 runs it46's efull_nodet path: asymptotic pairwise
# U-statistics, log or sqrt chosen per taxon by leave-one-out BH counts. Its weak cells are all there: dev mid R17 FDR
# 0.124, mid R19.null P(any) 0.16. bank2 (mid R17 / R17.null / R19.null, 10 reps per tuning template; local):
# - the abundance-given-presence kernel is NOT the fix it46 suggested: cleaner tail on R19.null hmp_tongue (p < 0.01
#   at 0.79x nominal vs log 1.88x) but not on twinsuk (1.50x vs 1.30x), its null bulk is shifted (median p 0.38), and
#   it has 40-70% of log's power on R17 (TP 28 vs 57 over 20 cells);
# - the log / sqrt selection adds false calls where both options discover: R17 TP / FP over 20 cells, selection 60 / 8,
#   log alone 57 / 6, sqrt alone 45 / 3 (the per-taxon mixture BH'd jointly calls taxa neither option calls);
# - thinning each pair harder than rho 0.5 trims the null tail at little cost: R19.null hmp_tongue p < 0.01 log 1.88x
#   (rho 0.5) -> 1.74x (0.35) -> 1.60x (0.25); R17 TP / FP log 57 / 6 (0.5), 63 / 6 (0.35), 55 / 5 (0.25).
# Hence two arms, both on top of eperm_c5q and both only on the efull_nodet path (nothing else changes):
#   c6a: log alone (no selection) at the design's rho;   c6b: log alone at rho 0.35.
.efull_log_rho <- function(counts, meta, formula, tested_term, rho) {      # efull's log test, thinning factor given
  d <- .uk_design(meta, formula, tested_term); g1 <- d$g1; depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts)
  rows <- .with_seed(7L, if (nrow(counts) > 300L) sort(sample.int(nrow(counts), 300L)) else seq_len(nrow(counts)))
  cs <- counts[rows, , drop = FALSE]; cs <- cs[rowSums(cs > 0) >= 3, , drop = FALSE]
  med <- function(gm) { u <- .u_stat_cov(cs, depth, g1, NULL, gm, "log", rho = rho); stats::median(u$U / sqrt(u$v), na.rm = TRUE) }
  flo <- med(-2); fhi <- med(2)
  gm <- if (is.finite(flo) && is.finite(fhi) && sign(flo) != sign(fhi)) stats::uniroot(med, c(-2, 2), f.lower = flo, f.upper = fhi, tol = 2e-3)$root else 0
  u <- .u_stat_cov(counts, depth, g1, NULL, gm, "log", rho = rho); ok <- rowSums(counts > 0) >= 3 & u$v > 0
  xx <- stats::qlogis(pmin(pmax(rowMeans(counts > 0), 0.01), 0.99)); mm <- .squeeze(ifelse(ok, u$v, NA), u$df, x = xx)
  z <- ifelse(ok, stats::qnorm(stats::pt(u$U / sqrt(mm$v), mm$df)), NA)
  list(feature = rownames(counts), p = 2 * stats::pnorm(-abs(z)), est = u$U / log(2))
}
.c6_fit <- function(counts, meta, formula, tested_term, conf_rho = NA_real_) {
  v <- .epg_fit(counts, meta, formula, tested_term, paths = "cont", adet = TRUE, nodet_conf = TRUE, pk = "pos", minsw = 2L)
  if (identical(v$path, "efull_nodet")) {
    if (is.na(conf_rho)) v$p_full <- v$p_log else {
      w <- .efull_log_rho(counts, meta, formula, tested_term, conf_rho); i <- match(v$feature, w$feature); v$p_full <- w$p[i] }
    v$path <- "efull_conflog" }
  v
}
register_candidate("eperm_c6a", function(counts, meta, formula, tested_term) {
  v <- .c6_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_full, q = .storey_gq(v$p_full), estimate = v$est)
}, notes = "it49: eperm_c5q with log alone (no log / sqrt selection) under designed depth confounding")
register_candidate("eperm_c6b", function(counts, meta, formula, tested_term) {
  v <- .c6_fit(counts, meta, formula, tested_term, conf_rho = 0.35); data.frame(feature = v$feature, p = v$p_full, q = .storey_gq(v$p_full), estimate = v$est)
}, notes = "it49: eperm_c6a with each pair thinned to 0.35 of its common depth under designed depth confounding")
