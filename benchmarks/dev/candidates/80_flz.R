# --- it55: covariate designs through Freedman-Lane with the plain partial-correlation statistic (Phase 2, target 1) -----
# srv18 rejected Freedman-Lane for covariates (R21 +0.2, B_conf07 16.30 -> 14.85) -- with the HC0 t, which it54 showed
# costs power whenever the exposure effect is multiplicative. Re-try: residual scores (on 1 + Z) against x residualised
# on 1 + Z, plain correlation, permuted Freedman-Lane style; c10 everywhere else. Freedman-Lane is approximately exact.
# Result (bank3 B_conf07, TP / false calls per cell): hmp_tongue 5.8 / 0.6 (c5q 5.8 / 0.8, logistic presence 5.6 / 0.4),
# twinsuk 33.0 / 5.0 (c5q 30.0 / 3.0, logistic 32.8 / 3.0) -- one cell with 16 false calls (c5q 6). B_conf07's
# confounder also moves 10% of null taxa multiplicatively; a linear adjustment of the scores removes it only
# approximately, and the plain statistic no longer discounts the samples where that residual is largest. +10% power for
# a higher FDR: not promoted. The covariate designs stay on efull.
register_candidate("eperm_c11", function(counts, meta, formula, tested_term) {
  .with_opts({ v <- .epg_fit(counts, meta, formula, tested_term, paths = c("cont", "cov"), adet = TRUE, nodet_conf = TRUE, pk = "pos", minsw = 2L)
    data.frame(feature = v$feature, p = v$p_full, q = .storey_gq(v$p_full), estimate = v$est) }, pursue.flstat = "corz")
}, notes = "it55: eperm_c10 + covariate designs through Freedman-Lane with the partial-correlation statistic")
