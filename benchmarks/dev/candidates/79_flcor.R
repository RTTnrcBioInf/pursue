# --- it54: the continuous-exposure permutation statistic -- plain correlation, not the HC0 t (Phase 2, target 2: B_cont) ---
# bank3, twinsuk B_cont r04 (m 2212): c5q 0 TP; ZicoSeq, LDM 31. The sqrt kernel's scores are as informative as
# ZicoSeq's input (t vs a plain regression of sqrt RA on exposure: r 0.98), and its permutation null is right (sd 1.03).
# What loses is the studentisation: under a multiplicative exposure effect the residual spread grows with |x|, so the
# HC0 standard error of a real association is inflated by the high-leverage samples, while under permutation it is
# not -- the observed t shrinks relative to its null. On that cell, 20 000 permutations of the same scores: plain
# correlation 20 TP / 1 FP, HC0 t 0 / 0; null taxa's FPR at p < 0.01 identical (0.0136 both).
# With no covariates the permutation of a plain correlation is exact under exchangeability, and the compositional
# centre (thinning by exp(gamma x)) restores exchangeability of null taxa under a compositional shift, so nothing of
# validity is lost; the HC0 t stays where covariates enter (Freedman-Lane is then only approximate anyway).
# options(pursue.flstat = "cor"), read by .fl_stat in 70_epg.R; memo caches cleared around the call.
.with_opts <- function(expr, ...) { old <- options(...); .clear_fit_memos(); on.exit({ options(old); .clear_fit_memos() }); expr }
register_candidate("eperm_c10", function(counts, meta, formula, tested_term) {
  .with_opts({ v <- .epg_fit(counts, meta, formula, tested_term, paths = "cont", adet = TRUE, nodet_conf = TRUE, pk = "pos", minsw = 2L)
    data.frame(feature = v$feature, p = v$p_full, q = .storey_gq(v$p_full), estimate = v$est) }, pursue.flstat = "cor")
}, notes = "it54: eperm_c5q with a plain-correlation permutation statistic for continuous exposure without covariates")
