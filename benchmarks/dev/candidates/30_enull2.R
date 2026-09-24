# 30_enull2.R -- iteration 3: the estimated-scale empirical null (05_enull.R) on top of existing fits.
# Shares their memoised fits, so these cost almost nothing beyond the underlying candidate.
.en_wrap <- function(v, b, se) { e <- enull2(b, se); data.frame(feature = v$feature, p = e$p, estimate = b - e$delta) }
register_candidate("sharedlink_en", function(counts, meta, formula, tested_term) {
  v <- .sl_fit(counts, meta, formula, tested_term); .en_wrap(v, v$b, v$se)
}, notes = "it3: shared learned link, Wald test against an empirical null with estimated centre AND scale")
register_candidate("nbglm_en", function(counts, meta, formula, tested_term) {
  v <- .nb_fit(counts, meta, formula, tested_term); .en_wrap(v, v$b, v$se)
}, notes = "it3: NB GLM model-based SE, empirical null with estimated centre and scale")
register_candidate("nbglm_hc3_en", function(counts, meta, formula, tested_term) {
  v <- .nb_fit(counts, meta, formula, tested_term); .en_wrap(v, v$b, v$seh)
}, notes = "it3: NB GLM sandwich SE, empirical null with estimated centre and scale")
