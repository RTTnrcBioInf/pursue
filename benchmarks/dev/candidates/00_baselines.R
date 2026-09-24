# Baselines every candidate is measured against. Not R&D candidates themselves.

# pursue02 and pursue02_presence read different arms of the same fit: memoise it per cell so the
# suite does not pay for PURSUE twice. (Each forked worker has its own copy; candidates within a
# cell run sequentially, so the last-input check is exact.)
.pursue_memo <- new.env()
.pursue_fit <- function(counts, meta, formula, tested_term) {
  if (!is.null(.pursue_memo$key) && identical(.pursue_memo$key, list(counts, meta, formula, tested_term))) return(.pursue_memo$val)
  v <- method_pursue(counts, meta, formula, tested_term)
  .pursue_memo$key <- list(counts, meta, formula, tested_term); .pursue_memo$val <- v; v
}
register_candidate("pursue02", function(counts, meta, formula, tested_term) {
  r <- .pursue_fit(counts, meta, formula, tested_term); r <- r[r$arm == "combined", ]
  data.frame(feature = r$feature, p = r$p)
}, notes = "PURSUE 0.2 as shipped: Cauchy combination of abundance and presence arms")

register_candidate("pursue02_presence", function(counts, meta, formula, tested_term) {
  r <- .pursue_fit(counts, meta, formula, tested_term); r <- r[r$arm == "presence", ]
  data.frame(feature = r$feature, p = r$p)
}, notes = "PURSUE 0.2 presence arm alone: stabilised ZINB, fixed depth offset, LRT on structural presence")

# Detection GLM; `depth` = none | covariate. LRT on the tested term's columns, any formula.
.detect_glm <- function(counts, meta, formula, tested_term, depth = c("none", "covariate")) {
  depth <- match.arg(depth)
  X <- stats::model.matrix(formula, meta); a <- attr(X, "assign")
  tc <- which(a == which(attr(stats::terms(formula), "term.labels") == tested_term))
  if (depth == "covariate") { ld <- log(meta$depth); X <- cbind(X, ld = ld - mean(ld)) }
  D <- counts > 0; n <- ncol(D)
  p <- vapply(seq_len(nrow(D)), function(j) {
    d <- as.integer(D[j, ]); if (sum(d) < 3 || sum(d) > n - 3) return(NA_real_)
    f1 <- tryCatch(suppressWarnings(stats::glm.fit(X, d, family = stats::binomial())), error = function(e) NULL)
    f0 <- tryCatch(suppressWarnings(stats::glm.fit(X[, -tc, drop = FALSE], d, family = stats::binomial())), error = function(e) NULL)
    if (is.null(f1) || is.null(f0)) return(NA_real_)
    stats::pchisq(max(f0$deviance - f1$deviance, 0), length(tc), lower.tail = FALSE)
  }, numeric(1))
  data.frame(feature = rownames(counts), p = p)
}
register_candidate("logistic", function(counts, meta, formula, tested_term)
  .detect_glm(counts, meta, formula, tested_term, "none"),
  notes = "detected ~ design; no depth. Breaks under depth confounding on msq/mid")
register_candidate("logistic_depth", function(counts, meta, formula, tested_term)
  .detect_glm(counts, meta, formula, tested_term, "covariate"),
  notes = "detected ~ design + centred log-depth. Calibrated on house/msq/mid; no power under strong confounding")

register_candidate("lm_logtss", function(counts, meta, formula, tested_term) {
  r <- method_lm_logtss(counts, meta, formula, tested_term); data.frame(feature = r$feature, p = r$p)
}, notes = "benchmark elementary comparator: LM on log(TSS + pseudocount)")
