# -----------------------------------------------------------------------------
# regimes.R -- the Axis A regime grid (protocol section 4.3): a reference regime plus
# one-factor-at-a-time sweeps. make_regimes() is the single source of truth; the
# committed benchmarks/regimes.tsv is written from it.
# -----------------------------------------------------------------------------

reference_regime <- function() {
  data.frame(
    n_per_group = 50L, m = 500L, da_frac = 0.10, balance = "balanced", effect = "medium",
    signal_type = "mixed", depth_conf = 1, conf_phi = 0, exposure = "binary",
    design = "independent", bloom = "none", hetero = "no", stringsAsFactors = FALSE)
}

regime_sweeps <- function() list(
  n_per_group = c(20L, 50L, 100L, 200L),
  m           = c(200L, 500L, 1000L),
  da_frac     = c(0, 0.05, 0.10, 0.25, 0.40),
  balance     = c("balanced", "80_20", "100_0"),
  effect      = c("small", "medium", "large", "graded"),
  signal_type = c("abundance", "prevalence", "mixed"),
  depth_conf  = c(1, 2, 4, 9),
  conf_phi    = c(0, 0.4, 0.7),
  exposure    = c("binary", "continuous"),
  design      = c("independent", "repeated"),
  bloom       = c("none", "x20"),
  hetero      = c("no", "yes"))

make_regimes <- function(n_rep_alt = 20L, n_rep_null = 50L) {
  ref <- reference_regime(); sw <- regime_sweeps()
  rows <- list(cbind(regime_id = "R00", factor = "reference", level = "reference", ref, stringsAsFactors = FALSE))
  k <- 0L
  for (f in names(sw)) for (lv in sw[[f]]) {
    if (identical(as.character(lv), as.character(ref[[f]]))) next
    r <- ref; r[[f]] <- if (is.numeric(ref[[f]])) as.numeric(lv) else as.character(lv)
    if (is.integer(ref[[f]])) r[[f]] <- as.integer(lv)
    k <- k + 1L
    rows[[length(rows) + 1L]] <- cbind(regime_id = sprintf("R%02d", k), factor = f, level = as.character(lv), r, stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  out$is_null <- out$da_frac == 0
  out$n_rep <- ifelse(out$is_null, n_rep_null, n_rep_alt)
  rownames(out) <- NULL
  out
}

effect_sizes <- function(effect) switch(effect,
  small  = list(lfc2 = 0.5, logor = log(2),  graded = FALSE),
  medium = list(lfc2 = 1.0, logor = log(4),  graded = FALSE),
  large  = list(lfc2 = 2.0, logor = log(8),  graded = FALSE),
  graded = list(lfc2 = c(0.25, 2), logor = log(c(1.5, 8)), graded = TRUE))

if (sys.nframe() == 0L) {                 # Rscript regimes.R -> writes regimes.tsv beside benchmarks/
  out <- make_regimes()
  utils::write.table(out, file.path(dirname(dirname(dirname(normalizePath(sys.frames()[[1]]$ofile %||% ".")))), "regimes.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
}
