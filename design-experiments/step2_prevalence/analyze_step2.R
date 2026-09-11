#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
f <- if (length(args) >= 1) args[1] else "step2_results.csv"
r <- read.csv(f, stringsAsFactors = FALSE)
prev_methods  <- c("naive", "depthcov", "er_lm", "occ_depth", "occ_full", "occ_cal", "occ_trunc", "occ_ztnb", "zip")
abund_methods <- c("abund_lm", "abund_lm_depth", "abund_trunc", "abund_ztnb", "lm_pseudo", "lm_pseudo_depth", "lm_pseudo_mc")
wide <- function(d, val) {
  a <- aggregate(as.formula(paste(val, "~ method + c_depth + n")), data = d, FUN = function(x) mean(x, na.rm = TRUE))
  t <- reshape(a, idvar = c("c_depth", "n"), timevar = "method", direction = "wide")
  names(t) <- sub(paste0("^", val, "\\."), "", names(t)); t[order(t$c_depth, t$n), ]
}
show <- function(title, d, val) { cat("\n====", title, "====\n"); print(wide(d, val), digits = 2, row.names = FALSE) }

pv <- subset(r, method %in% prev_methods); pv$method <- factor(pv$method, levels = prev_methods)
cat("################ PREVALENCE ARM (test of structural presence) ################\n")
show("Type I on NULL taxa (p<0.05)",                pv, "rej_null")
show("ABUNDANCE SHADOW: rejection on abund_only taxa (should be ~0.05)", pv, "rej_abund_only")
show("Power on prev_only taxa (p<0.05)",            pv, "rej_prev_only")
show("Power on both",                               pv, "rej_both")
show("FDR at BH 0.05 (truth = prev_only|both)",     pv, "fdr")
show("TPR at BH 0.05",                              pv, "tpr")
show("n_fail (of n_taxa)",                          pv, "n_fail")

ab <- subset(r, method %in% abund_methods); ab$method <- factor(ab$method, levels = abund_methods)
cat("\n\n################ ABUNDANCE ARM (test of abundance-when-present) ################\n")
show("Type I on NULL taxa (p<0.05)",                ab, "rej_null")
show("PREVALENCE SHADOW: rejection on prev_only taxa (should be ~0.05)", ab, "rej_prev_only")
show("Power on abund_only taxa (p<0.05)",           ab, "rej_abund_only")
show("FDR at BH 0.05 (truth = abund_only|both)",    ab, "fdr")
show("TPR at BH 0.05",                              ab, "tpr")

cat("\n==== by dispersion (c_depth=4, n=100): type I ====\n")
a <- aggregate(rej_null ~ method + disp, data = subset(r, c_depth == 4 & n == 100), FUN = mean)
print(reshape(a, idvar = "method", timevar = "disp", direction = "wide"), digits = 2, row.names = FALSE)
