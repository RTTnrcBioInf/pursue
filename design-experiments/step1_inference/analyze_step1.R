#!/usr/bin/env Rscript
# Aggregate step-1 results into decision tables.
args <- commandArgs(trailingOnly = TRUE)
f <- if (length(args) >= 1) args[1] else "step1_results.csv"
r <- read.csv(f, stringsAsFactors = FALSE)

mse <- function(x) sd(x) / sqrt(length(x))
fmt <- function(m, s) sprintf("%.3f±%.3f", m, s)

order_methods <- c("ols", "modt", "fl_199", "fl_999", "fl_9999",
                   "pvw_199", "pvw_999", "pvt_199", "pvt_999", "mm_199", "mm_999")
r$method <- factor(r$method, levels = intersect(order_methods, unique(r$method)))

cat("\n==================== A. NULL CELLS: calibration (FPR at 0.05, mean ± MC se) ====================\n")
nul <- subset(r, cell == "null")
a <- aggregate(fpr05 ~ method + error + n, data = nul, FUN = mean)
tab <- reshape(a, idvar = c("error", "n"), timevar = "method", direction = "wide")
names(tab) <- sub("^fpr05\\.", "", names(tab))
tab <- tab[order(tab$error, tab$n), ]
print(tab, digits = 3, row.names = FALSE)

cat("\n==================== A2. NULL CELLS: BH rejections at q=0.05 per 500 taxa (should be ~0) ====================\n")
a <- aggregate(n_rej ~ method + error + n, data = nul, FUN = mean)
tab <- reshape(a, idvar = c("error", "n"), timevar = "method", direction = "wide")
names(tab) <- sub("^n_rej\\.", "", names(tab))
print(tab[order(tab$error, tab$n), ], digits = 2, row.names = FALSE)

cat("\n==================== B. ALT CELLS: observed FDR at BH 0.05 ====================\n")
alt <- subset(r, cell == "alt")
a <- aggregate(fdr ~ method + error + n, data = alt, FUN = mean)
tab <- reshape(a, idvar = c("error", "n"), timevar = "method", direction = "wide")
names(tab) <- sub("^fdr\\.", "", names(tab))
print(tab[order(tab$error, tab$n), ], digits = 3, row.names = FALSE)

cat("\n==================== C. ALT CELLS: TPR at BH 0.05 ====================\n")
a <- aggregate(tpr ~ method + error + n, data = alt, FUN = mean)
tab <- reshape(a, idvar = c("error", "n"), timevar = "method", direction = "wide")
names(tab) <- sub("^tpr\\.", "", names(tab))
print(tab[order(tab$error, tab$n), ], digits = 3, row.names = FALSE)

cat("\n==================== D. Resolution: min attainable p (alt cells, median over reps) ====================\n")
a <- aggregate(p_min ~ method + n, data = alt, FUN = median)
tab <- reshape(a, idvar = "n", timevar = "method", direction = "wide")
names(tab) <- sub("^p_min\\.", "", names(tab))
print(format(tab, digits = 2, scientific = TRUE), row.names = FALSE)

cat("\n==================== E. Summary across error models (alt cells): mean FDR / mean TPR ====================\n")
a <- aggregate(cbind(fdr, tpr) ~ method + n, data = alt, FUN = mean)
a$fdr <- round(a$fdr, 3); a$tpr <- round(a$tpr, 3)
print(reshape(a, idvar = "n", timevar = "method", direction = "wide"), row.names = FALSE)

cat("\n==================== F. Worst-case FDR over error models (alt cells) ====================\n")
w <- aggregate(fdr ~ method + error + n, data = alt, FUN = mean)
w <- aggregate(fdr ~ method + n, data = w, FUN = max)
print(reshape(w, idvar = "n", timevar = "method", direction = "wide"), digits = 3, row.names = FALSE)

cat("\nlimma prior df (d0) by scenario, null cells:\n")
print(aggregate(d0 ~ error + n, data = subset(nul, method == "modt"), FUN = median), digits = 3, row.names = FALSE)
