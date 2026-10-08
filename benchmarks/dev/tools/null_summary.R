#!/usr/bin/env Rscript
# null_summary.R -- calibration on null-only dev-suite runs: per setting and candidate, P(any false discovery)
# (BH q 0.05; nominal <= 0.05 under a global null), mean false positives, and the null FPR at p < 0.05 with its
# binomial excess over 0.05 (one-sided p, cells pooled -- read it next to P(any), since false calls cluster).
#   Rscript benchmarks/dev/tools/null_summary.R benchmarks/dev/results/<label>/cells.csv [P(any) limit, default 0.075]
a <- commandArgs(TRUE); d <- utils::read.csv(a[1], stringsAsFactors = FALSE); lim <- if (length(a) > 1) as.numeric(a[2]) else 0.075
d <- d[d$n_pos == 0 & is.na(d$error), ]
s <- do.call(rbind, lapply(split(d, list(d$setting, d$candidate), drop = TRUE), function(x) data.frame(setting = x$setting[1], candidate = x$candidate[1],
  cells = nrow(x), p_any = mean(x$fp > 0), mean_fp = mean(x$fp), fpr05 = sum(x$null_p05) / sum(x$null_n),
  excess_p = stats::pbinom(sum(x$null_p05) - 1, sum(x$null_n), 0.05, lower.tail = FALSE))))
s <- s[order(s$setting, s$candidate), ]; print(format(s, digits = 3), row.names = FALSE)
o <- do.call(rbind, lapply(split(d, d$candidate), function(x) data.frame(candidate = x$candidate[1], cells = nrow(x), p_any = mean(x$fp > 0),
  worst_setting_p_any = max(tapply(x$fp > 0, x$setting, mean)), fpr05 = sum(x$null_p05) / sum(x$null_n))))
cat("\noverall:\n"); print(format(o, digits = 3), row.names = FALSE)
cat("\nP(any) <= ", lim, " overall: ", paste(o$candidate[o$p_any <= lim], collapse = ", "), "\n", sep = "")
