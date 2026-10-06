#!/usr/bin/env Rscript
# merge_DE.R -- per-(task, method) outputs of run_DE.sh -> one table per task plus a summary.
#   Rscript hpc/merge_DE.R results/axisDE results/summary/axisDE
a <- commandArgs(TRUE); src <- a[1]; dst <- a[2]; dir.create(dst, recursive = TRUE, showWarnings = FALSE)
f <- list.files(src, pattern = "^[DE]__.*\\.csv$", full.names = TRUE)
parts <- strsplit(sub("\\.csv$", "", basename(f)), "__")
task <- vapply(parts, function(p) paste(p[1], p[2], p[length(p)], sep = "__"), "")
options(width = 200)
for (t in sort(unique(task))) {
  x <- do.call(rbind, lapply(f[task == t], read.csv, stringsAsFactors = FALSE))
  write.csv(x, file.path(dst, paste0(t, ".csv")), row.names = FALSE)
  cat("\n==", t, "==\n")
  if (grepl("biotruth$", t)) {
    cat("contrast:", unique(x$contrast), "\n")
    k <- intersect(c("method", "arm", "n_tested", "n_calls", "n_calls_annotated", "frac_expected_direction", "enrichment_or", "enrichment_p", "runtime_s"), names(x))
    y <- x[order(-x$n_calls), k]; num <- vapply(y, is.numeric, TRUE); y[num] <- lapply(y[num], signif, 3); print(y, row.names = FALSE)
  } else if (grepl("spikein$", t)) {
    y <- do.call(rbind, lapply(split(x, list(x$method, x$arm), drop = TRUE), function(d) data.frame(method = d$method[1], arm = d$arm[1],
      groupings = nrow(d), mean_calls = mean(d$n_calls), spikein_calls = sum(d$spikein_calls), groupings_with_spikein_call = sum(d$spikein_calls > 0))))
    print(y[order(y$spikein_calls, -y$mean_calls), ], row.names = FALSE)
  } else {
    y <- do.call(rbind, lapply(split(x, list(x$method, x$arm), drop = TRUE), function(d) data.frame(method = d$method[1], arm = d$arm[1],
      splits = nrow(d), nhits = mean(d$nhits), replication_pct = mean(d$replication_pct, na.rm = TRUE), conflict_pct = mean(d$conflict_pct, na.rm = TRUE))))
    num <- vapply(y, is.numeric, TRUE); y[num] <- lapply(y[num], signif, 3); print(y[order(-y$replication_pct), ], row.names = FALSE)
  }
}
