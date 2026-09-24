#!/usr/bin/env Rscript
# compare.R -- rank candidates from several dev-suite runs together. Cells are identical across
# runs (same dev seed, setting, template, replicate), so results merge without re-running.
#   Rscript benchmarks/dev/compare.R --labels it1,it2 [--candidates a,b,c] [--brief]
# Where a candidate appears in several runs, the LATER label wins.
suppressMessages(library(optparse))
op <- OptionParser(option_list = list(make_option("--labels"), make_option("--candidates", default = "all"),
  make_option("--brief", action = "store_true", default = FALSE), make_option("--bench-root", default = NULL)))
opt <- parse_args(op)
root <- normalizePath(if (!is.null(opt$`bench-root`)) opt$`bench-root` else Sys.getenv("PURSUE_BENCH_ROOT", unset = "benchmarks"))
source(file.path(root, "dev", "score.R"))
labs <- strsplit(opt$labels, ",")[[1]]
D <- NULL
for (l in labs) {
  f <- file.path(root, "dev", "results", l, "cells.csv")
  if (file.exists(f)) x <- utils::read.csv(f, stringsAsFactors = FALSE) else {
    # run still in progress (or interrupted): read its per-cell checkpoints instead
    ck <- list.files(file.path(root, "dev", "results", l, "cells"), pattern = "\\.csv$", full.names = TRUE)
    if (!length(ck)) stop("no results for label ", l)
    x <- do.call(rbind, lapply(ck, utils::read.csv, stringsAsFactors = FALSE))
    cat(sprintf("[%s: in progress, %d cells checkpointed]\n", l, length(ck)))
  }
  # later label wins per (candidate, setting): runs on different settings of the same candidate merge
  if (!is.null(D)) D <- D[!(paste(D$candidate, D$setting) %in% unique(paste(x$candidate, x$setting))), ]
  D <- rbind(D, x)
}
if (opt$candidates != "all") D <- D[D$candidate %in% strsplit(opt$candidates, ",")[[1]], ]
sm <- summarise_cells(D); sb <- scoreboard(sm)
cat(sprintf("runs: %s | %d candidates | %d cell-candidate rows\n\n", paste(labs, collapse = ", "), length(unique(D$candidate)), nrow(D)))
print_report(sm, sb, per_setting = !opt$brief)
