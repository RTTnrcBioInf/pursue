#!/usr/bin/env Rscript
# Expand the benchmark design into task lists, one line per cell, consumed by the SLURM
# array scripts (line number = SLURM_ARRAY_TASK_ID).
#   Rscript hpc/make_tasklist.R [--pool evaluation|tuning|all] [--simulators house,msq,mid,sd2,sps] [--methods all]
suppressMessages(library(optparse))
op <- OptionParser(option_list = list(
  make_option("--pool", default = "evaluation"), make_option("--simulators", default = "house,msq,mid,sd2,sps"),
  make_option("--methods", default = "all"), make_option("--out", default = "hpc/tasks"),
  make_option("--regimes", default = "all", help = "comma-separated regime ids for Axis A (default all)"),
  make_option("--max-rep", type = "integer", default = 0L, help = "cap replicates per cell (0 = protocol values)")))
opt <- parse_args(op)
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), ".."))
bench <- file.path(root, "benchmarks"); dir.create(file.path(root, opt$out), recursive = TRUE, showWarnings = FALSE)
tpl <- read.delim(file.path(bench, "templates.tsv"), comment.char = "#", stringsAsFactors = FALSE)
pools <- if (opt$pool == "all") c("evaluation", "tuning") else opt$pool
tA <- tpl$id[tpl$pool %in% pools]
reg <- read.delim(file.path(bench, "regimes.tsv"), stringsAsFactors = FALSE)
if (opt$regimes != "all") reg <- reg[reg$regime_id %in% strsplit(opt$regimes, ",")[[1]], ]
if (opt$`max-rep` > 0) reg$n_rep <- pmin(reg$n_rep, opt$`max-rep`)
sims <- strsplit(opt$simulators, ",")[[1]]
line <- function(axis, sim, t, r, rep) sprintf("--axis %s --simulator %s --template %s --regime %s --replicate %d --methods %s", axis, sim, t, r, rep, opt$methods)
A <- unlist(lapply(sims, function(s) lapply(tA, function(t) lapply(seq_len(nrow(reg)), function(i) sapply(seq_len(reg$n_rep[i]), function(k) line("A", s, t, reg$regime_id[i], k))))))
writeLines(A, file.path(root, opt$out, "axisA.txt"))
sp <- read.delim(file.path(bench, "implant_specs.tsv"), comment.char = "#", stringsAsFactors = FALSE)
if (opt$`max-rep` > 0) sp$n_rep <- pmin(sp$n_rep, opt$`max-rep`)
B <- unlist(lapply(tA, function(t) lapply(seq_len(nrow(sp)), function(i) sapply(seq_len(sp$n_rep[i]), function(k) line("B", "implant", t, sp$spec_id[i], k)))))
writeLines(B, file.path(root, opt$out, "axisB.txt"))
tC <- tpl$id[tpl$pool %in% c(pools, "axisE", "axisD")]
nC_rep <- if (opt$`max-rep` > 0) min(50L, opt$`max-rep`) else 50L
C <- unlist(lapply(tC, function(t) lapply(c("full", "stratified", "depth_quintile", "cluster"), function(s) sapply(seq_len(nC_rep), function(k) line("C", "nullperm", t, s, k)))))
writeLines(C, file.path(root, opt$out, "axisC.txt"))
R <- unlist(lapply(sims, function(s) sapply(tA, function(t) sprintf("%s %s", s, t))))
writeLines(R, file.path(root, opt$out, "realism.txt"))
cat(sprintf("axis A: %d cells | axis B: %d | axis C: %d | realism: %d\n", length(A), length(B), length(C), length(R)))
cat("submit with: sbatch --array=1-N%200 hpc/slurm/axisA.sbatch  (N from above; %200 caps concurrent tasks)\n")
