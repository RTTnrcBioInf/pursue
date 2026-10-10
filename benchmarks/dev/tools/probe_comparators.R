#!/usr/bin/env Rscript
# probe_comparators.R -- what makes LDM and ZicoSeq stronger than PURSUE on msq / sd2 (bank1)? Two outputs:
#   1. the installed sources of LDM and GUniFrac (ZicoSeq), deparsed function by function, for reading;
#   2. LDM on every banked cell with each scale's p, F and q kept apart (p.otu.freq, p.otu.tran, p.otu.omni),
#      exactly as the benchmark wrapper calls it, plus ZicoSeq's raw and adjusted p.
#   Rscript benchmarks/dev/tools/probe_comparators.R benchmarks/dev/bank/bank1 benchmarks/dev/results/probe1
a <- commandArgs(TRUE); bank <- a[1]; out <- a[2]; dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (pkg in c("LDM", "GUniFrac")) {
  if (!requireNamespace(pkg, quietly = TRUE)) { cat(pkg, "not installed\n"); next }
  ns <- asNamespace(pkg); fns <- sort(ls(ns, all.names = TRUE))
  txt <- unlist(lapply(fns, function(f) { x <- get(f, envir = ns); if (is.function(x)) c(paste0("#### ", pkg, "::", f), deparse(x), "") }))
  writeLines(c(paste("#", pkg, as.character(utils::packageVersion(pkg))), txt), file.path(out, paste0(pkg, "_source.R.txt")))
  cat(pkg, as.character(utils::packageVersion(pkg)), length(fns), "objects\n")
}
fs <- list.files(bank, pattern = "^(sd2|msq).*\\.rds$", full.names = TRUE)
res <- parallel::mclapply(fs, function(f) {
  b <- readRDS(f); fo <- if (is.character(b$formula)) stats::as.formula(b$formula) else b$formula
  y <- t(b$counts); storage.mode(y) <- "numeric"; assign("y", y, envir = globalenv())
  o <- tryCatch(suppressMessages(LDM::ldm(stats::as.formula(paste("y ~", b$tested_term)), data = b$meta, fdr.nominal = 0.05, seed = 1, n.perm.max = 20000, verbose = FALSE)),
                error = function(e) e)
  if (inherits(o, "error")) return(list(file = basename(f), error = conditionMessage(o)))
  keep <- grep("otu|global", names(o), value = TRUE)
  z <- tryCatch(suppressMessages(GUniFrac::ZicoSeq(meta.dat = b$meta, feature.dat = b$counts, grp.name = b$tested_term, feature.dat.type = "count", prev.filter = 0,
         mean.abund.filter = 0, max.abund.filter = 0, min.prop = 0, is.winsor = TRUE, outlier.pct = 0.03, is.post.sample = TRUE, post.sample.no = 25,
         link.func = list(function(x) sign(x) * sqrt(abs(x))), stats.combine.func = max, perm.no = 1999, ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
         is.fwer = FALSE, verbose = FALSE, return.feature.dat = TRUE)), error = function(e) e)
  small <- function(L) L[vapply(L, function(x) utils::object.size(x) < 2e5, logical(1))]       # keep the result small (no permutation matrices)
  zk <- if (inherits(z, "error")) list(error = conditionMessage(z)) else small(z[setdiff(names(z), "feature.dat")])
  list(file = basename(f), features = rownames(b$counts), ldm = small(o[keep]), ldm_names = names(o), zico = zk, zico_names = if (inherits(z, "error")) NULL else names(z))
}, mc.cores = as.integer(Sys.getenv("CORES", "8")), mc.preschedule = FALSE)
saveRDS(res, file.path(out, "probe.rds"), compress = "xz")
cat(sum(vapply(res, function(r) is.null(r$error), logical(1))), "of", length(res), "cells probed ->", out, "\n")
