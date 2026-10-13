#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# spread_suite.R -- spread-only implants on the tuning templates (dev/tools/spread.R): does a method call taxa whose
# location is unchanged when one group (or one end of a continuous exposure) is more variable?
#
#   Rscript benchmarks/dev/spread_suite.R --candidates eperm_c10,ref_zicoseq --settings H_unb,H_cont --reps 20 \
#           --cores 24 --label srv27
#
# Same candidates, templates, prevalence filter and dev master seed as devsuite.R. There is no location signal, so
# every BH q <= 0.05 call is false under the benchmark's truth. Per cell and candidate: calls on the altered taxa and on
# the untouched ones, and p < 0.05 counts on each (FPR). Checkpointed per cell; summary.csv / summary.txt at the end.
# -----------------------------------------------------------------------------
suppressMessages(library(optparse))
op <- OptionParser(option_list = list(
  make_option("--candidates", default = "eperm_c10"),
  make_option("--settings", default = "H_unb,H_bal,H_all,H_rd,H_cont"),
  make_option("--templates", default = "hmp_tongue,twinsuk_stool"),
  make_option("--reps", type = "integer", default = 5L),
  make_option("--cores", type = "integer", default = 2L),
  make_option("--label", default = format(Sys.time(), "%Y%m%d-%H%M")),
  make_option("--seed", type = "integer", default = 9001L)))
opt <- parse_args(op)
root <- normalizePath(Sys.getenv("PURSUE_BENCH_ROOT", unset = "benchmarks"), mustWork = TRUE)
Sys.setenv(PURSUE_BENCH_ROOT = root)
if (!nzchar(Sys.getenv("PURSUE_DATA_ROOT"))) Sys.setenv(PURSUE_DATA_ROOT = file.path(root, "data"))
for (f in c("R/engine/templates.R", "R/engine/regimes.R", "R/engine/metrics.R", "R/simulators/sim_house.R",
            "R/simulators/implant.R", "R/simulators/dispatch.R", "R/methods/elementary.R"))
  source(file.path(root, f))
source(file.path(root, "dev", "tools", "spread.R"))

.cands <- new.env()
register_candidate <- function(name, fn, notes = "") assign(name, list(fn = fn, notes = notes), envir = .cands)
for (f in sort(list.files(file.path(root, "dev", "candidates"), pattern = "\\.R$", full.names = TRUE))) source(f)
# two cheap location references beside 90_comparators.R's: Wilcoxon on proportions and a t-test on log proportions
local({ e <- new.env(); sys.source(file.path(root, "R/methods/elementary.R"), envir = e)
  for (m in c("wilcoxon_tss", "lm_logtss")) local({ fn <- paste0("method_", m)
    register_candidate(paste0("ref_", m), function(counts, meta, formula, tested_term) {
      r <- e[[fn]](counts, meta, formula, tested_term); r <- r[r$arm %in% c("single", "combined"), ]
      if (all(!is.finite(r$p))) stop("no p-values: ", paste(unique(r$status), collapse = ","))
      data.frame(feature = r$feature, p = r$p, q = r$q, estimate = r$estimate) }) }) })
cand <- strsplit(opt$candidates, ",")[[1]]
miss <- setdiff(cand, ls(.cands)); if (length(miss)) stop("unknown candidates: ", paste(miss, collapse = ", "))
sets <- strsplit(opt$settings, ",")[[1]]
miss <- setdiff(sets, names(SPREAD_SPECS)); if (length(miss)) stop("unknown settings: ", paste(miss, collapse = ", "))
tpls <- strsplit(opt$templates, ",")[[1]]
TPL <- setNames(lapply(tpls, function(id) readRDS(file.path(root, "devdata", paste0(id, ".rds")))), tpls)
cell_seed <- function(...) { k <- paste(..., sep = "|"); as.integer((opt$seed * 1e5 + sum(utf8ToInt(k) * seq_along(utf8ToInt(k)))) %% .Machine$integer.max) }

jobs <- expand.grid(set = sets, tp = tpls, rep = seq_len(opt$reps), stringsAsFactors = FALSE)
cat(sprintf("spread suite '%s': %d candidates x %d settings x %d templates x %d reps = %d cells, %d cores\n",
            opt$label, length(cand), length(sets), length(tpls), opt$reps, nrow(jobs), opt$cores))
odir <- file.path(root, "dev", "results", opt$label); cdir <- file.path(odir, "cells")
dir.create(cdir, recursive = TRUE, showWarnings = FALSE); progress <- file.path(odir, "progress.log")

run_job <- function(k) {
  st <- jobs$set[k]; tid <- jobs$tp[k]; rep <- jobs$rep[k]
  ck <- file.path(cdir, sprintf("%s__%s__r%02d.csv", st, tid, rep))
  if (file.exists(ck)) { prev <- tryCatch(utils::read.csv(ck, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(prev) && all(cand %in% prev$candidate)) return(prev[prev$candidate %in% cand, ]) }
  sim <- implant_spread(TPL[[tid]], SPREAD_SPECS[[st]], cell_seed("spread", st, tid, rep))
  if (!is.null(sim$unsupported)) return(NULL)
  keep <- rowMeans(sim$counts > 0) >= 0.10 & rowSums(sim$counts > 0) >= 3
  ct <- sim$counts[keep, , drop = FALSE]; sp <- sim$truth$spread[match(rownames(ct), sim$truth$feature)]
  meta <- sim$meta[colnames(ct), , drop = FALSE]
  out <- lapply(cand, function(cn) {
    t0 <- proc.time()[["elapsed"]]
    r <- tryCatch(.cands[[cn]]$fn(ct, meta, sim$formula, sim$tested_term), error = function(e) e)
    el <- proc.time()[["elapsed"]] - t0
    base <- data.frame(setting = st, template = tid, rep = rep, candidate = cn, m = nrow(ct), n_sp = sum(sp == 1L),
                       mean_shift = sim$loc_diag[["mean_shift"]], median_shift = sim$loc_diag[["median_shift"]],
                       abs_mean_shift = sim$loc_diag[["abs_mean_shift"]], stringsAsFactors = FALSE)
    if (inherits(r, "error")) return(cbind(base, tested = 0, rej = NA, rej_sp = NA, rej_un = NA, p05_sp = NA, p05_un = NA,
                                           nt_sp = NA, nt_un = NA, secs = el, error = substr(conditionMessage(r), 1, 100)))
    p <- r$p[match(rownames(ct), r$feature)]; ok <- is.finite(p); q <- rep(NA_real_, length(p))
    if (is.null(r$q)) q[ok] <- stats::p.adjust(p[ok], "BH") else { qo <- r$q[match(rownames(ct), r$feature)]; q[ok] <- ifelse(is.finite(qo[ok]), qo[ok], 1) }
    rj <- !is.na(q) & q <= 0.05
    cbind(base, tested = sum(ok), rej = sum(rj), rej_sp = sum(rj & sp == 1L), rej_un = sum(rj & sp == 0L),
          p05_sp = sum(ok & sp == 1L & p < 0.05), p05_un = sum(ok & sp == 0L & p < 0.05),
          nt_sp = sum(ok & sp == 1L), nt_un = sum(ok & sp == 0L), secs = el, error = NA_character_)
  })
  out <- do.call(rbind, out)
  tmp <- paste0(ck, ".tmp", Sys.getpid()); utils::write.csv(out, tmp, row.names = FALSE); file.rename(tmp, ck)
  cat(sprintf("%s %s %s r%d\n", format(Sys.time(), "%H:%M:%S"), st, tid, rep), file = progress, append = TRUE)
  out
}
t_start <- proc.time()[["elapsed"]]
res <- parallel::mclapply(seq_len(nrow(jobs)), run_job, mc.cores = opt$cores, mc.preschedule = FALSE)
D <- do.call(rbind, Filter(function(x) is.data.frame(x), res))
write.csv(D, file.path(odir, "cells.csv"), row.names = FALSE)

summ <- function(D) {
  D <- D[is.na(D$error), ]
  s <- do.call(rbind, lapply(split(D, list(D$setting, D$candidate), drop = TRUE), function(x) data.frame(
    setting = x$setting[1], candidate = x$candidate[1], cells = nrow(x), p_any = mean(x$rej > 0), calls = mean(x$rej),
    calls_spread = mean(x$rej_sp), calls_other = mean(x$rej_un), fpr05_spread = sum(x$p05_sp) / sum(x$nt_sp),
    fpr05_other = sum(x$p05_un) / sum(x$nt_un), secs = median(x$secs), stringsAsFactors = FALSE)))
  s[order(s$setting, -s$fpr05_spread), ] }
sm <- summ(D); write.csv(sm, file.path(odir, "summary.csv"), row.names = FALSE)
dg <- aggregate(cbind(m, n_sp, mean_shift, median_shift, abs_mean_shift) ~ setting + template, D[!duplicated(paste(D$setting, D$template, D$rep)), ], mean)
options(width = 200); sm$secs <- round(sm$secs, 1)
sink(file.path(odir, "summary.txt"), split = TRUE)
cat(sprintf("%d cells in %.1f min. No location signal: every call is false. fpr05 = share of taxa at p < 0.05 (nominal 0.05).\n\n",
            nrow(jobs), (proc.time()[["elapsed"]] - t_start) / 60))
print(format(sm, digits = 3), row.names = FALSE)
cat("\nlocation diagnostic on the altered taxa (absolute scale, sd units; binary: case - control, continuous: correlation with z):\n")
print(format(dg, digits = 3), row.names = FALSE)
sink()
