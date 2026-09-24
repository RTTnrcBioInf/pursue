#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# devsuite.R -- the inner loop of PURSUE R&D (claude/rd-charter.md).
#
# Scores candidate DA procedures on the TUNING templates only (hmp_tongue, twinsuk_stool, in
# benchmarks/devdata/), never on the evaluation templates the benchmark reports, with seeds
# disjoint from the benchmark's. Every candidate sees exactly the same cells.
#
#   Rscript benchmarks/dev/devsuite.R --candidates pursue02,logistic_depth --sims house,implant \
#           --reps 5 --cores 2 --label baseline
#
# Candidates live in benchmarks/dev/candidates/*.R and call
#   register_candidate("name", function(counts, meta, formula, tested_term) ..., notes = "...")
# `counts` is features x samples, already prevalence-filtered exactly as the benchmark does;
# `meta$depth` is each sample's library size. Return data.frame(feature, p [, estimate]):
# one row per feature, p the candidate's single DA call. Architecture is unconstrained.
#
# Scoring (charter): calibration is a hard gate, then power. Against truth_abs, per cell and
# candidate: rejections at BH q <= 0.05, true / false positives, and the false positive rate
# among truly null features at p < 0.05 (the calibration measure that does not depend on power).
# -----------------------------------------------------------------------------
suppressMessages(library(optparse))
op <- OptionParser(option_list = list(
  make_option("--candidates", default = "all"),
  make_option("--sims", default = "house,implant", help = "house,implant,msq,mid"),
  make_option("--templates", default = "hmp_tongue,twinsuk_stool"),
  make_option("--settings", default = "all", help = "comma-separated setting ids, or all"),
  make_option("--reps", type = "integer", default = 5L),
  make_option("--cores", type = "integer", default = 2L),
  make_option("--label", default = format(Sys.time(), "%Y%m%d-%H%M")),
  make_option("--seed", type = "integer", default = 9001L, help = "dev master seed; the benchmark uses 1"),
  make_option("--bench-root", default = NULL)))
opt <- parse_args(op)
if (!is.null(opt$`bench-root`)) Sys.setenv(PURSUE_BENCH_ROOT = opt$`bench-root`)
root <- normalizePath(Sys.getenv("PURSUE_BENCH_ROOT", unset = "benchmarks"), mustWork = TRUE)
Sys.setenv(PURSUE_BENCH_ROOT = root)
if (!nzchar(Sys.getenv("PURSUE_DATA_ROOT"))) Sys.setenv(PURSUE_DATA_ROOT = file.path(root, "data"))
for (f in c("R/engine/templates.R", "R/engine/regimes.R", "R/engine/metrics.R", "R/simulators/sim_house.R",
            "R/simulators/implant.R", "R/simulators/dispatch.R", "R/methods/elementary.R"))
  source(file.path(root, f))
cache_dir <- file.path(dirname(root), "cache")

# ---- candidates ----------------------------------------------------------------------------
.cands <- new.env()
register_candidate <- function(name, fn, notes = "") assign(name, list(fn = fn, notes = notes), envir = .cands)
for (f in sort(list.files(file.path(root, "dev", "candidates"), pattern = "\\.R$", full.names = TRUE))) source(f)
cand <- if (opt$candidates == "all") sort(ls(.cands)) else strsplit(opt$candidates, ",")[[1]]
miss <- setdiff(cand, ls(.cands)); if (length(miss)) stop("unknown candidates: ", paste(miss, collapse = ", "))

# ---- settings ------------------------------------------------------------------------------
regimes <- read.delim(file.path(root, "regimes.tsv"), stringsAsFactors = FALSE, comment.char = "#")
spec <- function(signal = "mixed", da = 0.10, balance = "balanced")
  list(n_per_group = 50, da_frac = da, effect = "medium", signal_type = signal, balance = balance,
       conf_phi = 0, exposure = "binary")
S <- list()
for (r in c("R00", "R06", "R08", "R11", "R13", "R15", "R16", "R17", "R19", "R22"))
  S[[length(S) + 1L]] <- list(id = paste0("house:", r), sim = "house", regime = r)
for (r in c("R00", "R06", "R11", "R17", "R19")) for (s in c("msq", "mid"))
  S[[length(S) + 1L]] <- list(id = paste0(s, ":", r), sim = s, regime = r)
imp <- list(B_ref = spec(), B_null = spec(da = 0), B_prev = spec("prevalence"), B_abund = spec("abundance"),
            B_allup = spec(balance = "100_0"))
for (b in names(imp)) S[[length(S) + 1L]] <- list(id = paste0("implant:", b), sim = "implant", spec = imp[[b]])
sims <- strsplit(opt$sims, ",")[[1]]
S <- Filter(function(s) s$sim %in% sims, S)
if (opt$settings != "all") S <- Filter(function(s) s$id %in% strsplit(opt$settings, ",")[[1]], S)
tpls <- strsplit(opt$templates, ",")[[1]]
TPL <- setNames(lapply(tpls, function(id) readRDS(file.path(root, "devdata", paste0(id, ".rds")))), tpls)

make_cell <- function(s, tpl, seed) {
  if (s$sim == "implant") implant(tpl, s$spec, seed)
  else simulate_cell_data(s$sim, tpl, regimes[regimes$regime_id == s$regime, ], seed, cache_dir)
}
cell_seed <- function(...) { k <- paste(..., sep = "|"); as.integer((opt$seed * 1e5 + sum(utf8ToInt(k) * seq_along(utf8ToInt(k)))) %% .Machine$integer.max) }

jobs <- expand.grid(si = seq_along(S), tp = tpls, rep = seq_len(opt$reps), stringsAsFactors = FALSE)
cat(sprintf("devsuite '%s': %d candidates x %d settings x %d templates x %d reps = %d cells, %d cores\n",
            opt$label, length(cand), length(S), length(tpls), opt$reps, nrow(jobs), opt$cores))
cat("candidates:", paste(cand, collapse = ", "), "\n\n")

score <- function(p, truth) {
  ok <- is.finite(p); q <- rep(NA_real_, length(p)); q[ok] <- p.adjust(p[ok], "BH")
  rej <- !is.na(q) & q <= 0.05
  c(tested = sum(ok), n_pos = sum(truth[ok] == 1), rej = sum(rej), tp = sum(rej & truth == 1),
    fp = sum(rej & truth == 0), null_n = sum(ok & truth == 0), null_p05 = sum(ok & truth == 0 & p < 0.05))
}

run_job <- function(k) {
  s <- S[[jobs$si[k]]]; tid <- jobs$tp[k]; rep <- jobs$rep[k]
  sim <- tryCatch(make_cell(s, TPL[[tid]], cell_seed(s$id, tid, rep)), error = function(e) e)
  if (inherits(sim, "error") || !is.null(sim$unsupported)) return(NULL)
  keep <- rowMeans(sim$counts > 0) >= 0.10 & rowSums(sim$counts > 0) >= 3
  ct <- sim$counts[keep, , drop = FALSE]
  truth <- sim$truth$truth_abs[match(rownames(ct), sim$truth$feature)]; truth[is.na(truth)] <- 0L
  meta <- sim$meta[colnames(ct), , drop = FALSE]
  if (is.null(meta$depth)) meta$depth <- colSums(sim$counts)
  out <- lapply(cand, function(cn) {
    t0 <- proc.time()[["elapsed"]]
    r <- tryCatch(.cands[[cn]]$fn(ct, meta, sim$formula, sim$tested_term), error = function(e) e)
    el <- proc.time()[["elapsed"]] - t0
    if (inherits(r, "error")) return(data.frame(setting = s$id, template = tid, rep = rep, candidate = cn,
      tested = 0, n_pos = sum(truth), rej = 0, tp = 0, fp = 0, null_n = 0, null_p05 = 0, secs = el,
      error = substr(conditionMessage(r), 1, 100), stringsAsFactors = FALSE))
    p <- r$p[match(rownames(ct), r$feature)]
    data.frame(setting = s$id, template = tid, rep = rep, candidate = cn, t(score(p, truth)), secs = el,
               error = NA_character_, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}
t_start <- proc.time()[["elapsed"]]
res <- parallel::mclapply(seq_len(nrow(jobs)), run_job, mc.cores = opt$cores, mc.preschedule = FALSE)
D <- do.call(rbind, Filter(Negate(is.null), res))
odir <- file.path(root, "dev", "results", opt$label); dir.create(odir, recursive = TRUE, showWarnings = FALSE)
write.csv(D, file.path(odir, "cells.csv"), row.names = FALSE)

# ---- summary: pooled over templates and reps -------------------------------------------------
sm <- do.call(rbind, lapply(split(D, list(D$setting, D$candidate), drop = TRUE), function(x) data.frame(
  setting = x$setting[1], candidate = x$candidate[1], cells = nrow(x),
  n_pos = mean(x$n_pos), TP = mean(x$tp), FP = mean(x$fp),
  FDR = sum(x$fp) / max(sum(x$rej), 1), FPR = sum(x$null_p05) / max(sum(x$null_n), 1),
  secs = mean(x$secs), errors = sum(!is.na(x$error)), stringsAsFactors = FALSE)))
write.csv(sm, file.path(odir, "summary.csv"), row.names = FALSE)

# ---- scoreboard: the charter's gate, then power ---------------------------------------------
# Gate: false positive rate among null features <= 0.075 in EVERY setting (1.5x nominal: room
# for Monte Carlo noise at a few reps), and pooled FDR <= 0.10 in every setting with signal.
# Power: mean true positives per cell over the signal settings, and relative to the best
# calibrated candidate in each setting (so easy settings do not dominate the average).
sig <- sm$n_pos > 0
best <- vapply(split(sm$TP[sig], sm$setting[sig]), max, numeric(1))
sm$rel <- ifelse(sig, sm$TP / pmax(unname(best[sm$setting]), 1e-9), NA)
sb <- do.call(rbind, lapply(split(sm, sm$candidate), function(x) data.frame(
  candidate = x$candidate[1],
  worst_FPR = max(x$FPR, na.rm = TRUE), worst_FPR_at = x$setting[which.max(x$FPR)],
  worst_FDR = if (any(x$n_pos > 0)) max(x$FDR[x$n_pos > 0]) else NA,
  worst_FDR_at = if (any(x$n_pos > 0)) x$setting[x$n_pos > 0][which.max(x$FDR[x$n_pos > 0])] else NA,
  mean_TP = mean(x$TP[x$n_pos > 0]), rel_power = mean(x$rel, na.rm = TRUE),
  secs_per_cell = mean(x$secs), errors = sum(x$errors), stringsAsFactors = FALSE)))
sb$calibrated <- sb$worst_FPR <= 0.075 & (is.na(sb$worst_FDR) | sb$worst_FDR <= 0.10) & sb$errors == 0
sb <- sb[order(!sb$calibrated, -sb$rel_power), ]
write.csv(sb, file.path(odir, "scoreboard.csv"), row.names = FALSE)

fmt <- function(x, d = 3) formatC(x, format = "f", digits = d)
cat(sprintf("%d cells in %.2f min\n\n=== scoreboard: calibration gate, then power ===\n", nrow(jobs),
            (proc.time()[["elapsed"]] - t_start) / 60))
print(data.frame(candidate = sb$candidate, gate = ifelse(sb$calibrated, "PASS", "fail"),
                 worst_FPR = fmt(sb$worst_FPR), at = sb$worst_FPR_at, worst_FDR = fmt(sb$worst_FDR), at_ = sb$worst_FDR_at,
                 rel_power = fmt(sb$rel_power, 2), mean_TP = fmt(sb$mean_TP, 2), s_cell = fmt(sb$secs_per_cell, 1)),
      row.names = FALSE)
cat("\n=== per setting: TP / FDR / FPR ===\n")
w <- reshape(sm[, c("setting", "candidate", "TP", "FDR", "FPR")], idvar = "setting", timevar = "candidate", direction = "wide")
for (s in sort(unique(sm$setting))) {
  x <- sm[sm$setting == s, ]; x <- x[order(-x$TP), ]
  cat(sprintf("\n%s  (DA per cell %.1f)\n", s, max(x$n_pos)))
  print(data.frame(candidate = x$candidate, TP = fmt(x$TP, 2), FP = fmt(x$FP, 2), FDR = fmt(x$FDR), FPR = fmt(x$FPR),
                   s = fmt(x$secs, 1)), row.names = FALSE)
}
cat("\nwrote ", odir, "\n", sep = "")
