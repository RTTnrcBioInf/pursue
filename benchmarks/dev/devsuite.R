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
  make_option("--suite", default = "core", help = "core (the 20 original settings) | ext (5 stress settings) | full"),
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
# ext: stress settings added 2026-09-24 -- small n (R01), 40% DA (R09), a confounder that also moves
# 10% of null features (R21), a x20 bloom of a common feature (R24: an absolute-vs-relative test),
# repeated measures with subject-level exposure (R23), unbalanced groups with 3x variance in cases
# (R25: every full-benchmark method fails it -- a variance change is not DA under the location truth)
EXT <- c("R01", "R09", "R21", "R23", "R24", "R25")
for (r in EXT) S[[length(S) + 1L]] <- list(id = paste0("house:", r), sim = "house", regime = r, ext = TRUE)
# a hard bloom on top of house R00 (tools/stress.R): the absolute-vs-relative test R24 is too mild for
S[[length(S) + 1L]] <- list(id = "bloom:x4", sim = "house", regime = "R00", ext = TRUE, bloom = 4)
for (r in c("R00", "R06", "R11", "R17", "R19")) for (s in c("msq", "mid"))
  S[[length(S) + 1L]] <- list(id = paste0(s, ":", r), sim = s, regime = r)
imp <- list(B_ref = spec(), B_null = spec(da = 0), B_prev = spec("prevalence"), B_abund = spec("abundance"),
            B_allup = spec(balance = "100_0"))
for (b in names(imp)) S[[length(S) + 1L]] <- list(id = paste0("implant:", b), sim = "implant", spec = imp[[b]])
sims <- strsplit(opt$sims, ",")[[1]]
S <- Filter(function(s) s$sim %in% sims, S)
is_ext <- vapply(S, function(s) isTRUE(s$ext), logical(1))
S <- switch(opt$suite, core = S[!is_ext], ext = S[is_ext], full = S, stop("--suite must be core, ext or full"))
if (opt$settings != "all") S <- Filter(function(s) s$id %in% strsplit(opt$settings, ",")[[1]], S)
tpls <- strsplit(opt$templates, ",")[[1]]
TPL <- setNames(lapply(tpls, function(id) readRDS(file.path(root, "devdata", paste0(id, ".rds")))), tpls)

source(file.path(root, "dev", "tools", "stress.R"))
make_cell <- function(s, tpl, seed) {
  if (s$sim == "implant") return(implant(tpl, s$spec, seed))
  x <- simulate_cell_data(s$sim, tpl, regimes[regimes$regime_id == s$regime, ], seed, cache_dir)
  if (!is.null(s$bloom)) x <- apply_bloom(x, s$bloom, seed)
  x
}
cell_seed <- function(...) { k <- paste(..., sep = "|"); as.integer((opt$seed * 1e5 + sum(utf8ToInt(k) * seq_along(utf8ToInt(k)))) %% .Machine$integer.max) }

jobs <- expand.grid(si = seq_along(S), tp = tpls, rep = seq_len(opt$reps), stringsAsFactors = FALSE)
cat(sprintf("devsuite '%s': %d candidates x %d settings x %d templates x %d reps = %d cells, %d cores\n",
            opt$label, length(cand), length(S), length(tpls), opt$reps, nrow(jobs), opt$cores))
cat("candidates:", paste(cand, collapse = ", "), "\n\n")

# A candidate may return its own q (its FDR procedure is part of the method, as in the benchmark's
# .finish()); otherwise BH on its p. Tested = finite p either way, so a q-level filter cannot shrink
# the set of positives a candidate is scored against.
score <- function(p, truth, qown = NULL) {
  ok <- is.finite(p); q <- rep(NA_real_, length(p))
  if (is.null(qown)) q[ok] <- p.adjust(p[ok], "BH") else q[ok] <- ifelse(is.finite(qown[ok]), qown[ok], 1)
  rej <- !is.na(q) & q <= 0.05
  c(tested = sum(ok), n_pos = sum(truth[ok] == 1), rej = sum(rej), tp = sum(rej & truth == 1),
    fp = sum(rej & truth == 0), null_n = sum(ok & truth == 0), null_p05 = sum(ok & truth == 0 & p < 0.05))
}

# Every finished cell is checkpointed to <odir>/cells/ and reused on a re-run with the same label,
# so an interrupted run (a reclaimed sandbox, a killed job) loses only the cells in flight. A cell
# is keyed by setting, template and replicate; all candidates of the run are stored together, so a
# checkpoint is reused only when it holds every requested candidate.
odir <- file.path(root, "dev", "results", opt$label); cdir <- file.path(odir, "cells")
dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
progress <- file.path(odir, "progress.log")
run_job <- function(k) {
  s <- S[[jobs$si[k]]]; tid <- jobs$tp[k]; rep <- jobs$rep[k]
  ck <- file.path(cdir, sprintf("%s__%s__r%02d.csv", gsub("[:]", "-", s$id), tid, rep))
  if (file.exists(ck)) { prev <- tryCatch(utils::read.csv(ck, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(prev) && all(cand %in% prev$candidate)) return(prev[prev$candidate %in% cand, ]) }
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
    p <- r$p[match(rownames(ct), r$feature)]; qo <- if (!is.null(r$q)) r$q[match(rownames(ct), r$feature)] else NULL
    data.frame(setting = s$id, template = tid, rep = rep, candidate = cn, t(score(p, truth, qo)), secs = el,
               error = NA_character_, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  tmp <- paste0(ck, ".tmp", Sys.getpid()); utils::write.csv(out, tmp, row.names = FALSE); file.rename(tmp, ck)
  cat(sprintf("%s %s %s r%d\n", format(Sys.time(), "%H:%M:%S"), s$id, tid, rep), file = progress, append = TRUE)
  out
}
t_start <- proc.time()[["elapsed"]]
res <- parallel::mclapply(seq_len(nrow(jobs)), run_job, mc.cores = opt$cores, mc.preschedule = FALSE)
D <- do.call(rbind, Filter(Negate(is.null), res))
write.csv(D, file.path(odir, "cells.csv"), row.names = FALSE)

# ---- summary, scoreboard and report: benchmarks/dev/score.R (shared with compare.R) ----------
source(file.path(root, "dev", "score.R"))
sm <- summarise_cells(D); write.csv(sm, file.path(odir, "summary.csv"), row.names = FALSE)
sb <- scoreboard(sm); write.csv(sb, file.path(odir, "scoreboard.csv"), row.names = FALSE)
cat(sprintf("%d cells in %.2f min\n\n", nrow(jobs), (proc.time()[["elapsed"]] - t_start) / 60))
print_report(sm, sb)
cat("\nwrote ", odir, "\n", sep = "")
