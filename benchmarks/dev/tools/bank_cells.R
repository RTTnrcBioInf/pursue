#!/usr/bin/env Rscript
# bank_cells.R -- save dev-suite cells (the simulated data, exactly as devsuite.R builds them: same seeds, same
# prevalence filter) together with per-feature p-values of chosen candidates, so a simulator that only the server
# can run (sd2: SparseDOSSA2 fits take hours) can be diagnosed and iterated on anywhere.
#   Rscript benchmarks/dev/tools/bank_cells.R --settings sd2:R08,msq:R00 --templates hmp_tongue --reps 5 \
#           --candidates eperm_cs1,ref_ldm,ref_zicoseq --out benchmarks/dev/bank/b1 --cores 16
# One RDS per cell: list(id, template, rep, counts, meta, formula (text; as.formula() it), tested_term, truth, p = features x candidates,
# q = the same for candidates returning their own q, eperm = per-option p of eperm_cs1 when it is run).
suppressMessages(library(optparse))
op <- OptionParser(option_list = list(make_option("--settings"), make_option("--templates", default = "hmp_tongue,twinsuk_stool"),
  make_option("--reps", type = "integer", default = 5L), make_option("--candidates", default = "eperm_cs1"),
  make_option("--out"), make_option("--cores", type = "integer", default = 2L), make_option("--seed", type = "integer", default = 9001L)))
opt <- parse_args(op)
root <- normalizePath(Sys.getenv("PURSUE_BENCH_ROOT", unset = "benchmarks"), mustWork = TRUE); Sys.setenv(PURSUE_BENCH_ROOT = root)
if (!nzchar(Sys.getenv("PURSUE_DATA_ROOT"))) Sys.setenv(PURSUE_DATA_ROOT = file.path(root, "data"))
for (f in c("R/engine/templates.R", "R/engine/regimes.R", "R/engine/metrics.R", "R/simulators/sim_house.R",
            "R/simulators/implant.R", "R/simulators/dispatch.R", "R/methods/elementary.R")) source(file.path(root, f))
source(file.path(root, "dev", "tools", "stress.R"))
.cands <- new.env(); register_candidate <- function(name, fn, notes = "") assign(name, list(fn = fn), envir = .cands)
for (f in sort(list.files(file.path(root, "dev", "candidates"), pattern = "\\.R$", full.names = TRUE))) source(f)
cache_dir <- file.path(dirname(root), "cache")
regimes <- read.delim(file.path(root, "regimes.tsv"), stringsAsFactors = FALSE, comment.char = "#")
spec <- function(signal = "mixed", da = 0.10, balance = "balanced") list(n_per_group = 50, da_frac = da, effect = "medium",
  signal_type = signal, balance = balance, conf_phi = 0, exposure = "binary")
imp <- list(B_ref = spec(), B_null = spec(da = 0), B_prev = spec("prevalence"), B_abund = spec("abundance"), B_allup = spec(balance = "100_0"),
            B_conf04 = utils::modifyList(spec(), list(conf_phi = 0.4)), B_conf07 = utils::modifyList(spec(), list(conf_phi = 0.7)),
            B_cont = utils::modifyList(spec("abundance"), list(exposure = "continuous")))
cell_seed <- function(...) { k <- paste(..., sep = "|"); as.integer((opt$seed * 1e5 + sum(utf8ToInt(k) * seq_along(utf8ToInt(k)))) %% .Machine$integer.max) }
make_cell <- function(id, tpl, seed) {                              # devsuite.R's make_cell, including the --extra forms
  nul <- grepl("\\.null$", id); base <- sub("\\.null$", "", id); sm <- sub(":.*", "", base); rg <- sub(".*:", "", base)
  if (sm == "implant") { sp <- imp[[rg]]; if (nul) sp$da_frac <- 0; return(implant(tpl, sp, seed)) }
  if (sm == "bloom") return(apply_bloom(simulate_cell_data("house", tpl, regimes[regimes$regime_id == "R00", ], seed, cache_dir), as.numeric(sub("x", "", rg)), seed))
  r <- regimes[regimes$regime_id == rg, ]; if (nul) { r$da_frac <- 0; r$is_null <- TRUE }
  simulate_cell_data(sm, tpl, r, seed, cache_dir)
}
cand <- strsplit(opt$candidates, ",")[[1]]; tpls <- strsplit(opt$templates, ",")[[1]]
TPL <- setNames(lapply(tpls, function(id) readRDS(file.path(root, "devdata", paste0(id, ".rds")))), tpls)
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
jobs <- expand.grid(id = strsplit(opt$settings, ",")[[1]], tp = tpls, rep = seq_len(opt$reps), stringsAsFactors = FALSE)
jobs <- jobs[!(grepl("^sps:", jobs$id) & jobs$tp != "hmp_tongue"), ]          # sps: hmp_tongue only (no SPsimSeq group in twinsuk)
run <- function(k) {
  id <- jobs$id[k]; tid <- jobs$tp[k]; rp <- jobs$rep[k]
  f <- file.path(opt$out, sprintf("%s__%s__r%02d.rds", gsub("[:]", "-", id), tid, rp)); if (file.exists(f)) return(f)
  sim <- tryCatch(make_cell(id, TPL[[tid]], cell_seed(id, tid, rp)), error = function(e) e)
  if (inherits(sim, "error") || !is.null(sim$unsupported)) return(NA_character_)
  keep <- rowMeans(sim$counts > 0) >= 0.10 & rowSums(sim$counts > 0) >= 3; ct <- sim$counts[keep, , drop = FALSE]
  meta <- sim$meta[colnames(ct), , drop = FALSE]; if (is.null(meta$depth)) meta$depth <- colSums(sim$counts)
  truth <- sim$truth[match(rownames(ct), sim$truth$feature), ]
  P <- matrix(NA_real_, nrow(ct), length(cand), dimnames = list(rownames(ct), cand)); Q <- P; secs <- setNames(rep(NA_real_, length(cand)), cand); ep <- NULL
  for (cn in cand) { t0 <- proc.time()[[3]]
    r <- tryCatch(.cands[[cn]]$fn(ct, meta, sim$formula, sim$tested_term), error = function(e) NULL); secs[cn] <- proc.time()[[3]] - t0
    if (is.null(r)) next
    P[, cn] <- r$p[match(rownames(ct), r$feature)]; if (!is.null(r$q)) Q[, cn] <- r$q[match(rownames(ct), r$feature)]
    if (cn == "eperm_cs1") { v <- .eperm_memo$val; g <- function(x) if (length(x) == nrow(ct)) x else rep(if (length(x) == 1L) x else NA, nrow(ct))
      ep <- data.frame(full = g(v$p_full), lead = g(v$p_ms), det = g(v$p_det), log = g(v$p_log), sqrt = g(v$p_sqrt), choice = g(v$choice),
                       gam = g(v$gam), rho = g(v$rho), perm = isTRUE(v$perm)) } }
  # the formula is stored as text: a formula object carries its environment, which here is the simulator's whole frame
  # (template, fitted parameters) -- bank1's msq twinsuk cells were 54 MB each because of it (2026-10-10)
  saveRDS(list(id = id, template = tid, rep = rp, seed = cell_seed(id, tid, rp), counts = ct, meta = meta, formula = paste(deparse(sim$formula), collapse = " "), tested_term = sim$tested_term,
               truth = truth, p = P, q = Q, secs = secs, eperm = ep), f, compress = "xz")
  cat(format(Sys.time(), "%H:%M:%S"), id, tid, rp, "\n"); f
}
out <- parallel::mclapply(seq_len(nrow(jobs)), run, mc.cores = opt$cores, mc.preschedule = FALSE)
cat(sum(!is.na(unlist(out))), "of", nrow(jobs), "cells banked in", opt$out, "\n")
