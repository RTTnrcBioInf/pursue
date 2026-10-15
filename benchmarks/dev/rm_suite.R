#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# rm_suite.R -- repeated-measures cells on real data (dev/tools/rm.R): subject-level exposure (P_subj) and within-subject
# exposure (P_pair), each with implanted signal and as an exact null (.null). Development template: hmp_tongue (the
# tuning template) with HMP subject labels from devdata/hmp16s_subjects.csv.
#
#   Rscript benchmarks/dev/rm_suite.R --candidates eperm_c12,ref_linda_re --settings P_subj,P_subj.null --reps 20 \
#           --null-reps 40 --cores 8 --label srv30
#
# Same candidates, prevalence filter (>= 10% and >= 3 samples) and dev master seed as devsuite.R. Per cell and candidate:
# calls at the method's own q <= 0.05 (BH where it reports none), true / false calls, and p < 0.05 on null taxa.
# Checkpointed per cell; cells.csv, summary.csv, summary.txt at the end.
# -----------------------------------------------------------------------------
suppressMessages(library(optparse))
op <- OptionParser(option_list = list(
  make_option("--candidates", default = "eperm_c12"),
  make_option("--settings", default = "P_subj,P_pair,P_subj.null,P_pair.null"),
  make_option("--template", default = "hmp_tongue"),
  make_option("--map", default = "devdata/hmp16s_subjects.csv"),
  make_option("--reps", type = "integer", default = 10L),
  make_option("--null-reps", type = "integer", default = 20L),
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
source(file.path(root, "dev", "tools", "rm.R"))

.cands <- new.env()
register_candidate <- function(name, fn, notes = "") assign(name, list(fn = fn, notes = notes), envir = .cands)
for (f in sort(list.files(file.path(root, "dev", "candidates"), pattern = "\\.R$", full.names = TRUE))) source(f)
cand <- strsplit(opt$candidates, ",")[[1]]
miss <- setdiff(cand, ls(.cands)); if (length(miss)) stop("unknown candidates: ", paste(miss, collapse = ", "))
sets <- strsplit(opt$settings, ",")[[1]]
miss <- setdiff(sets, names(RM_SPECS)); if (length(miss)) stop("unknown settings: ", paste(miss, collapse = ", "))

mp <- if (grepl("^/", opt$map)) opt$map else file.path(root, opt$map); if (!file.exists(mp)) stop("no subject map: ", mp, " (run dev/tools/hmp_subjects.R on the server)")
TPL <- rm_template(readRDS(file.path(root, "devdata", paste0(opt$template, ".rds"))), utils::read.csv(mp, stringsAsFactors = FALSE))
cat(sprintf("template %s: %d of the template's samples matched; %d samples after de-duplication, %d subjects, %d with >= 2 visits\n",
            TPL$id, TPL$n_matched, TPL$n_samples, TPL$n_subjects, TPL$n_subj2))
cell_seed <- function(...) { k <- paste(..., sep = "|"); as.integer((opt$seed * 1e5 + sum(utf8ToInt(k) * seq_along(utf8ToInt(k)))) %% .Machine$integer.max) }

jobs <- do.call(rbind, lapply(sets, function(s) data.frame(set = s, rep = seq_len(if (grepl("[.]null$", s)) opt$`null-reps` else opt$reps), stringsAsFactors = FALSE)))
cat(sprintf("rm suite '%s': %d candidates x %d cells, %d cores\n", opt$label, length(cand), nrow(jobs), opt$cores))
odir <- file.path(root, "dev", "results", opt$label); cdir <- file.path(odir, "cells")
dir.create(cdir, recursive = TRUE, showWarnings = FALSE); progress <- file.path(odir, "progress.log")

run_job <- function(k) {
  st <- jobs$set[k]; rep <- jobs$rep[k]
  ck <- file.path(cdir, sprintf("%s__%s__r%02d.csv", st, opt$template, rep))
  if (file.exists(ck)) { prev <- tryCatch(utils::read.csv(ck, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(prev) && all(cand %in% prev$candidate)) return(prev[prev$candidate %in% cand, ]) }
  sim <- rm_cell(TPL, RM_SPECS[[st]], cell_seed("rm", st, opt$template, rep))
  if (!is.null(sim$unsupported)) { message("unsupported: ", sim$unsupported); return(NULL) }
  keep <- rowMeans(sim$counts > 0) >= 0.10 & rowSums(sim$counts > 0) >= 3
  ct <- sim$counts[keep, , drop = FALSE]; pos <- sim$truth$truth_abs[match(rownames(ct), sim$truth$feature)] == 1L
  meta <- sim$meta[colnames(ct), , drop = FALSE]
  if (!is.null(prev <- if (file.exists(ck)) tryCatch(utils::read.csv(ck, stringsAsFactors = FALSE), error = function(e) NULL)))
    prev <- prev[prev$candidate %in% cand & is.na(prev$error), , drop = FALSE]
  todo <- setdiff(cand, prev$candidate)
  out <- lapply(todo, function(cn) {
    t0 <- proc.time()[["elapsed"]]
    r <- tryCatch(.cands[[cn]]$fn(ct, meta, sim$formula, sim$tested_term), error = function(e) e)
    el <- proc.time()[["elapsed"]] - t0
    base <- data.frame(setting = st, template = opt$template, rep = rep, candidate = cn, design = sim$design, m = nrow(ct),
                       n_pos = sum(pos), stringsAsFactors = FALSE)
    if (inherits(r, "error")) return(cbind(base, tested = 0, rej = NA, tp = NA, fp = NA, null_n = NA, null_p05 = NA, secs = el,
                                           error = substr(conditionMessage(r), 1, 120)))
    p <- r$p[match(rownames(ct), r$feature)]; ok <- is.finite(p); q <- rep(NA_real_, length(p))
    if (is.null(r$q)) q[ok] <- stats::p.adjust(p[ok], "BH") else { qo <- r$q[match(rownames(ct), r$feature)]; q[ok] <- ifelse(is.finite(qo[ok]), qo[ok], 1) }
    rj <- !is.na(q) & q <= 0.05
    cbind(base, tested = sum(ok), rej = sum(rj), tp = sum(rj & pos), fp = sum(rj & !pos), null_n = sum(ok & !pos),
          null_p05 = sum(ok & !pos & p < 0.05), secs = el, error = NA_character_)
  })
  out <- do.call(rbind, c(if (!is.null(prev) && nrow(prev)) list(prev), out))
  tmp <- paste0(ck, ".tmp", Sys.getpid()); utils::write.csv(out, tmp, row.names = FALSE); file.rename(tmp, ck)
  cat(sprintf("%s %s r%d\n", format(Sys.time(), "%H:%M:%S"), st, rep), file = progress, append = TRUE)
  out
}
t_start <- proc.time()[["elapsed"]]
res <- parallel::mclapply(seq_len(nrow(jobs)), run_job, mc.cores = opt$cores, mc.preschedule = FALSE)
D <- do.call(rbind, Filter(function(x) is.data.frame(x), res))
write.csv(D, file.path(odir, "cells.csv"), row.names = FALSE)

E <- D[!is.na(D$error), ]; D <- D[is.na(D$error), ]
sm <- do.call(rbind, lapply(split(D, list(D$setting, D$candidate), drop = TRUE), function(x) data.frame(
  setting = x$setting[1], candidate = x$candidate[1], cells = nrow(x), n_pos = mean(x$n_pos), TP = mean(x$tp), FP = mean(x$fp),
  FDR = if (sum(x$rej)) sum(x$fp) / sum(x$rej) else 0, p_any_false = mean(x$fp > 0), fpr05 = sum(x$null_p05) / sum(x$null_n),
  secs = stats::median(x$secs), stringsAsFactors = FALSE)))
sm <- sm[order(sm$setting, -sm$TP, sm$fpr05), ]
write.csv(sm, file.path(odir, "summary.csv"), row.names = FALSE)
options(width = 200); sm$secs <- round(sm$secs, 1)
sink(file.path(odir, "summary.txt"), split = TRUE)
cat(sprintf("%s: %d samples, %d subjects (%d with >= 2 visits). %d cells in %.1f min.\n", TPL$id, TPL$n_samples, TPL$n_subjects, TPL$n_subj2,
            nrow(jobs), (proc.time()[["elapsed"]] - t_start) / 60))
cat("TP / FP per cell at the method's own q <= 0.05; FDR pooled; p_any_false = share of cells with any false call;\n",
    "fpr05 = share of null taxa at p < 0.05 (nominal 0.05). On .null settings every call is false.\n\n", sep = "")
for (s in unique(sm$setting)) { cat("==", s, "\n"); print(format(sm[sm$setting == s, -1], digits = 3), row.names = FALSE); cat("\n") }
if (nrow(E)) { cat("errors (cells):\n"); print(table(E$candidate, substr(E$error, 1, 70))) }
sink()
