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
  make_option("--template", default = "hmp_tongue", help = "template(s) for the P_* settings, comma-separated: a dev template (devdata/<id>.rds) or a registry id (load_template)"),
  make_option("--pick-random", default = "", help = "templates whose subjects have many samples: two per subject at random (else the two earliest visits)"),
  make_option("--out", default = "", help = "output directory (default benchmarks/dev/results/<label>)"),
  make_option("--house-templates", default = "hmp_tongue,twinsuk_stool", help = "templates for the house settings (HW*)"),
  make_option("--map", default = "devdata/hmp16s_subjects.csv"),
  make_option("--reps", type = "integer", default = 10L),
  make_option("--null-reps", type = "integer", default = 20L),
  make_option("--cores", type = "integer", default = 2L),
  make_option("--label", default = format(Sys.time(), "%Y%m%d-%H%M")),
  make_option("--seed", type = "integer", default = 9001L),
  make_option("--dry", action = "store_true", default = FALSE, help = "load the templates, report subjects and visits, and stop (no cells)")))
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

is_house <- function(s) identical(RM_SPECS[[s]]$gen, "house")
ptp <- strsplit(opt$template, ",")[[1]]; prand <- strsplit(opt$`pick-random`, ",")[[1]]
TPL <- list(); if (!all(vapply(sets, is_house, logical(1)))) for (id in ptp) {
  dv <- file.path(root, "devdata", paste0(id, ".rds"))
  x <- if (file.exists(dv)) readRDS(dv) else load_template(id)          # conf2 (protocol 13): registry templates, held-out data
  if (!is.null(x$meta$X.SampleID)) {                                      # HMP (MicrobeDS): subjects and visits from the sample map
    mp <- if (grepl("^/", opt$map)) opt$map else file.path(root, opt$map); if (!file.exists(mp)) stop("no subject map: ", mp, " (run dev/tools/hmp_subjects.R on the server)")
    t1 <- rm_template(x, utils::read.csv(mp, stringsAsFactors = FALSE))
  } else if (all(c("subject", "visit") %in% names(x$meta))) {
    t1 <- list(counts = x$counts, meta = x$meta[, c("subject", "visit")], id = id, n_matched = ncol(x$counts), n_samples = ncol(x$counts),
               n_subjects = length(unique(x$meta$subject)), n_subj2 = sum(table(x$meta$subject) >= 2L))
  } else stop("template ", id, ": no subject labels (neither X.SampleID for the HMP map nor subject / visit columns)")
  t1$id <- id; t1$pick <- if (id %in% prand) "random" else "earliest"; TPL[[id]] <- t1
  cat(sprintf("template %s: %d samples matched; %d after de-duplication, %d subjects, %d with >= 2 samples; pairs: %s\n",
              id, t1$n_matched, t1$n_samples, t1$n_subjects, t1$n_subj2, t1$pick)) }
if (opt$dry) {                                                            # conf2 pre-launch check: design only, no counts analysed
  need <- max(c(0L, vapply(sets[!vapply(sets, is_house, logical(1))], function(s) RM_SPECS[[s]]$n_subj, 1L)))
  bad <- character(0)
  for (t1 in TPL) { k <- as.integer(table(t1$meta$subject))
    cat(sprintf("  %s: %d taxa; samples per subject min %d, median %g, max %d; %s\n", t1$id, nrow(t1$counts), min(k), stats::median(k), max(k),
                if (t1$n_subj2 >= need) sprintf("OK (%d >= %d subjects with 2 samples)", t1$n_subj2, need) else { bad <- c(bad, t1$id); sprintf("TOO FEW (%d < %d)", t1$n_subj2, need) })) }
  cat(if (length(bad)) paste("rm dry check FAILED:", paste(bad, collapse = ", ")) else "rm dry check: all templates usable", "\n")
  quit(save = "no", status = as.integer(length(bad) > 0))
}
htp <- strsplit(opt$`house-templates`, ",")[[1]]
HT <- setNames(lapply(htp, function(id) readRDS(file.path(root, "devdata", paste0(id, ".rds")))), htp)
cell_seed <- function(...) { k <- paste(..., sep = "|"); as.integer((opt$seed * 1e5 + sum(utf8ToInt(k) * seq_along(utf8ToInt(k)))) %% .Machine$integer.max) }

jobs <- do.call(rbind, lapply(sets, function(s) expand.grid(set = s, tp = if (is_house(s)) htp else ptp,
  rep = seq_len(if (grepl("[.]null$", s)) opt$`null-reps` else opt$reps), stringsAsFactors = FALSE)))
cat(sprintf("rm suite '%s': %d candidates x %d cells, %d cores\n", opt$label, length(cand), nrow(jobs), opt$cores))
odir <- if (nzchar(opt$out)) opt$out else file.path(root, "dev", "results", opt$label); cdir <- file.path(odir, "cells")
dir.create(cdir, recursive = TRUE, showWarnings = FALSE); progress <- file.path(odir, "progress.log")

run_job <- function(k) {
  st <- jobs$set[k]; rep <- jobs$rep[k]; tid <- jobs$tp[k]
  ck <- file.path(cdir, sprintf("%s__%s__r%02d.csv", st, tid, rep))
  if (file.exists(ck)) { prev <- tryCatch(utils::read.csv(ck, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(prev) && all(cand %in% prev$candidate)) return(prev[prev$candidate %in% cand, ]) }
  sim <- if (is_house(st)) simulate_house_within(HT[[tid]], RM_SPECS[[st]], cell_seed("rm", st, tid, rep)) else
    rm_cell(TPL[[tid]], RM_SPECS[[st]], cell_seed("rm", st, tid, rep))
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
    base <- data.frame(setting = st, template = tid, rep = rep, candidate = cn, design = sim$design, m = nrow(ct),
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
summ <- function(D, by) do.call(rbind, lapply(split(D, D[by], drop = TRUE), function(x) cbind(x[1, by, drop = FALSE], data.frame(
  cells = nrow(x), n_pos = mean(x$n_pos), TP = mean(x$tp), FP = mean(x$fp),
  FDR = if (sum(x$rej)) sum(x$fp) / sum(x$rej) else 0, p_any_false = mean(x$fp > 0), fpr05 = sum(x$null_p05) / sum(x$null_n),
  secs = stats::median(x$secs), stringsAsFactors = FALSE))))
sm <- summ(D, c("setting", "candidate")); sm <- sm[order(sm$setting, -sm$TP, sm$fpr05), ]
write.csv(sm, file.path(odir, "summary.csv"), row.names = FALSE)
multi <- length(unique(D$template)) > 1L
if (multi) { st <- summ(D, c("setting", "template", "candidate")); st <- st[order(st$setting, st$template, -st$TP, st$fpr05), ]
  write.csv(st, file.path(odir, "summary_by_template.csv"), row.names = FALSE) }
options(width = 200); sm$secs <- round(sm$secs, 1)
sink(file.path(odir, "summary.txt"), split = TRUE)
for (t in TPL) cat(sprintf("%s: %d samples, %d subjects (%d with >= 2 samples); pairs: %s\n", t$id, t$n_samples, t$n_subjects, t$n_subj2, t$pick))
cat(sprintf("%d cells in %.1f min.\n", nrow(jobs), (proc.time()[["elapsed"]] - t_start) / 60))
cat("TP / FP per cell at the method's own q <= 0.05; FDR pooled; p_any_false = share of cells with any false call;\n",
    "fpr05 = share of null taxa at p < 0.05 (nominal 0.05). On .null settings every call is false.\n\n", sep = "")
for (s in unique(sm$setting)) { cat("==", s, if (multi) "(all templates)", "\n"); print(format(sm[sm$setting == s, -1], digits = 3), row.names = FALSE); cat("\n") }
if (multi) for (s in unique(st$setting)) for (t in unique(st$template)) { x <- st[st$setting == s & st$template == t, ]
  if (nrow(x)) { cat("==", s, "|", t, "\n"); print(format(x[, -(1:2)], digits = 3), row.names = FALSE); cat("\n") } }
if (nrow(E)) { cat("errors (cells):\n"); print(table(E$candidate, substr(E$error, 1, 70))) }
sink()
