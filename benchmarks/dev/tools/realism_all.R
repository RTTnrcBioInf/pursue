#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# realism_all.R -- the protocol's realism gate (section 4.4; R/engine/realism.R, unchanged) for every simulator on the
# evaluation templates, as the protocol specifies, and on the two tuning templates. The gate was planned as a SLURM array
# (hpc/slurm/realism.sbatch); hpc/run_local.sh cannot run its task lines, so it was never computed on calc01.
#
# Beside the gate's three numbers (random-forest AUROC real vs simulated, median per-taxon KS, Frobenius distance of the
# Spearman matrices) it records why a simulator is separable, on the very table the gate scored (same seed):
#   cor50   mean |Spearman| among the 50 most prevalent taxa (co-occurrence / community structure)
#   sd_rich sd of per-sample richness (between-sample heterogeneity)
#   rich    median richness;  top1  median top-taxon share (dominance);  zero  fraction of zero cells
# Real-data values use all of the template's samples.
#
# Tuning templates are read from benchmarks/devdata/ (as the dev suite reads them, so the cached sparseDOSSA2 fits match);
# evaluation templates through load_template(), as the full run did. sparseDOSSA2 and MIDASim reuse cache/.
# Resumable: one csv per (simulator, template) in results/realism/, finished files kept; unsupported pairs and failures
# (<file>.err) are retried on the next run. Summary: benchmarks/dev/results/realism/{summary.txt, realism.csv}.
#   Rscript benchmarks/dev/tools/realism_all.R --cores 6
# -----------------------------------------------------------------------------
suppressMessages(library(optparse))
opt <- parse_args(OptionParser(option_list = list(
  make_option("--sims", type = "character", default = "sd2,sps,mid,msq,house"),
  make_option("--templates", type = "character", default = "risk_stool,hmp_stool,hmp_gingiva,hmp_skin_ear,hmp_vagina,hmp_tongue,twinsuk_stool"),
  make_option("--cores", type = "integer", default = 6L))))
root <- normalizePath(Sys.getenv("PURSUE_BENCH_ROOT", unset = "benchmarks"), mustWork = TRUE); Sys.setenv(PURSUE_BENCH_ROOT = root)
for (f in c("R/engine/templates.R", "R/engine/regimes.R", "R/engine/metrics.R", "R/simulators/sim_house.R",
            "R/simulators/implant.R", "R/simulators/dispatch.R", "R/engine/realism.R")) source(file.path(root, f))
cache <- file.path(dirname(root), "cache"); out <- file.path(dirname(root), "results", "realism")
rdir <- file.path(root, "dev", "results", "realism"); for (d in c(cache, out, rdir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# keep the simulated table the gate scores (realism_gate calls simulate_cell_data by name). The gate sets m to the
# template's feature count; sparseDOSSA2 always keeps the template's own features and so declares the m sweep unsupported,
# which made it refuse every realism task (first run, 2026-10-15). m is not a regime factor here: give sd2 the reference m.
.last <- new.env(); .sim0 <- simulate_cell_data
simulate_cell_data <- function(simulator, template, regime, seed, cache_dir = NULL) {
  if (identical(simulator, "sd2")) regime$m <- reference_regime()$m
  s <- .sim0(simulator, template, regime, seed, cache_dir); .last$s <- s; s }

treg <- read_template_registry()
tpl_load <- function(id) {
  p <- file.path(root, "devdata", paste0(id, ".rds"))
  if (identical(treg$pool[treg$id == id], "tuning") && file.exists(p)) { x <- readRDS(p); x$id <- id; x } else load_template(id) }
desc <- function(ct) {
  ct <- ct[, colSums(ct) > 0, drop = FALSE]; R <- sweep(ct, 2L, colSums(ct), "/")
  o <- order(-rowMeans(ct > 0))[seq_len(min(50L, nrow(ct)))]
  cm <- suppressWarnings(stats::cor(t(log1p(R[o, , drop = FALSE] * 1e4)), method = "spearman"))
  c(cor50 = mean(abs(cm[upper.tri(cm)]), na.rm = TRUE), sd_rich = stats::sd(colSums(ct > 0)), rich = stats::median(colSums(ct > 0)),
    top1 = stats::median(apply(R, 2L, max)), zero = mean(ct == 0)) }

sims <- strsplit(opt$sims, ",")[[1]]; tpls <- strsplit(opt$templates, ",")[[1]]
tasks <- expand.grid(tpl = tpls, sim = sims, stringsAsFactors = FALSE)[, c("sim", "tpl")]   # heavy simulators first
one <- function(k) {
  s <- tasks$sim[k]; t <- tasks$tpl[k]; f <- file.path(out, sprintf("%s__%s.csv", s, t)); fe <- paste0(f, ".err")
  if (file.exists(f) && is.finite(utils::read.csv(f)$auroc[1])) return(sprintf("kept  %s %s", s, t))   # unsupported rows are retried
  t0 <- Sys.time()
  res <- tryCatch({
    if (!simulator_available(s)) stop("simulator package not installed")
    tp <- tpl_load(t); .last$s <- NULL
    g <- realism_gate(tp, s, n_sim = 200L, seed = 1L, cache_dir = cache)
    if (!is.null(.last$s$counts)) { dr <- desc(tp$counts); ds <- desc(.last$s$counts)
      g <- cbind(g, as.data.frame(as.list(stats::setNames(c(dr, ds), c(paste0("real_", names(dr)), paste0("sim_", names(ds))))))) }
    g$pool <- treg$pool[treg$id == t]; g$seconds <- round(as.numeric(difftime(Sys.time(), t0, units = "secs"))); g
  }, error = function(e) e)
  if (inherits(res, "error")) { writeLines(conditionMessage(res), fe); return(sprintf("FAIL  %s %s: %s", s, t, conditionMessage(res))) }
  utils::write.csv(res, f, row.names = FALSE); if (file.exists(fe)) file.remove(fe)
  if (!is.finite(res$auroc[1])) return(sprintf("unsup %s %s (%s)", s, t, res$note[1]))
  sprintf("ok    %s %s  AUROC %.3f  %ds", s, t, res$auroc[1], res$seconds[1]) }
msg <- parallel::mclapply(seq_len(nrow(tasks)), function(k) { m <- one(k); message(m); m }, mc.cores = opt$cores, mc.preschedule = FALSE)
msg <- vapply(msg, function(m) if (is.character(m)) m else paste("worker died:", conditionMessage(attr(m, "condition") %||% simpleError("?"))), "")

# summary
files <- file.path(out, sprintf("%s__%s.csv", tasks$sim, tasks$tpl)); have <- file.exists(files)
rows <- lapply(files[have], utils::read.csv, stringsAsFactors = FALSE)
cols <- unique(unlist(lapply(rows, names)))
d <- do.call(rbind, lapply(rows, function(r) { for (c in setdiff(cols, names(r))) r[[c]] <- NA; r[, cols] }))
sink(file.path(rdir, "summary.txt"))
cat(format(Sys.time()), "realism gate (protocol 4.4), 200 real vs 200 simulated samples, seed 1\n")
cat("AUROC > 0.9 = low realism (down-weighted 0.5 in the summary ranking). real_* use all template samples.\n\n")
for (m in msg) cat(m, "\n"); cat("\n")
if (!is.null(d) && nrow(d)) {
  num <- vapply(d, is.numeric, logical(1)); d[num] <- lapply(d[num], round, 3)
  d <- d[order(match(d$simulator, c("house", "msq", "mid", "sd2", "sps")), d$pool, d$template), ]
  utils::write.csv(d, file.path(rdir, "realism.csv"), row.names = FALSE)
  show <- intersect(c("simulator", "template", "pool", "auroc", "ks_median", "frob", "real_cor50", "sim_cor50", "real_sd_rich", "sim_sd_rich",
                      "real_rich", "sim_rich", "real_top1", "sim_top1", "real_zero", "sim_zero"), names(d))
  print(d[, show], row.names = FALSE)
  e <- d[d$pool == "evaluation" & is.finite(d$auroc), ]
  if (nrow(e)) { cat("\nmean over evaluation templates\n")
    a <- stats::aggregate(e[, intersect(c("auroc", "ks_median", "frob", "sim_cor50", "real_cor50", "sim_sd_rich", "real_sd_rich"), names(e))],
                          list(simulator = e$simulator), mean, na.rm = TRUE)
    a[-1] <- lapply(a[-1], round, 3); a$flagged <- tapply(e$auroc > 0.9, e$simulator, sum)[a$simulator]; print(a, row.names = FALSE) }
}
sink()
cat(readLines(file.path(rdir, "summary.txt")), sep = "\n")
