#!/usr/bin/env Rscript
# Warm the per-template simulator caches ONCE, before submitting any array job.
#
# MIDASim needs a setup object per (template, m) and sparseDOSSA2 a fit per template. Both are
# slow -- the sparseDOSSA2 fit is minutes to hours on a full template. Without a warm cache
# every array task discovers the miss at the same moment and refits in parallel, which is the
# single most expensive mistake available here. Run this on one node first; the array jobs
# then only read.
#
#   Rscript hpc/prep_caches.R                      # evaluation pool, m = 500 and 1000
#   Rscript hpc/prep_caches.R --pool all --simulators mid
#   Rscript hpc/prep_caches.R --templates hmp_stool,hmp_vagina --timeout 7200
suppressMessages(library(optparse))
opt <- parse_args(OptionParser(option_list = list(
  make_option("--pool", type = "character", default = "evaluation"),
  make_option("--templates", type = "character", default = "all"),
  make_option("--simulators", type = "character", default = "mid,sd2"),
  make_option("--m", type = "character", default = "500,1000", help = "feature counts MIDASim needs a setup for"),
  make_option("--cache", type = "character", default = "cache"),
  make_option("--timeout", type = "integer", default = 21600L, help = "seconds per (simulator, template)"),
  make_option("--bench-root", type = "character", default = NULL))))
root <- if (!is.null(opt$`bench-root`)) opt$`bench-root` else
  normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), "..", "benchmarks"), mustWork = TRUE)
Sys.setenv(PURSUE_BENCH_ROOT = root)
for (f in c("R/engine/templates.R", "R/engine/regimes.R", "R/simulators/sim_house.R", "R/simulators/implant.R",
            "R/simulators/nullperm.R", "R/simulators/dispatch.R")) source(file.path(root, f))

cache <- if (grepl("^(/|[A-Za-z]:)", opt$cache)) opt$cache else file.path(root, "..", opt$cache)
dir.create(cache, recursive = TRUE, showWarnings = FALSE); cache <- normalizePath(cache)
reg <- read_template_registry()
ids <- if (opt$templates != "all") strsplit(opt$templates, ",")[[1]] else
  reg$id[reg$pool %in% (if (opt$pool == "all") unique(reg$pool) else strsplit(opt$pool, ",")[[1]])]
sims <- strsplit(opt$simulators, ",")[[1]]
ms <- as.integer(strsplit(opt$m, ",")[[1]])
cat("R:", R.version.string, "| library:", .libPaths()[1], "\n")
if (getRversion() < "4.5")
  cat("!! R < 4.5 -- did you forget `conda activate pursue-bench`?\n")
cat("cache:", cache, "\ntemplates:", paste(ids, collapse = ", "), "\nsimulators:", paste(sims, collapse = ", "), "\n\n")

rows <- list()
for (s in sims) {
  if (!simulator_available(s)) { cat(sprintf("%-5s skipped: package not installed\n", s)); next }
  for (id in ids) {
    tpl <- tryCatch(load_template(id), error = function(e) NULL)
    if (is.null(tpl)) { cat(sprintf("%-5s %-22s template failed to load\n", s, id)); next }
    for (m in if (s == "mid") ms else NA_integer_) {
      rg <- reference_regime(); rg$n_per_group <- 5L; if (!is.na(m)) rg$m <- m
      lab <- if (is.na(m)) id else paste0(id, " m=", m)
      t0 <- Sys.time()
      st <- tryCatch({ setTimeLimit(elapsed = opt$timeout, transient = TRUE); on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
                       invisible(simulate_cell_data(s, tpl, rg, 1L, cache)); "ok" },
                     error = function(e) paste("FAILED:", conditionMessage(e)))
      setTimeLimit(elapsed = Inf)
      el <- as.numeric(Sys.time() - t0, units = "secs")
      cat(sprintf("%-5s %-30s %8.1fs  %s\n", s, lab, el, st))
      rows[[length(rows) + 1L]] <- data.frame(simulator = s, template = id, m = m, seconds = el, status = st, stringsAsFactors = FALSE)
    }
  }
}
if (length(rows)) {
  d <- do.call(rbind, rows)
  write.csv(d, file.path(dirname(root), "hpc", "prep_caches.csv"), row.names = FALSE)
  cat("\n", sum(d$status == "ok"), "/", nrow(d), " warmed. Cache files:\n", sep = "")
  print(data.frame(file = basename(list.files(cache, pattern = "\\.rds$")),
                   MB = round(file.size(list.files(cache, pattern = "\\.rds$", full.names = TRUE)) / 1e6, 1)), row.names = FALSE)
  cat("wrote hpc/prep_caches.csv\n")
} else cat("nothing to do\n")
