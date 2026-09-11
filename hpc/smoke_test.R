#!/usr/bin/env Rscript
# Smoke test for the HPC environment: loads every template, runs every simulator on a
# tiny regime, runs every method on a tiny dataset, and prints what works. Run this
# before submitting any array job; the tables it writes (hpc/smoke_*.csv) are the record
# of which wrappers are verified on this cluster.
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), "..", "benchmarks"), mustWork = TRUE)
Sys.setenv(PURSUE_BENCH_ROOT = root); setwd(root)
for (f in c("R/engine/templates.R", "R/engine/regimes.R", "R/engine/metrics.R", "R/engine/io.R", "R/engine/realism.R",
            "R/simulators/sim_house.R", "R/simulators/implant.R", "R/simulators/nullperm.R", "R/simulators/dispatch.R",
            "R/methods/elementary.R", "R/methods/external.R", "R/methods/registry.R")) source(f)
cat("== templates ==\n"); reg <- read_template_registry(); tt <- list()
for (i in seq_len(nrow(reg))) { r <- tryCatch({ t <- load_template(reg$id[i]); sprintf("%d features x %d samples, depth median %d", nrow(t$counts), ncol(t$counts), as.integer(median(colSums(t$counts)))) },
  error = function(e) paste("FAIL:", conditionMessage(e))); cat(sprintf("  %-22s %s\n", reg$id[i], r)); tt[[reg$id[i]]] <- r }
write.csv(data.frame(template = names(tt), status = unlist(tt)), file.path(root, "..", "hpc", "smoke_templates.csv"), row.names = FALSE)

tpl <- load_template("hmp_tongue", max_samples = 120)
regime <- reference_regime(); regime$n_per_group <- 15L; regime$m <- 80L
cat("\n== simulators (reference regime, 15/group, 80 features) ==\n"); ss <- list()
for (s in simulator_registry()$id) {
  r <- tryCatch({ if (!simulator_available(s)) "not installed" else { o <- simulate_cell_data(s, tpl, regime, 1L, file.path(root, "cache"))
      if (!is.null(o$unsupported)) paste("unsupported:", paste(o$unsupported, collapse = ",")) else sprintf("ok: %d x %d, %d DA", nrow(o$counts), ncol(o$counts), sum(o$truth$truth_abs)) } },
    error = function(e) paste("FAIL:", conditionMessage(e))); cat(sprintf("  %-8s %s\n", s, r)); ss[[s]] <- r }
write.csv(data.frame(simulator = names(ss), status = unlist(ss)), file.path(root, "..", "hpc", "smoke_simulators.csv"), row.names = FALSE)

sim <- simulate_house(tpl, regime, 2L); keep <- rowMeans(sim$counts > 0) >= 0.1; ct <- sim$counts[keep, ]
cat("\n== methods (house data, ", nrow(ct), " features x ", ncol(ct), " samples) ==\n", sep = ""); ms <- list()
for (m in method_registry()$id) {
  run <- run_method(m, ct, sim$meta, sim$formula, sim$tested_term, timeout_s = 900)
  st <- paste(unique(run$result$status), collapse = "/"); cat(sprintf("  %-18s %7.1fs  %s\n", m, run$runtime_s, st)); ms[[m]] <- c(st, run$runtime_s, run$version) }
write.csv(data.frame(method = names(ms), status = sapply(ms, `[`, 1), runtime_s = sapply(ms, `[`, 2), version = sapply(ms, `[`, 3)),
          file.path(root, "..", "hpc", "smoke_methods.csv"), row.names = FALSE)
cat("\nwrote hpc/smoke_templates.csv, smoke_simulators.csv, smoke_methods.csv\n")
