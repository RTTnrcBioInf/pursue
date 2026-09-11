#!/usr/bin/env Rscript
# Smoke test for the HPC environment: loads every template, runs every simulator on a tiny
# regime, runs every method on a tiny dataset, and records what works. Run this before
# submitting any array job; the tables it writes (hpc/smoke_*.csv) are the record of which
# wrappers are verified on this cluster, and hpc/smoke_errors.log holds the full messages.
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), "..", "benchmarks"), mustWork = TRUE)
hpc <- normalizePath(file.path(root, "..", "hpc"))
Sys.setenv(PURSUE_BENCH_ROOT = root); setwd(root)
for (f in c("R/engine/templates.R", "R/engine/regimes.R", "R/engine/metrics.R", "R/engine/io.R", "R/engine/realism.R",
            "R/simulators/sim_house.R", "R/simulators/implant.R", "R/simulators/nullperm.R", "R/simulators/dispatch.R",
            "R/methods/elementary.R", "R/methods/external.R", "R/methods/registry.R")) source(f)

errlog <- file.path(hpc, "smoke_errors.log")
cat("smoke test", format(Sys.time()), "\nR", R.version.string, "\n\n", file = errlog)
trunc1 <- function(x, n = 120L) { x <- gsub("[\r\n]+", " ", paste(x, collapse = " ")); if (nchar(x) > n) paste0(substr(x, 1L, n), " [...]") else x }
note <- function(what, e) { cat("=== ", what, "\n", conditionMessage(e), "\n\n", sep = "", file = errlog, append = TRUE)
                           paste("error:", trunc1(conditionMessage(e))) }

cat("== templates ==\n")
reg <- read_template_registry(); rows <- list()
for (i in seq_len(nrow(reg))) {
  id <- reg$id[i]
  r <- tryCatch({ t <- load_template(id)
      list(status = "ok", n_features = nrow(t$counts), n_samples = ncol(t$counts), depth_median = as.integer(median(colSums(t$counts)))) },
    error = function(e) list(status = note(paste("template", id), e), n_features = NA, n_samples = NA, depth_median = NA))
  cat(sprintf("  %-22s %-6s %6s features x %5s samples, depth median %s\n", id, r$status, r$n_features, r$n_samples, r$depth_median))
  rows[[i]] <- data.frame(template = id, pool = reg$pool[i], status = r$status, n_features = r$n_features,
                          n_samples = r$n_samples, depth_median = r$depth_median, stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, rows), file.path(hpc, "smoke_templates.csv"), row.names = FALSE)

tpl <- load_template("hmp_tongue", max_samples = 120)
regime <- reference_regime(); regime$n_per_group <- 15L; regime$m <- 80L
cat("\n== simulators (reference regime, 15/group, 80 features) ==\n")
sreg <- simulator_registry(); rows <- list()
for (i in seq_len(nrow(sreg))) {
  s <- sreg$id[i]
  st <- tryCatch({ if (!simulator_available(s)) paste0("not_installed (", sreg$package[i], ")") else {
        o <- simulate_cell_data(s, tpl, regime, 1L, file.path(root, "cache"))
        if (!is.null(o$unsupported)) paste("unsupported:", paste(o$unsupported, collapse = ",")) else
          sprintf("ok: %d x %d, %d DA", nrow(o$counts), ncol(o$counts), sum(o$truth$truth_abs)) } },
    error = function(e) note(paste("simulator", s), e))
  cat(sprintf("  %-8s %s\n", s, st))
  rows[[i]] <- data.frame(simulator = s, package = sreg$package[i], status = st,
                          version = tryCatch(as.character(utils::packageVersion(sreg$package[i])), error = function(e) NA),
                          stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, rows), file.path(hpc, "smoke_simulators.csv"), row.names = FALSE)

sim <- simulate_house(tpl, regime, 2L); keep <- rowMeans(sim$counts > 0) >= 0.1; ct <- sim$counts[keep, ]
cat("\n== methods (house data, ", nrow(ct), " features x ", ncol(ct), " samples) ==\n", sep = "")
mreg <- method_registry(); rows <- list()
for (i in seq_len(nrow(mreg))) {
  m <- mreg$id[i]
  run <- tryCatch(run_method(m, ct, sim$meta, sim$formula, sim$tested_term, timeout_s = 900),
                  error = function(e) list(result = data.frame(status = note(paste("method", m), e)), runtime_s = NA, version = NA))
  st <- trunc1(paste(unique(run$result$status), collapse = "/"))
  if (identical(st, "not_installed")) st <- paste0("not_installed (", mreg$package[i], ")")
  n_p <- if ("p" %in% names(run$result)) sum(is.finite(run$result$p)) else NA
  cat(sprintf("  %-18s %7.1fs  %s\n", m, if (is.na(run$runtime_s)) 0 else run$runtime_s, st))
  rows[[i]] <- data.frame(method = m, package = mreg$package[i], status = st, n_p_finite = n_p,
                          runtime_s = run$runtime_s, version = run$version, stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, rows), file.path(hpc, "smoke_methods.csv"), row.names = FALSE)
cat("\nwrote hpc/smoke_{templates,simulators,methods}.csv and smoke_errors.log\n")
