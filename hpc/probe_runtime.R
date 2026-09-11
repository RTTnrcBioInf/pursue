#!/usr/bin/env Rscript
# Measure how each method's runtime scales with feature count and sample size, so the compute
# budget for the benchmark is a measurement rather than an extrapolation.
#
# Motivation: the 2026-09-11 smoke test put fastEmu at 70.9 s on 60 features x 30 samples while
# every other method sat between 0.02 s and 9.4 s. Whether that is tolerable at the reference
# regime (500 features, 100 samples) depends entirely on the exponent -- linear in features is
# ~33 min per cell, quadratic is ~4.5 h -- and the full grid is ~17 300 cells. One is a large
# but payable bill; the other is two CPU-years for a single method.
#
#   Rscript hpc/probe_runtime.R                       # ~40 min, all methods
#   Rscript hpc/probe_runtime.R --methods fastemu,ancombc2,locom --timeout 1800
#
# Writes hpc/probe_runtime.csv and prints the fitted log-log slope in m per method, plus the
# projected seconds at the reference regime.
suppressMessages(library(optparse))
opt <- parse_args(OptionParser(option_list = list(
  make_option("--m", type = "character", default = "30,60,90,120", help = "feature counts"),
  make_option("--n", type = "character", default = "30,60", help = "total samples"),
  make_option("--methods", type = "character", default = "all"),
  make_option("--template", type = "character", default = "hmp_tongue"),
  make_option("--timeout", type = "integer", default = 900L, help = "seconds per method per cell"),
  make_option("--ref-m", type = "integer", default = 500L), make_option("--ref-n", type = "integer", default = 100L),
  make_option("--bench-root", type = "character", default = NULL))))
root <- if (!is.null(opt$`bench-root`)) opt$`bench-root` else
  normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), "..", "benchmarks"), mustWork = TRUE)
Sys.setenv(PURSUE_BENCH_ROOT = root); hpc <- normalizePath(file.path(root, "..", "hpc"))
for (f in c("R/engine/templates.R", "R/engine/regimes.R", "R/simulators/sim_house.R", "R/simulators/implant.R",
            "R/simulators/nullperm.R", "R/simulators/dispatch.R", "R/methods/elementary.R",
            "R/methods/external.R", "R/methods/registry.R")) source(file.path(root, f))

cat("R:", R.version.string, "| library:", .libPaths()[1], "\n")
ms <- as.integer(strsplit(opt$m, ",")[[1]]); ns <- as.integer(strsplit(opt$n, ",")[[1]])
mids <- if (opt$methods == "all") method_registry()$id else strsplit(opt$methods, ",")[[1]]
tpl <- load_template(opt$template, max_samples = max(ns) * 2L)
cat("grid: m =", paste(ms, collapse = ","), " n =", paste(ns, collapse = ","),
    " methods =", length(mids), " timeout =", opt$timeout, "s\n\n")

rows <- list()
for (n in ns) for (m in ms) {
  rg <- reference_regime(); rg$n_per_group <- as.integer(n / 2); rg$m <- m
  sim <- simulate_house(tpl, rg, 1L)
  keep <- rowMeans(sim$counts > 0) >= 0.1; ct <- sim$counts[keep, , drop = FALSE]
  cat(sprintf("-- m=%d n=%d (%d features kept)\n", m, n, nrow(ct)))
  for (mid in mids) {
    if (!method_available(mid)) next
    run <- run_method(mid, ct, sim$meta, sim$formula, sim$tested_term, timeout_s = opt$timeout)
    ok <- any(is.finite(run$result$p))
    cat(sprintf("   %-18s %8.2fs %s\n", mid, run$runtime_s, if (ok) "" else "(no p-values)"))
    rows[[length(rows) + 1L]] <- data.frame(method = mid, m = nrow(ct), n = n, seconds = run$runtime_s,
                                            usable = ok, stringsAsFactors = FALSE)
    write.csv(do.call(rbind, rows), file.path(hpc, "probe_runtime.csv"), row.names = FALSE)
  }
}

d <- do.call(rbind, rows)
cat("\n== scaling ==\nslope_m is the log-log exponent in feature count: 1 = linear, 2 = quadratic.\n")
cat("projected = seconds for one cell at the reference regime (m =", opt$`ref-m`, ", n =", opt$`ref-n`, ").\n\n")
fits <- lapply(split(d, d$method), function(x) {
  # sub-millisecond timings carry no signal and log(0) is undefined; a method too fast to
  # measure is also a method nobody needs a budget for.
  x <- x[is.finite(x$seconds) & x$seconds > 0.005, ]
  if (nrow(x) < 3 || length(unique(x$m)) < 3) return(NULL)
  # only fit the n term when n actually varies, otherwise it is aliased with the intercept
  f <- if (length(unique(x$n)) > 1) stats::lm(log(seconds) ~ log(m) + log(n), data = x)
       else stats::lm(log(seconds) ~ log(m), data = x)
  cf <- coef(f)
  pr <- tryCatch(exp(unname(stats::predict(f, data.frame(m = opt$`ref-m`, n = opt$`ref-n`)))[1]),
                 error = function(e) NA_real_, warning = function(w) NA_real_)
  data.frame(method = x$method[1], n_points = nrow(x),
             slope_m = unname(cf["log(m)"]),
             slope_n = if ("log(n)" %in% names(cf)) unname(cf["log(n)"]) else NA_real_,
             max_observed_s = max(x$seconds), projected_s = pr, stringsAsFactors = FALSE)
})
fit <- do.call(rbind, Filter(Negate(is.null), fits))
if (is.null(fit) || !nrow(fit)) {
  cat("Not enough usable timings to fit a slope: need >= 3 distinct feature counts per method,\n",
      "each taking more than 5 ms. Widen --m, or restrict --methods to the slow ones.\n", sep = "")
} else {
  fit <- fit[order(-fit$projected_s, na.last = TRUE), ]
  fit$projected_min <- round(fit$projected_s / 60, 1)
  fit$axisA_cpu_h <- round(fit$projected_s * 13750 / 3600)
  print(format(fit[c("method", "n_points", "slope_m", "slope_n", "max_observed_s", "projected_min", "axisA_cpu_h")],
               digits = 3), row.names = FALSE)
  write.csv(fit, file.path(hpc, "probe_scaling.csv"), row.names = FALSE)
}
cat("\nwrote hpc/probe_runtime.csv and hpc/probe_scaling.csv\n")
cat("axisA_cpu_h projects one method across all 13750 Axis A cells at the reference regime.\n")
