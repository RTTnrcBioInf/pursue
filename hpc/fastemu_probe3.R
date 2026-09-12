#!/usr/bin/env Rscript
# fastemu_probe3.R -------------------------------------------------------------
# Probe 2 (2026-09-12) settled the reference set: `reference_set_size` DEFAULTS to 30, which is
# exactly what our wrapper passes, so the argument has never done anything -- estimates matched
# to 4 dp and timings to 1.4%. Shrinking the set works but weakly (16x smaller buys 3.2x), so a
# large fixed cost sits outside it. Probe 2 also printed `fastEmu_refit = FALSE` as a default.
#
# That is the last plausible knob: if the null refits behind each score test can be done the
# fastEmu way rather than the radEmu way, the fixed cost is what moves. This probe times it at
# n = 100 (the R00 sample size, not probe 2's 60) and at two feature counts, and checks that it
# does not change the estimates -- a "speedup" that moves the coefficients is a different model,
# not a faster one.
#
# Runtime: roughly 15 minutes. Safe to run alongside the grid; it is single-threaded.
# ------------------------------------------------------------------------------
if (!requireNamespace("fastEmu", quietly = TRUE)) { cat("fastEmu not installed.\n"); quit(status = 0) }
fn <- get("fastEmuFit", envir = asNamespace("fastEmu")); fml <- names(formals(fn))
dat_arg <- if ("data" %in% fml) "data" else "covariate_data"
cat("fastEmu ", as.character(packageVersion("fastEmu")),
    " | fastEmu_refit in signature: ", "fastEmu_refit" %in% fml, "\n\n", sep = "")

sim <- function(n, J, seed = 1L) {
  set.seed(seed)
  Y <- matrix(rnbinom(n * J, mu = 20, size = 2), n, J,
              dimnames = list(paste0("s", seq_len(n)), paste0("t", seq_len(J))))
  Y[, seq_len(5)] <- matrix(rnbinom(n * 5, mu = 60, size = 2), n, 5)
  Y[rowSums(Y) == 0, 1] <- 1L
  list(Y = Y, d = data.frame(g = factor(rep(c("a", "b"), each = n / 2)), row.names = rownames(Y)))
}
res <- list()
timed <- function(key, J, extra = list(), n = 100L) {
  s <- sim(n, J)
  a <- c(list(Y = s$Y, formula = ~ g, run_score_tests = TRUE), setNames(list(s$d), dat_arg),
         list(test_kj = data.frame(k = 2L, j = seq_len(J))), extra)
  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch({ invisible(capture.output(z <- suppressMessages(do.call(fn, a)))); z },
                error = function(e) e)
  el <- proc.time()[["elapsed"]] - t0
  if (inherits(f, "error")) {
    cat(sprintf("  %-38s J=%-4d ERROR: %s\n", key, J, trimws(substr(conditionMessage(f), 1, 70))))
    return(invisible(NULL))
  }
  res[[paste0(key, "_", J)]] <<- list(sec = el, est = f$coef$estimate)
  cat(sprintf("  %-38s J=%-4d %7.1f s\n", key, J, el)); invisible(NULL)
}

cat("A. n = 100, J = 80\n")
timed("defaults",                           80L)
timed("fastEmu_refit = TRUE",               80L, list(fastEmu_refit = TRUE))
timed("fastEmu_refit = TRUE, ref size 5",   80L, list(fastEmu_refit = TRUE, reference_set_size = 5L))

cat("\nB. n = 100, J = 160  (does anything that helped still help?)\n")
timed("defaults",             160L)
timed("fastEmu_refit = TRUE", 160L, list(fastEmu_refit = TRUE))

cat("\nC. speedups and whether the estimand moved (J = 80)\n")
b <- res[["defaults_80"]]
for (k in grep("^fastEmu_refit.*_80$", names(res), value = TRUE)) {
  r <- res[[k]]
  d <- if (length(r$est) == length(b$est)) max(abs(r$est - b$est), na.rm = TRUE) else NA_real_
  cat(sprintf("   %-38s %.2fx faster | max |estimate diff| = %s\n",
              sub("_80$", "", k), b$sec / r$sec,
              if (is.na(d)) "length mismatch" else sprintf("%.4f", d)))
}

cat("\nD. projection to a grid cell (m = 500, n = 100)\n")
for (k in c("defaults", "fastEmu_refit = TRUE")) {
  a <- res[[paste0(k, "_80")]]; z <- res[[paste0(k, "_160")]]
  if (is.null(a) || is.null(z)) next
  e <- log(z$sec / a$sec) / log(2)
  cat(sprintf("   %-24s exponent %.2f -> ~%.0f min per cell; 17 300 cells = %.0f CPU-hours\n",
              k, e, z$sec * (500 / 160)^e / 60, z$sec * (500 / 160)^e / 3600 * 17300))
}
cat("\nThe other sixteen methods together are budgeted near 3 300 CPU-hours. Anything here that\n")
cat("does not land in the same order of magnitude means fastEmu runs on a declared slice of the\n")
cat("grid, not the whole of it.\n")
