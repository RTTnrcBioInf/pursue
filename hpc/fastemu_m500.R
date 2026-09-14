#!/usr/bin/env Rscript
# fastemu_m500.R ---------------------------------------------------------------
# Stop extrapolating and measure the number the budget actually turns on.
#
# We have two measured points at n = 100 -- 76.3 s at J = 80 and 167.5 s at J = 160 -- giving a
# local exponent of 1.13 and a projection of ~10 min for a real m = 500 cell. We also have an
# earlier direct measurement of 72.5 min at m = 500. Those cannot both be right: reaching 72.5
# min from 167.5 s needs an exponent of 2.86, not 1.13. One fit at the real size settles it.
#
# Runtime: somewhere between 10 and 75 minutes, which is exactly the question.
# ------------------------------------------------------------------------------
if (!requireNamespace("fastEmu", quietly = TRUE)) { cat("fastEmu not installed.\n"); quit(status = 0) }
fn <- get("fastEmuFit", envir = asNamespace("fastEmu")); fml <- names(formals(fn))
dat_arg <- if ("data" %in% fml) "data" else "covariate_data"

n <- 100L; J <- 500L
set.seed(1)
Y <- matrix(rnbinom(n * J, mu = 20, size = 2), n, J,
            dimnames = list(paste0("s", seq_len(n)), paste0("t", seq_len(J))))
Y[, seq_len(50)] <- matrix(rnbinom(n * 50, mu = 60, size = 2), n, 50)
Y[rowSums(Y) == 0, 1] <- 1L
d <- data.frame(g = factor(rep(c("a", "b"), each = n / 2)), row.names = rownames(Y))

cat(sprintf("fastEmu %s | one fit at n = %d, m = %d, defaults, started %s\n",
            as.character(packageVersion("fastEmu")), n, J, format(Sys.time(), "%H:%M:%S")))
a <- list(Y = Y, formula = ~ g, run_score_tests = TRUE,
          test_kj = data.frame(k = 2L, j = seq_len(J)))
a[[dat_arg]] <- d
t0 <- proc.time()[["elapsed"]]
f <- tryCatch({ invisible(capture.output(z <- suppressMessages(do.call(fn, a)))); z },
              error = function(e) e)
el <- proc.time()[["elapsed"]] - t0

if (inherits(f, "error")) {
  cat("ERROR after ", round(el / 60, 1), " min: ", conditionMessage(f), "\n", sep = "")
} else {
  cat(sprintf("\n  %.1f min for one m = 500 cell\n", el / 60))
  cat(sprintf("  reference set actually used: %d taxa\n",
              length(if (!is.null(z$reference_set)) z$reference_set else integer(0))))
  cat(sprintf("\n  Full grid (17 300 cells): %.0f CPU-hours, %.0f days at 44 cores.\n",
              el / 3600 * 17300, el / 3600 * 17300 / 44 / 24))
  cat(sprintf("  Proposed slice (~300 cells): %.0f CPU-hours, %.1f hours at 44 cores.\n",
              el / 3600 * 300, el / 3600 * 300 / 44))
  cat("\n  The sixteen other methods measured 1461 CPU-hours for all of axis A.\n")
}
