#!/usr/bin/env Rscript
# Is fastEmu's data-driven reference set really collapsing to 2 taxa at large m, or did the
# m = 500 script measure the length of a list? `reference_set_size` defaults to 30, so a set of
# 2 would mean the estimand is "relative to two features" -- a fact about the method, not a
# detail. Fit at two sizes and print the structure rather than a length.
if (!requireNamespace("fastEmu", quietly = TRUE)) { cat("fastEmu not installed\n"); quit(status = 0) }
fn <- get("fastEmuFit", envir = asNamespace("fastEmu")); fml <- names(formals(fn))
dat_arg <- if ("data" %in% fml) "data" else "covariate_data"
for (J in c(100L, 400L)) {
  n <- 100L; set.seed(1)
  Y <- matrix(rnbinom(n * J, mu = 20, size = 2), n, J,
              dimnames = list(paste0("s", seq_len(n)), paste0("t", seq_len(J))))
  Y[, seq_len(5)] <- matrix(rnbinom(n * 5, mu = 60, size = 2), n, 5)
  Y[rowSums(Y) == 0, 1] <- 1L
  a <- list(Y = Y, formula = ~ g, run_score_tests = TRUE, test_kj = data.frame(k = 2L, j = seq_len(J)))
  a[[dat_arg]] <- data.frame(g = factor(rep(c("a", "b"), each = n / 2)), row.names = rownames(Y))
  t0 <- proc.time()[["elapsed"]]
  z <- tryCatch({ invisible(capture.output(o <- suppressMessages(do.call(fn, a)))); o },
                error = function(e) e)
  cat(sprintf("\n--- m = %d (%.1f min) ---\n", J, (proc.time()[["elapsed"]] - t0) / 60))
  if (inherits(z, "error")) { cat("ERROR:", conditionMessage(z), "\n"); next }
  cat("str(fit$reference_set):\n"); str(z$reference_set)
  cat("str(fit$reference_set_names):\n"); str(z$reference_set_names)
  cat(sprintf("taxa in reference set: %d of %d\n", length(unlist(z$reference_set)), J))
}
