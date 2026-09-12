#!/usr/bin/env Rscript
# fastemu_probe.R --------------------------------------------------------------
# Settles two things about the *installed* fastEmu before we commit grid budget:
#   (1) which reference-set argument this version actually accepts, and
#   (2) how its cost scales in the number of features.
#
# Background: the benchmark wrapper passes `reference_set_size = 30` only if that
# name appears in formals(fastEmuFit). On 2026-09-12 fastEmu and radEmu timed
# identically at matched feature counts (228.78 s vs 227.29 s), which means the
# reference set is doing no work -- either the argument is named something else,
# or it is being accepted and ignored. If fastEmu costs what radEmu costs, its
# share of the full grid is ~43 000 CPU-hours and it comes out of the run.
#
# Runtime: a few minutes. Run it on a login shell, not inside the grid.
# ------------------------------------------------------------------------------

ok <- requireNamespace("fastEmu", quietly = TRUE) && requireNamespace("radEmu", quietly = TRUE)
if (!ok) { cat("fastEmu and/or radEmu not installed -- nothing to probe.\n"); quit(status = 0) }

cat("radEmu  ", as.character(packageVersion("radEmu")),  "\n", sep = "")
cat("fastEmu ", as.character(packageVersion("fastEmu")), "\n\n", sep = "")

fn  <- get("fastEmuFit", envir = asNamespace("fastEmu"))
fml <- names(formals(fn))
cat("fastEmuFit formals:\n  ", paste(fml, collapse = ", "), "\n\n", sep = "")

ref_args <- grep("ref", fml, value = TRUE, ignore.case = TRUE)
cat("reference-set arguments in the signature: ",
    if (length(ref_args)) paste(ref_args, collapse = ", ") else "NONE",
    "\n", sep = "")
cat("exported functions mentioning 'emu':\n  ",
    paste(grep("emu", ls("package:fastEmu"), value = TRUE, ignore.case = TRUE), collapse = ", "),
    "\n\n", sep = "")

# radEmu renamed covariate_data -> data somewhere in 2.x; mirror the wrapper.
dat_arg  <- if ("data" %in% fml) "data" else "covariate_data"
has_tkj  <- "test_kj" %in% fml

sim <- function(n, J, seed = 1L) {
  set.seed(seed)
  Y <- matrix(rnbinom(n * J, mu = 20, size = 2), n, J,
              dimnames = list(paste0("s", seq_len(n)), paste0("t", seq_len(J))))
  Y[, seq_len(5)] <- matrix(rnbinom(n * 5, mu = 60, size = 2), n, 5)   # a little real signal
  Y[rowSums(Y) == 0, 1] <- 1L
  d <- data.frame(g = factor(rep(c("a", "b"), each = n / 2)), row.names = rownames(Y))
  list(Y = Y, d = d, J = J)
}

fit_last <- NULL
timed <- function(label, J, extra = list(), n = 60L) {
  s <- sim(n, J)
  a <- c(list(Y = s$Y, formula = ~ g, run_score_tests = TRUE),
         setNames(list(s$d), dat_arg), extra)
  if (has_tkj) a$test_kj <- data.frame(k = 2L, j = seq_len(J))
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch({ invisible(capture.output(f <- suppressMessages(do.call(fn, a)))); f },
                  error = function(e) e)
  el <- proc.time()[["elapsed"]] - t0
  if (inherits(res, "error")) {
    cat(sprintf("  %-34s J=%-4d  ERROR: %s\n", label, J,
                trimws(substr(conditionMessage(res), 1, 90))))
    return(invisible(NA_real_))
  }
  fit_last <<- res
  cat(sprintf("  %-34s J=%-4d  %7.1f s\n", label, J, el))
  invisible(el)
}

cat("A. default settings, scaling in J\n")
t40 <- timed("defaults", 40L)
t80 <- timed("defaults", 80L)
if (!is.na(t40) && !is.na(t80) && t40 > 0)
  cat(sprintf("   -> t(80)/t(40) = %.2f  (exponent %.2f; 1 = linear, 2 = quadratic)\n\n",
              t80 / t40, log(t80 / t40) / log(2)))

cat("B. what the returned object says about the reference set\n")
if (!is.null(fit_last)) {
  cat("   names(fit): ", paste(names(fit_last), collapse = ", "), "\n", sep = "")
  for (nm in grep("ref|constrain", names(fit_last), value = TRUE, ignore.case = TRUE)) {
    v <- fit_last[[nm]]
    cat(sprintf("   %s: %s\n", nm,
                paste(utils::head(as.character(unlist(v)), 12), collapse = " ")))
  }
} else cat("   (no successful default fit to inspect)\n")
cat("\n")

cat("C. with an explicit reference set, J = 80\n")
tried <- FALSE
if ("reference_set_size" %in% fml) {
  tried <- TRUE; tr <- timed("reference_set_size = 10", 80L, list(reference_set_size = 10L))
  if (!is.na(tr) && !is.na(t80)) cat(sprintf("   -> speedup vs defaults: %.2fx\n", t80 / tr))
}
if ("reference_set" %in% fml) {
  tried <- TRUE; tr <- timed("reference_set = 1:10", 80L, list(reference_set = seq_len(10)))
  if (!is.na(tr) && !is.na(t80)) cat(sprintf("   -> speedup vs defaults: %.2fx\n", t80 / tr))
}
if (!tried) cat("   no reference-set argument in the signature -- see line above for what is.\n")

cat("\nRead it like this:\n")
cat("  speedup ~1x, or no reference-set argument  -> fastEmu costs what radEmu costs; drop it.\n")
cat("  speedup >>1x                               -> wire that argument name into the wrapper.\n")
cat("  exponent near 2 with no working reference  -> the 500-feature cells are the whole budget.\n")
