#!/usr/bin/env Rscript
# fastemu_probe.R --------------------------------------------------------------
# What does fastEmu's reference set actually buy us, on the installed version?
#
# Run 1 (2026-09-12) answered the first half: fastEmu 2.0.1 / radEmu 2.3.2.0, and
# `reference_set` and `reference_set_size` are both real arguments -- so the wrapper's
# `reference_set_size = 30` IS being passed, and the earlier "fastEmu costs what radEmu costs"
# result cannot be explained by a dropped argument. That leaves two possibilities, which this
# run separates:
#   (a) the DEFAULT reference set is already small, so passing 30 changes nothing and fastEmu
#       is simply not much cheaper than radEmu at our feature counts; or
#   (b) the reference set works and the earlier comparison was confounded.
# Section A prints the defaults, section C times a fixed problem at reference sets of
# J, 30 and 5, and section D projects the per-cell cost to m = 500.
#
# Runtime: roughly 10 minutes. Run it on a login shell, not inside the grid.
# ------------------------------------------------------------------------------

if (!requireNamespace("fastEmu", quietly = TRUE) || !requireNamespace("radEmu", quietly = TRUE)) {
  cat("fastEmu and/or radEmu not installed -- nothing to probe.\n"); quit(status = 0)
}
cat("radEmu  ", as.character(packageVersion("radEmu")),  "\n", sep = "")
cat("fastEmu ", as.character(packageVersion("fastEmu")), "\n\n", sep = "")

fn  <- get("fastEmuFit", envir = asNamespace("fastEmu"))
fo  <- formals(fn); fml <- names(fo)

cat("A. defaults of the reference-set arguments\n")
for (nm in intersect(c("reference_set", "reference_set_size", "fastEmu_refit", "refit",
                       "penalize", "return_wald_p", "compute_cis"), fml)) {
  # an argument with no default is the empty symbol; evaluating it throws, so catch rather
  # than test -- `missing()` on a non-argument is not reliable inside every R version.
  txt <- tryCatch(paste(deparse(fo[[nm]]), collapse = " "), error = function(e) "<no default>")
  if (!nzchar(trimws(txt))) txt <- "<no default>"
  cat(sprintf("   %-20s = %s\n", nm, txt))
}
cat("   exported: ", paste(grep("emu", getNamespaceExports("fastEmu"), value = TRUE,
                                ignore.case = TRUE), collapse = ", "), "\n\n", sep = "")

dat_arg <- if ("data" %in% fml) "data" else "covariate_data"
has_tkj <- "test_kj" %in% fml

sim <- function(n, J, seed = 1L) {
  set.seed(seed)
  Y <- matrix(rnbinom(n * J, mu = 20, size = 2), n, J,
              dimnames = list(paste0("s", seq_len(n)), paste0("t", seq_len(J))))
  Y[, seq_len(5)] <- matrix(rnbinom(n * 5, mu = 60, size = 2), n, 5)
  Y[rowSums(Y) == 0, 1] <- 1L
  list(Y = Y, d = data.frame(g = factor(rep(c("a", "b"), each = n / 2)), row.names = rownames(Y)))
}

fits <- new.env()
timed <- function(label, J, extra = list(), n = 60L, keep = NULL) {
  s <- sim(n, J)
  a <- c(list(Y = s$Y, formula = ~ g, run_score_tests = TRUE), setNames(list(s$d), dat_arg), extra)
  if (has_tkj) a$test_kj <- data.frame(k = 2L, j = seq_len(J))
  t0 <- proc.time()[["elapsed"]]
  res <- tryCatch({ invisible(capture.output(f <- suppressMessages(do.call(fn, a)))); f },
                  error = function(e) e)
  el <- proc.time()[["elapsed"]] - t0
  if (inherits(res, "error")) {
    cat(sprintf("   %-26s J=%-4d  ERROR: %s\n", label, J,
                trimws(substr(conditionMessage(res), 1, 80))))
    return(invisible(NA_real_))
  }
  if (!is.null(keep)) assign(keep, res, envir = fits)
  cat(sprintf("   %-26s J=%-4d  %7.1f s\n", label, J, el))
  invisible(el)
}

cat("B. defaults, scaling in the number of features\n")
t40 <- timed("defaults", 40L)
t80 <- timed("defaults", 80L, keep = "def80")
if (!is.na(t40) && !is.na(t80) && t40 > 0)
  cat(sprintf("   -> t(80)/t(40) = %.2f, exponent %.2f  (1 = linear, 2 = quadratic)\n\n",
              t80 / t40, log(t80 / t40) / log(2)))

cat("C. reference-set size at fixed J = 80\n")
tall <- timed("reference_set_size = 80", 80L, list(reference_set_size = 80L))
t30  <- timed("reference_set_size = 30", 80L, list(reference_set_size = 30L), keep = "r30")
t05  <- timed("reference_set_size = 5",  80L, list(reference_set_size = 5L))
base <- if (!is.na(tall)) tall else t80
for (x in list(c(30, t30), c(5, t05)))
  if (!is.na(x[2]) && !is.na(base) && x[2] > 0)
    cat(sprintf("   -> size %2d is %.2fx the cost of the full reference set\n", x[1], x[2] / base))
cat("\n")

cat("D. what the fitted objects say\n")
if (!is.null(fits$def80)) {
  cat("   names(fit): ", paste(names(fits$def80), collapse = ", "), "\n", sep = "")
  for (nm in grep("ref|constrain", names(fits$def80), value = TRUE, ignore.case = TRUE))
    cat(sprintf("   default fit $%s: %s\n", nm,
                paste(utils::head(as.character(unlist(fits$def80[[nm]])), 12), collapse = " ")))
}
if (!is.null(fits$def80) && !is.null(fits$r30) &&
    !is.null(fits$def80$coef) && !is.null(fits$r30$coef)) {
  a <- fits$def80$coef$estimate; b <- fits$r30$coef$estimate
  if (length(a) == length(b))
    cat(sprintf("   estimates default vs size 30: max |diff| = %.4f  (0 = the argument did nothing)\n",
                max(abs(a - b), na.rm = TRUE)))
}
cat("\n")

cat("E. projection to a grid cell (m = 500 features)\n")
if (!is.na(t40) && !is.na(t80) && t40 > 0) {
  e <- log(t80 / t40) / log(2)
  for (tt in list(c("defaults", t80), c("size 30", t30)))
    if (!is.na(as.numeric(tt[2])))
      cat(sprintf("   %-9s  ~%.0f min per cell at m = 500  (exponent %.2f, 60 samples)\n",
                  tt[1], as.numeric(tt[2]) * (500 / 80)^e / 60, e))
  cat("   Grid is ~4 400 axis-A/B cells. Multiply.\n")
}
cat("\nDecision rule: if `size 30` is not at least ~4x cheaper than the full reference set,\n")
cat("fastEmu is radEmu with extra steps at our feature counts and comes out of the run.\n")
