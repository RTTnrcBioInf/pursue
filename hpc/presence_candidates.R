#!/usr/bin/env Rscript
# presence_candidates.R --------------------------------------------------------
# What should PURSUE 0.3 test? Three results pin it down:
#
#   1. The ablation (hpc/pursue_ablation.txt): switching off the ridge, the dispersion
#      shrinkage or both leaves the presence arm at 50-60% of plain logistic regression's
#      power. Tuning is not the problem.
#   2. Axis A R17-R19: under depth confounding (depth ratio x2/x4/x9 between groups) the
#      presence arm's FDR stays flat at 0.013-0.046, while logistic_presence -- which has no
#      depth term at all -- climbs to 0.658 averaged over simulators. The depth handling WORKS.
#   3. Axis C, real data with permuted labels: the presence arm's false positive rate is
#      0.018 against a nominal 0.05, with KS 0.35. Its null p-values pile up near 1. That is
#      the signature of a weakly identified zero-inflation logit: for any feature with no
#      structural zeros the group effect on it is unidentifiable, the LRT is ~0, and p ~ 1.
#
# So: keep the depth handling, drop the latent mixture. The candidate that does exactly that
# is a binomial GLM on DETECTION with a complementary log-log link and a log-depth offset.
# Under Poisson sampling it is not an approximation: P(X > 0) = 1 - exp(-lambda * N), so
# cloglog P(X > 0) = log(lambda) + log(N). The group coefficient is then the log fold change in
# relative abundance, estimated from the zeros -- which is where the abundance information
# lives at 16S sparsity, and exactly what the current abundance arm (E[X | X > 0], nearly
# independent of lambda) throws away. One test might replace both arms.
#
# Candidates, all scored on the SAME yardstick -- any differentially abundant feature
# (truth_abs), BH at 0.05 -- so their true- and false-positive counts compare directly:
#   pursue_presence   PURSUE 0.2 ZINB structural-presence arm (the current product)
#   pursue_combined   PURSUE 0.2 as shipped
#   logistic          detected ~ group                        (the one beating PURSUE; no depth)
#   logistic_depth    detected ~ group + log-depth             (logistic, made depth-aware)
#   cloglog_offset    cloglog: detected ~ group, offset log N  (the exact Poisson-sampling model)
#   cloglog_depth     cloglog: detected ~ group + log-depth    (log-depth free, to absorb overdispersion)
#
# Settings: power without confounding (axis-B implants B_ref, B_prev), robustness under depth
# confounding (house R00 / R17 / R18 / R19), calibration under the global null (house R06).
#
#   Rscript hpc/presence_candidates.R [--reps 10] [--template hmp_stool]
# Runtime: roughly 30-45 minutes, almost all of it the PURSUE fits.
# ------------------------------------------------------------------------------
suppressMessages(library(optparse))
op <- OptionParser(option_list = list(
  make_option("--reps", type = "integer", default = 10L),
  make_option("--template", default = "hmp_stool"),
  make_option("--out", default = "hpc/presence_candidates.csv"),
  make_option("--bench-root", default = NULL)))
opt <- parse_args(op)
if (!is.null(opt$`bench-root`)) Sys.setenv(PURSUE_BENCH_ROOT = opt$`bench-root`)
root <- Sys.getenv("PURSUE_BENCH_ROOT", unset = "benchmarks")
root <- tryCatch(normalizePath(root, mustWork = TRUE), error = function(e)
  stop("benchmarks directory not found at '", root, "' -- run from the repo root", call. = FALSE))
Sys.setenv(PURSUE_BENCH_ROOT = root)
if (!nzchar(Sys.getenv("PURSUE_DATA_ROOT"))) Sys.setenv(PURSUE_DATA_ROOT = file.path(root, "data"))
for (f in c("R/engine/templates.R", "R/engine/regimes.R", "R/engine/metrics.R",
            "R/simulators/sim_house.R", "R/simulators/implant.R", "R/simulators/dispatch.R",
            "R/methods/elementary.R", "R/methods/external.R", "R/methods/registry.R")) source(file.path(root, f))

tpl <- tryCatch(load_template(opt$template), error = function(e)
  stop("could not load template '", opt$template, "': ", conditionMessage(e), call. = FALSE))
cat(sprintf("template %s: %d features x %d samples\n\n", opt$template, nrow(tpl$counts), ncol(tpl$counts)))
regimes <- read.delim(file.path(root, "regimes.tsv"), stringsAsFactors = FALSE, comment.char = "#")

implant_spec <- function(signal) list(n_per_group = 50, da_frac = 0.10, effect = "medium", signal_type = signal,
                                      balance = "balanced", conf_phi = 0, exposure = "binary")
settings <- list(
  list(id = "B_ref  (power, mixed)",        make = function(s) implant(tpl, implant_spec("mixed"), s)),
  list(id = "B_prev (power, prevalence)",   make = function(s) implant(tpl, implant_spec("prevalence"), s)),
  list(id = "R00    (house, no confound)",  make = function(s) simulate_cell_data("house", tpl, regimes[regimes$regime_id == "R00", ], s)),
  list(id = "R17    (depth x2)",            make = function(s) simulate_cell_data("house", tpl, regimes[regimes$regime_id == "R17", ], s)),
  list(id = "R18    (depth x4)",            make = function(s) simulate_cell_data("house", tpl, regimes[regimes$regime_id == "R18", ], s)),
  list(id = "R19    (depth x9)",            make = function(s) simulate_cell_data("house", tpl, regimes[regimes$regime_id == "R19", ], s)),
  list(id = "R06    (global null)",         make = function(s) simulate_cell_data("house", tpl, regimes[regimes$regime_id == "R06", ], s)))

# One detection GLM per feature; LRT on the group term.
detect_test <- function(D, grp, ldepth, link, depth_mode) {
  n <- ncol(D); g <- as.numeric(grp); cd <- ldepth - mean(ldepth)
  fam <- stats::binomial(link = link)
  vapply(seq_len(nrow(D)), function(j) {
    d <- as.integer(D[j, ]); if (sum(d) < 3 || sum(d) > n - 3) return(NA_real_)
    X1 <- switch(depth_mode, none = cbind(1, g), covariate = cbind(1, g, cd), offset = cbind(1, g))
    X0 <- X1[, -2, drop = FALSE]
    off <- if (depth_mode == "offset") ldepth else rep(0, n)
    f1 <- tryCatch(suppressWarnings(stats::glm.fit(X1, d, family = fam, offset = off)), error = function(e) NULL)
    f0 <- tryCatch(suppressWarnings(stats::glm.fit(X0, d, family = fam, offset = off)), error = function(e) NULL)
    if (is.null(f1) || is.null(f0) || !f1$converged || !f0$converged) return(NA_real_)
    stats::pchisq(max(f0$deviance - f1$deviance, 0), 1, lower.tail = FALSE)
  }, numeric(1))
}

score <- function(p, truth) {
  ok <- is.finite(p); q <- rep(NA_real_, length(p)); q[ok] <- p.adjust(p[ok], "BH")
  rej <- !is.na(q) & q <= 0.05
  c(tested = sum(ok), n_pos = sum(truth[ok] == 1), rej = sum(rej),
    tp = sum(rej & truth == 1), fp = sum(rej & truth == 0),
    fpr = if (any(ok & truth == 0)) mean(p[ok & truth == 0] < 0.05) else NA)
}

rows <- list()
for (si in seq_along(settings)) for (rep in seq_len(opt$reps)) {
  st <- settings[[si]]
  sim <- tryCatch(st$make(5000L * si + rep), error = function(e) { message(st$id, " r", rep, ": ", conditionMessage(e)); NULL })
  if (is.null(sim) || !is.null(sim$unsupported)) next
  keep <- rowMeans(sim$counts > 0) >= 0.10 & rowSums(sim$counts > 0) >= 3
  ct <- sim$counts[keep, , drop = FALSE]
  truth <- sim$truth$truth_abs[match(rownames(ct), sim$truth$feature)]; truth[is.na(truth)] <- 0L
  depth <- if (!is.null(sim$meta$depth)) sim$meta$depth else colSums(sim$counts)
  grp <- sim$meta[[sim$tested_term]]; if (is.factor(grp)) grp <- as.integer(grp) - 1L
  D <- ct > 0; ld <- log(depth)

  P <- list(
    logistic       = detect_test(D, grp, ld, "logit",   "none"),
    logistic_depth = detect_test(D, grp, ld, "logit",   "covariate"),
    cloglog_offset = detect_test(D, grp, ld, "cloglog", "offset"),
    cloglog_depth  = detect_test(D, grp, ld, "cloglog", "covariate"))
  pr <- tryCatch(method_pursue(ct, sim$meta, sim$formula, sim$tested_term), error = function(e) NULL)
  if (!is.null(pr)) {
    P$pursue_presence <- pr$p[pr$arm == "presence"][match(rownames(ct), pr$feature[pr$arm == "presence"])]
    P$pursue_combined <- pr$p[pr$arm == "combined"][match(rownames(ct), pr$feature[pr$arm == "combined"])]
  }
  for (v in names(P)) rows[[length(rows) + 1L]] <-
    data.frame(setting = st$id, rep = rep, variant = v, t(score(P[[v]], truth)), stringsAsFactors = FALSE)
  cat(sprintf("  %-28s rep %2d done\n", st$id, rep))
}
D <- do.call(rbind, rows)
if (is.null(D)) { cat("no rows produced\n"); quit(status = 1) }
write.csv(D, opt$out, row.names = FALSE)

# FDR pooled over replicates (sum FP / sum rejections) is steadier than a mean of per-cell ratios
S <- do.call(rbind, lapply(split(D, list(D$setting, D$variant), drop = TRUE), function(x)
  data.frame(setting = x$setting[1], variant = x$variant[1], TP = mean(x$tp), FP = mean(x$fp),
             FDR = sum(x$fp) / max(sum(x$rej), 1), FPR = mean(x$fpr, na.rm = TRUE),
             n_pos = mean(x$n_pos), stringsAsFactors = FALSE)))
ord <- c("pursue_combined", "pursue_presence", "logistic", "logistic_depth", "cloglog_offset", "cloglog_depth")
for (s in unique(S$setting)) {
  a <- S[S$setting == s, ]; a <- a[order(match(a$variant, ord)), ]
  cat("\n", s, "   (DA features per cell: ", round(max(a$n_pos), 1), ")\n", sep = "")
  print(data.frame(variant = a$variant, TP = round(a$TP, 2), FP = round(a$FP, 2), FDR = round(a$FDR, 3),
                   FPR_nulls = round(a$FPR, 4)), row.names = FALSE)
}
cat("\nWhat a 0.3 presence test has to do, all at once:\n")
cat("  B_ref / B_prev   TP close to `logistic`                     (the power PURSUE is missing)\n")
cat("  R17 / R18 / R19  FDR stays near 0.05 as depth confounding grows (what `logistic` cannot do)\n")
cat("  R06              FPR_nulls near 0.05, not 0.018             (calibrated, not conservative)\n")
cat("\nwrote ", opt$out, "\n", sep = "")
