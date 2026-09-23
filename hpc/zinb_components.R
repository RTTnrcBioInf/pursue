#!/usr/bin/env Rscript
# zinb_components.R ------------------------------------------------------------
# presence_candidates.R overturned the plan to replace PURSUE's presence arm with a detection
# GLM. On house R17/R18/R19 (depth ratio x2/x4/x9 between groups):
#
#                       R00 TP   R17 TP   R18 TP   R19 TP     R19 FDR
#   logistic + log-depth  7.0      4.4      1.0      0.2       0.500
#   PURSUE presence arm   3.3      4.0      3.9      4.3       0.085
#
# The detection GLM collapses under confounding; PURSUE's ZINB keeps its power. The reason is
# in the model. PURSUE's ZINB is
#     logit psi_i = X_i' alpha                  (structural presence)
#     log  mu_i   = log N_i + X_i' gamma        (abundance when present, NB with dispersion theta)
# and it gets the detection-depth relationship right BY CONSTRUCTION: P(X=0) follows from
# NB(0 | mu, theta) with mu proportional to depth, which saturates correctly under overdispersion.
# It estimates that relationship from the magnitudes of the nonzero counts, so it can tell a
# depth effect from a group effect even when the two are collinear. A detection GLM has only
# the binary pattern to work with, and there the two cannot be separated.
#
# And the group term is in BOTH components -- but PURSUE 0.2 only tests alpha. It fits gamma,
# the group effect on abundance, with a depth offset, from every sample including the zeros, and
# then discards it; its reported abundance result comes from a separate LM on detected cells
# only, which recovers 0.3% of true positives.
#
# This script takes PURSUE's presence test verbatim (PURSUE:::fit_presence_arm, same theta),
# then refits the same ZINB at the same theta with constraints, adding:
#   zinb_mu      LRT on gamma for the tested term  (abundance, zeros included, depth offset)
#   zinb_joint   LRT on alpha and gamma together   (the "any effect" test, df = 2k)
# and, as a check that the refit is faithful, its own psi test, whose p-values should track
# PURSUE's own presence p-values almost exactly.
#
# Settings as presence_candidates.R, plus B_abund (abundance-only truth), where the presence
# test has nothing to find and the abundance test has everything to find.
#
#   Rscript hpc/zinb_components.R [--reps 10] [--cores 16] [--template hmp_stool]
# Runtime: ~2.5x a PURSUE fit per cell; with --cores 16 roughly 15-30 minutes.
# ------------------------------------------------------------------------------
suppressMessages(library(optparse))
op <- OptionParser(option_list = list(
  make_option("--reps", type = "integer", default = 10L),
  make_option("--cores", type = "integer", default = 8L),
  make_option("--template", default = "hmp_stool"),
  make_option("--out", default = "hpc/zinb_components.csv"),
  make_option("--bench-root", default = NULL)))
opt <- parse_args(op)
if (!is.null(opt$`bench-root`)) Sys.setenv(PURSUE_BENCH_ROOT = opt$`bench-root`)
root <- Sys.getenv("PURSUE_BENCH_ROOT", unset = "benchmarks")
root <- tryCatch(normalizePath(root, mustWork = TRUE), error = function(e)
  stop("benchmarks directory not found at '", root, "' -- run from the repo root", call. = FALSE))
Sys.setenv(PURSUE_BENCH_ROOT = root)
if (!nzchar(Sys.getenv("PURSUE_DATA_ROOT"))) Sys.setenv(PURSUE_DATA_ROOT = file.path(root, "data"))
for (f in c("R/engine/templates.R", "R/engine/regimes.R", "R/engine/metrics.R",
            "R/simulators/sim_house.R", "R/simulators/implant.R", "R/simulators/dispatch.R")) source(file.path(root, f))
stopifnot(requireNamespace("PURSUE", quietly = TRUE))
cat("PURSUE ", as.character(packageVersion("PURSUE")), ", ", opt$cores, " cores\n", sep = "")

tpl <- tryCatch(load_template(opt$template), error = function(e)
  stop("could not load template '", opt$template, "': ", conditionMessage(e), call. = FALSE))
cat(sprintf("template %s: %d features x %d samples\n\n", opt$template, nrow(tpl$counts), ncol(tpl$counts)))
regimes <- read.delim(file.path(root, "regimes.tsv"), stringsAsFactors = FALSE, comment.char = "#")
spec <- function(signal) list(n_per_group = 50, da_frac = 0.10, effect = "medium", signal_type = signal,
                              balance = "balanced", conf_phi = 0, exposure = "binary")
house <- function(id) function(s) simulate_cell_data("house", tpl, regimes[regimes$regime_id == id, ], s)
settings <- list(
  list(id = "B_abund (abundance-only truth)", make = function(s) implant(tpl, spec("abundance"), s)),
  list(id = "B_prev  (prevalence-only truth)", make = function(s) implant(tpl, spec("prevalence"), s)),
  list(id = "B_ref   (mixed truth)",          make = function(s) implant(tpl, spec("mixed"), s)),
  list(id = "R00     (house, no confound)",   make = house("R00")),
  list(id = "R17     (depth x2)",             make = house("R17")),
  list(id = "R18     (depth x4)",             make = house("R18")),
  list(id = "R19     (depth x9)",             make = house("R19")),
  list(id = "R06     (global null)",          make = house("R06")))

# ---- the same ZINB as PURSUE:::fit_presence_arm, with a mask on each component -----------
# `gfix` holds the value a constrained abundance coefficient is held at: 0 for the raw test,
# the empirical-null centre delta for the centred one.
zinb_make <- function(xj, lN, X, theta, tau = 5) {
  D <- xj > 0; p <- ncol(X)
  nll <- function(par, fz, fc, gfix) {
    a <- numeric(p); a[fz] <- par[seq_len(sum(fz))]
    g <- gfix; g[fc] <- par[sum(fz) + seq_len(sum(fc))]
    psi <- stats::plogis(drop(X %*% a)); mu <- exp(lN + drop(X %*% g))
    l0 <- exp(theta * (log(theta) - log(theta + mu)))
    P0 <- (1 - psi) + psi * l0
    ll <- sum(log(pmax(P0[!D], 1e-300))) +
          sum(log(pmax(psi[D], 1e-300)) + stats::dnbinom(xj[D], size = theta, mu = mu[D], log = TRUE))
    -ll + sum(a^2) / (2 * tau^2)
  }
  best <- function(fz, fc, gfix, starts) {
    lo <- c(rep(-12, sum(fz)), ifelse(which(fc) == 1L, -30, -10)); hi <- c(rep(12, sum(fz)), rep(10, sum(fc)))
    out <- NULL
    for (st in starts) {
      f <- tryCatch(stats::nlminb(st, nll, fz = fz, fc = fc, gfix = gfix, lower = lo, upper = hi), error = function(e) NULL)
      if (!is.null(f) && is.finite(f$objective) && (is.null(out) || f$objective < out$objective)) out <- f
    }
    out
  }
  list(nll = nll, best = best, D = D, p = p)
}
lrt_p <- function(null, full_obj, df) {
  if (is.null(null)) return(NA_real_)
  s <- 2 * (null$objective - full_obj)
  if (!is.finite(s) || s < -1e-6) return(NA_real_)
  stats::pchisq(max(s, 0), df, lower.tail = FALSE)
}
wig <- function(v, i, d) { v[i] <- v[i] + d; v }

# Pass 1, per feature: full fit; the presence test; the RAW abundance and joint tests (gamma held
# at 0); and gamma-hat with its SE, which the cross-feature centring needs.
zinb_pass1 <- function(xj, lN, X, tcols, theta, tau = 5) {
  n <- length(xj); D <- xj > 0; n1 <- sum(D); p <- ncol(X); k <- length(tcols)
  if (n1 < 3L || n1 > n - 3L || !is.finite(theta)) return(NULL)
  Z <- zinb_make(xj, lN, X, theta, tau)
  a0 <- rep(0, p); a0[1L] <- stats::qlogis(min(max(mean(D), 0.05), 0.95))
  g0 <- rep(0, p); g0[1L] <- mean(log(xj[D]) - lN[D])
  if (n1 >= p + 3L) { cf <- tryCatch(stats::lm.fit(X[D, , drop = FALSE], log(xj[D]) - lN[D])$coefficients, error = function(e) NULL)
    if (!is.null(cf) && all(is.finite(cf))) g0 <- unname(cf) }
  allz <- rep(TRUE, p); allc <- rep(TRUE, p); z0 <- rep(0, p)
  s0 <- c(a0, g0)
  full <- Z$best(allz, allc, z0, list(s0, wig(s0, tcols[1L], -1.5), wig(s0, tcols[1L], 1.5)))
  if (is.null(full)) return(NULL)
  af <- full$par[seq_len(p)]; gf <- full$par[p + seq_len(p)]
  nz <- allz; nz[tcols] <- FALSE; nc <- allc; nc[tcols] <- FALSE
  n_psi <- Z$best(nz, allc, z0, list(c(af[nz], gf), wig(c(af[nz], gf), sum(nz) + 1L, 1)))
  n_mu  <- Z$best(allz, nc, z0, list(c(af, gf[nc]), wig(c(af, gf[nc]), p + 1L, 1)))
  n_bo  <- Z$best(nz, nc,   z0, list(c(af[nz], gf[nc]), wig(c(af[nz], gf[nc]), sum(nz) + 1L, 1)))
  se <- rep(NA_real_, k)
  H <- tryCatch(stats::optimHess(full$par, Z$nll, fz = allz, fc = allc, gfix = z0), error = function(e) NULL)
  if (!is.null(H)) { V <- tryCatch(solve(H), error = function(e) NULL)
    if (!is.null(V)) { d <- diag(V)[p + tcols]; se <- ifelse(is.finite(d) & d > 0, sqrt(d), NA_real_) } }
  list(obj = full$objective, af = af, gf = gf, g_hat = gf[tcols], g_se = se,
       p_psi = lrt_p(n_psi, full$objective, k), p_mu_raw = lrt_p(n_mu, full$objective, k),
       p_joint_raw = lrt_p(n_bo, full$objective, 2L * k))
}

# Pass 2, per feature: the same constrained fits, with gamma held at the empirical-null centre
# delta instead of 0. The test becomes "does this feature differ from the TYPICAL feature" --
# PURSUE's own estimand -- which removes the compositional shift that the observed library size
# (itself raised by whatever is truly changing) otherwise pushes onto every null feature.
zinb_pass2 <- function(xj, lN, X, tcols, theta, r1, delta, tau = 5) {
  if (is.null(r1)) return(c(mu = NA_real_, joint = NA_real_))
  p <- ncol(X); k <- length(tcols); Z <- zinb_make(xj, lN, X, theta, tau)
  allz <- rep(TRUE, p); allc <- rep(TRUE, p); nz <- allz; nz[tcols] <- FALSE; nc <- allc; nc[tcols] <- FALSE
  gfix <- rep(0, p); gfix[tcols] <- delta
  af <- r1$af; gf <- r1$gf
  n_mu <- Z$best(allz, nc, gfix, list(c(af, gf[nc]), wig(c(af, gf[nc]), p + 1L, 1)))
  n_bo <- Z$best(nz, nc,   gfix, list(c(af[nz], gf[nc]), wig(c(af[nz], gf[nc]), sum(nz) + 1L, 1)))
  c(mu = lrt_p(n_mu, r1$obj, k), joint = lrt_p(n_bo, r1$obj, 2L * k))
}

logistic_p <- function(D, grp) {
  vapply(seq_len(nrow(D)), function(j) { d <- as.integer(D[j, ]); if (sum(d) < 3 || sum(d) > length(d) - 3) return(NA_real_)
    f1 <- suppressWarnings(stats::glm.fit(cbind(1, grp), d, family = stats::binomial()))
    f0 <- suppressWarnings(stats::glm.fit(matrix(1, length(d)), d, family = stats::binomial()))
    stats::pchisq(max(f0$deviance - f1$deviance, 0), 1, lower.tail = FALSE) }, numeric(1))
}
score <- function(p, truth) {
  ok <- is.finite(p); q <- rep(NA_real_, length(p)); q[ok] <- p.adjust(p[ok], "BH"); rej <- !is.na(q) & q <= 0.05
  c(tested = sum(ok), n_pos = sum(truth[ok] == 1), rej = sum(rej), tp = sum(rej & truth == 1),
    fp = sum(rej & truth == 0), fpr = if (any(ok & truth == 0)) mean(p[ok & truth == 0] < 0.05) else NA)
}

rows <- list(); parity <- c(); deltas <- NULL
enull <- get("enull_center", envir = asNamespace("PURSUE"))   # PURSUE 0.2's own empirical-null centring
for (si in seq_along(settings)) for (rep in seq_len(opt$reps)) {
  st <- settings[[si]]
  sim <- tryCatch(st$make(7000L * si + rep), error = function(e) { message(st$id, " r", rep, ": ", conditionMessage(e)); NULL })
  if (is.null(sim) || !is.null(sim$unsupported)) next
  keep <- rowMeans(sim$counts > 0) >= 0.10 & rowSums(sim$counts > 0) >= 3
  ct <- sim$counts[keep, , drop = FALSE]; otu <- t(ct)
  truth <- sim$truth$truth_abs[match(rownames(ct), sim$truth$feature)]; truth[is.na(truth)] <- 0L
  meta <- sim$meta[rownames(otu), , drop = FALSE]
  depth <- rowSums(otu); lN <- log(depth)
  design <- PURSUE:::build_design(sim$formula, meta, sim$tested_term, depth)
  pr <- tryCatch(PURSUE:::fit_presence_arm(otu, depth, design, n_cores = opt$cores), error = function(e) NULL)
  if (is.null(pr)) { message(st$id, " r", rep, ": presence arm failed"); next }
  r1 <- parallel::mclapply(seq_len(ncol(otu)), function(j)
          tryCatch(zinb_pass1(otu[, j], lN, design$X, design$tested_cols, pr$pres_theta[j]), error = function(e) NULL),
        mc.cores = opt$cores)
  gh <- vapply(r1, function(r) if (is.null(r)) NA_real_ else r$g_hat[1L], numeric(1))
  gs <- vapply(r1, function(r) if (is.null(r)) NA_real_ else r$g_se[1L],  numeric(1))
  ce <- enull(gh, gs); delta <- ce$delta
  r2 <- parallel::mclapply(seq_len(ncol(otu)), function(j)
          tryCatch(zinb_pass2(otu[, j], lN, design$X, design$tested_cols, pr$pres_theta[j], r1[[j]], delta),
                   error = function(e) c(mu = NA_real_, joint = NA_real_)), mc.cores = opt$cores)
  r2 <- do.call(rbind, r2)
  pick <- function(nm) vapply(r1, function(r) if (is.null(r)) NA_real_ else r[[nm]], numeric(1))
  own_psi <- pick("p_psi")
  ok <- is.finite(own_psi) & is.finite(pr$pres_p)
  if (sum(ok) > 10) parity <- c(parity, stats::cor(-log10(own_psi[ok]), -log10(pr$pres_p[ok]), method = "spearman"))
  deltas <- rbind(deltas, data.frame(setting = st$id, rep = rep, delta = delta, pi0 = ce$pi0))
  grp <- meta[[sim$tested_term]]; grp <- if (is.factor(grp)) as.integer(grp) - 1L else grp
  P <- list(pursue_presence = pr$pres_p,
            zinb_mu_raw = pick("p_mu_raw"), zinb_mu = r2[, "mu"],
            zinb_joint_raw = pick("p_joint_raw"), zinb_joint = r2[, "joint"],
            logistic = logistic_p(t(otu > 0), grp))
  for (v in names(P)) rows[[length(rows) + 1L]] <-
    data.frame(setting = st$id, rep = rep, variant = v, t(score(P[[v]], truth)), stringsAsFactors = FALSE)
  cat(sprintf("  %-34s rep %2d done\n", st$id, rep))
}
D <- do.call(rbind, rows)
if (is.null(D)) { cat("no rows produced\n"); quit(status = 1) }
write.csv(D, opt$out, row.names = FALSE)

cat(sprintf("\nrefit fidelity: Spearman(own psi-test, PURSUE presence p) median %.3f, min %.3f over %d cells\n",
            stats::median(parity), min(parity), length(parity)))
cat("  (near 1 means the mu and joint tests below come from the same model PURSUE fits)\n")
S <- do.call(rbind, lapply(split(D, list(D$setting, D$variant), drop = TRUE), function(x)
  data.frame(setting = x$setting[1], variant = x$variant[1], TP = mean(x$tp), FP = mean(x$fp),
             FDR = sum(x$fp) / max(sum(x$rej), 1), FPR = mean(x$fpr, na.rm = TRUE), n_pos = mean(x$n_pos),
             stringsAsFactors = FALSE)))
ord <- c("pursue_presence", "zinb_mu_raw", "zinb_mu", "zinb_joint_raw", "zinb_joint", "logistic")
if (!is.null(deltas)) {
  dd <- aggregate(cbind(delta, pi0) ~ setting, deltas, mean)
  cat("\nempirical-null centre of the abundance effect, by setting (0 = no compositional shift):\n")
  print(data.frame(setting = dd$setting, delta = round(dd$delta, 3), pi0 = round(dd$pi0, 3)), row.names = FALSE)
}
for (s in unique(S$setting)) {
  a <- S[S$setting == s, ]; a <- a[order(match(a$variant, ord)), ]
  cat("\n", s, "   (DA features per cell: ", round(max(a$n_pos), 1), ")\n", sep = "")
  print(data.frame(variant = a$variant, TP = round(a$TP, 2), FP = round(a$FP, 2), FDR = round(a$FDR, 3),
                   FPR_nulls = round(a$FPR, 4)), row.names = FALSE)
}
cat("\n*_raw hold the abundance coefficient at 0; the unsuffixed versions hold it at the empirical-null\n")
cat("centre. A synthetic check showed the raw tests reaching FPR 0.18-0.31 whenever OTHER features\n")
cat("changed in abundance, while staying calibrated on a global null: compositional bias through the\n")
cat("depth offset. Centring is meant to remove exactly that.\n")
cat("\nWhat would make zinb_joint the 0.3 test:\n")
cat("  B_abund          it finds what the presence test cannot (the 0.2 abundance arm finds ~0)\n")
cat("  B_ref / R00      TP near logistic without logistic's behaviour under confounding\n")
cat("  R17 / R18 / R19  keeps PURSUE's power under depth confounding, FDR near 0.05\n")
cat("  R06              FPR_nulls near 0.05\n")
cat("\nwrote ", opt$out, "\n", sep = "")
