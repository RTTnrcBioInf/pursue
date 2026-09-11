#!/usr/bin/env Rscript
# =============================================================================
# Step 1 -- which inference mode for the abundance arm?
#
# Simulates the abundance-arm response directly (log-scale linear model with a
# confounder) so that ONLY the inference layer differs between methods. Every
# method sees the same per-taxon design and the same fitted coefficients.
#
#   ols    classical t, df = n_j - p
#   modt   limma moderated t (empirical-Bayes variance shrinkage)
#   fl_B   Freedman-Lane residual permutation, empirical p, B permutations
#   pvw_B  permutation-variance Wald (LOCOM2-style): z = beta / sd(beta*), N(0,1) tail,
#          beta* from Freedman-Lane (reduced-model residuals)
#   pvt_B  same, but beta* from ter Braak permutation (full-model residuals), so the
#          null SD is not inflated by the effect under H1
#   mm_B   moment-matched permutation null: Pearson type III fitted to the first
#          three moments of F* = t*^2 from B permutations (E-MANOVA / PERMANOVA-S)
#
# Outputs one long CSV of per-replicate metrics.
# =============================================================================

suppressMessages({ library(limma); library(parallel) })

args <- commandArgs(trailingOnly = TRUE)
out_csv   <- if (length(args) >= 1) args[1] else "step1_results.csv"
quick     <- length(args) >= 2 && args[2] == "quick"

set.seed(20260909)

# ---------------------------------------------------------------- settings ---
m        <- if (quick) 100 else 500          # taxa
pi1      <- 0.10                             # fraction DA in non-null cells
delta    <- 1.0                              # effect size in units of sigma_j
rho_xz   <- 0.5                              # confounder-exposure correlation
B_fl     <- if (quick) c(199, 999) else c(199, 999, 9999)
B_small  <- c(199, 999)                      # for pvw / mm
alpha_bh <- 0.05
n_rep_null <- if (quick) 3 else 20
n_rep_alt  <- if (quick) 3 else 15
n_grid   <- if (quick) c(20, 40) else c(10, 20, 40, 100)   # total n
errors   <- c("gauss", "t3", "contam", "skew", "hetero", "missing")
n_cores  <- max(1L, detectCores() - 0L)

# ------------------------------------------------------ error generators ---
rerr <- function(n, kind) {
  switch(kind,
    gauss   = rnorm(n),
    t3      = rt(n, df = 3) / sqrt(3),
    contam  = { o <- runif(n) < 0.10; (ifelse(o, rnorm(n, 0, 4), rnorm(n))) / sqrt(2.5) },
    skew    = (rgamma(n, shape = 2, rate = 1) - 2) / sqrt(2),
    hetero  = rnorm(n),           # scaled by group later
    missing = rnorm(n)
  )
}

# ------------------------------------------------------------- simulate ---
simulate_cell <- function(n, kind, null_cell) {
  # exposure x (binary) and confounder z, correlated
  u <- rnorm(n)
  unbalanced <- kind == "hetero"
  x <- if (unbalanced) as.integer(u > qnorm(0.70)) else as.integer(u > 0)
  z <- rho_xz * u + sqrt(1 - rho_xz^2) * rnorm(n)

  # per-taxon variances: scaled inverse chi-square, d0 = 4, s0 = 1 (limma's prior form)
  d0 <- 4; s0sq <- 1
  sigma2 <- s0sq * d0 / rchisq(m, df = d0)
  sigma  <- sqrt(sigma2)

  gamma <- rnorm(m, 0.5, 0.2)                     # confounder effect on every taxon
  is_da <- rep(FALSE, m)
  beta  <- numeric(m)
  if (!null_cell) {
    is_da[sample.int(m, round(pi1 * m))] <- TRUE
    beta[is_da] <- delta * sigma[is_da] * sample(c(-1, 1), sum(is_da), replace = TRUE)
  }

  Y <- matrix(NA_real_, n, m)
  for (j in seq_len(m)) {
    e <- rerr(n, kind)
    if (kind == "hetero") e <- e * ifelse(x == 1, sqrt(3), 1)   # small group, big variance
    Y[, j] <- gamma[j] * z + beta[j] * x + sigma[j] * e
  }
  if (kind == "missing") {
    miss <- matrix(runif(n * m) < 0.40, n, m)
    Y[miss] <- NA_real_
  }
  list(Y = Y, x = x, z = z, is_da = is_da)
}

# ------------------------------------------------- per-taxon test engine ---
pearson3_p <- function(Fobs, Fstar) {
  mu <- mean(Fstar); v <- var(Fstar)
  if (!is.finite(v) || v <= 0) return(1)
  g <- mean((Fstar - mu)^3) / v^1.5
  if (is.finite(g) && g > 0.05) {
    k <- 4 / g^2; th <- sqrt(v / k); c0 <- mu - k * th
    if (Fobs - c0 <= 0) return(1)
    pgamma(Fobs - c0, shape = k, scale = th, lower.tail = FALSE)
  } else {
    pnorm((Fobs - mu) / sqrt(v), lower.tail = FALSE)
  }
}

test_all <- function(y, x, z, perm_global) {
  ok <- is.finite(y)
  y <- y[ok]; xx <- x[ok]; zz <- z[ok]
  n <- length(y); p <- 3L; df <- n - p
  if (df < 2L || length(unique(xx)) < 2L) return(NULL)

  X_full <- cbind(1, zz, xx); X_red <- cbind(1, zz)
  qf <- qr(X_full); qr_ <- qr(X_red)
  coef <- qr.coef(qf, y); res <- qr.resid(qf, y)
  s2   <- sum(res^2) / df
  XtXi <- chol2inv(qr.R(qf))
  vxx  <- XtXi[3, 3]
  beta <- unname(coef[3])
  t_obs <- beta / sqrt(s2 * vxx)

  yhat_red <- qr.fitted(qr_, y); e_red <- qr.resid(qr_, y)
  w <- (XtXi %*% t(X_full))[3, ]                      # beta* = w' y*

  Bmax <- max(B_fl)
  perm <- if (all(ok)) perm_global else apply(perm_global[ok, , drop = FALSE], 2, rank)  # induced local perms
  Ystar <- yhat_red + e_red[perm]                     # n x Bmax (recycled columns)
  dim(Ystar) <- c(n, Bmax)
  beta_star <- as.vector(w %*% Ystar)
  rss_star  <- colSums(qr.resid(qf, Ystar)^2)
  t_star    <- beta_star / sqrt(rss_star / df * vxx)

  # ter Braak: permute FULL-model residuals around the full-model fit
  Bs <- max(B_small)
  yhat_full <- y - res
  Ytb <- yhat_full + res[perm[, seq_len(Bs)]]
  dim(Ytb) <- c(n, Bs)
  beta_tb <- as.vector(w %*% Ytb)

  out <- c(beta = beta, s2 = s2, df = df, t = t_obs,
           p_ols = 2 * pt(-abs(t_obs), df))
  for (B in B_fl) {
    ts <- t_star[seq_len(B)]
    out[paste0("p_fl_", B)] <- (1 + sum(abs(ts) >= abs(t_obs))) / (B + 1)
  }
  for (B in B_small) {
    bs <- beta_star[seq_len(B)]; ts <- t_star[seq_len(B)]
    zpv <- beta / sd(bs)
    out[paste0("p_pvw_", B)] <- 2 * pnorm(-abs(zpv))
    out[paste0("p_pvt_", B)] <- 2 * pnorm(-abs(beta / sd(beta_tb[seq_len(B)])))
    out[paste0("p_mm_", B)]  <- pearson3_p(t_obs^2, ts^2)
  }
  out
}

run_rep <- function(n, kind, null_cell, rep_id) {
  sim <- simulate_cell(n, kind, null_cell)
  perm_global <- replicate(max(B_fl), sample.int(n))     # one permutation set per replicate, shared by all taxa
  res <- lapply(seq_len(m), function(j) test_all(sim$Y[, j], sim$x, sim$z, perm_global))
  keep <- !vapply(res, is.null, logical(1))
  R <- do.call(rbind, res[keep]); is_da <- sim$is_da[keep]
  # limma moderated t on the same fits
  sq <- squeezeVar(R[, "s2"], df = R[, "df"])
  # recover vxx from t and s2: t = beta / sqrt(s2*vxx)  ->  vxx = (beta/t)^2 / s2
  vxx <- (R[, "beta"] / R[, "t"])^2 / R[, "s2"]
  t_mod <- R[, "beta"] / sqrt(sq$var.post * vxx)
  p_modt <- 2 * pt(-abs(t_mod), df = R[, "df"] + sq$df.prior)

  P <- cbind(ols = R[, "p_ols"], modt = p_modt,
             R[, grep("^p_(fl|pvw|pvt|mm)_", colnames(R)), drop = FALSE])
  colnames(P) <- sub("^p_", "", colnames(P))

  rows <- lapply(colnames(P), function(meth) {
    pv <- P[, meth]; q <- p.adjust(pv, "BH"); rej <- q <= alpha_bh
    tp <- sum(rej & is_da); fp <- sum(rej & !is_da); nd <- sum(is_da)
    pn <- pv[!is_da]
    data.frame(n = n, error = kind, cell = if (null_cell) "null" else "alt", rep = rep_id,
               method = meth, n_taxa = nrow(P),
               fdr = if (tp + fp > 0) fp / (tp + fp) else 0,
               n_rej = tp + fp,
               tpr = if (nd > 0) tp / nd else NA_real_,
               fpr05 = mean(pn < 0.05), fpr01 = mean(pn < 0.01),
               ks = suppressWarnings(ks.test(pn, "punif")$statistic),
               p_min = min(pv), d0 = sq$df.prior,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# ------------------------------------------------------------------ grid ---
grid <- expand.grid(n = n_grid, error = errors, cell = c("null", "alt"), stringsAsFactors = FALSE)
grid <- grid[!(grid$error == "missing" & grid$n < 40), ]     # two-part missingness needs n>=40
tasks <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  nr <- if (grid$cell[i] == "null") n_rep_null else n_rep_alt
  data.frame(grid[i, ], rep = seq_len(nr), stringsAsFactors = FALSE)
}))
cat(sprintf("m=%d taxa, %d tasks, %d cores\n", m, nrow(tasks), n_cores))

t0 <- Sys.time()
res <- mclapply(seq_len(nrow(tasks)), function(i) {
  tk <- tasks[i, ]
  set.seed(1e6 + i)
  tryCatch(run_rep(tk$n, tk$error, tk$cell == "null", tk$rep),
           error = function(e) { message("task ", i, " failed: ", conditionMessage(e)); NULL })
}, mc.cores = n_cores)
res <- do.call(rbind, res)
write.csv(res, out_csv, row.names = FALSE)
cat(sprintf("done in %.1f min -> %s (%d rows)\n",
            as.numeric(Sys.time() - t0, units = "mins"), out_csv, nrow(res)))
