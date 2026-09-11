#!/usr/bin/env Rscript
# =============================================================================
# Step 4 -- a stabilized joint zero-inflated fit for the prevalence arm
#
# Same generator as step 2 (structural presence x lognormal abundance x Poisson/NB
# sampling at depth N, depth confounded by c_depth) plus a pure-NB cell (no lognormal
# mixing) so misspecification can be told apart from small-sample effects.
#
# Joint model per taxon (jzi):
#   logit psi_i = a + b x_i                       (structural presence)
#   log mu_i    = log N_i + c + beta x_i          (abundance when present, NB(mu, theta))
#   P(X=0) = (1-psi) + psi NB(0|mu,theta);  P(X=x>0) = psi NB(x|mu,theta)
# Stabilizers:
#   * theta: per-taxon estimate from the zero-truncated NB on detected cells (closed form),
#     shrunk on the log scale toward a trend across taxa (loess of log theta on log mean
#     abundance) with prior weight d0 pseudo-observations, then FIXED for inference.
#   * zero part: ridge penalty (a^2 + b^2) / (2 tau^2), tau = 5 (weak; prevents separation).
#   * LRT for b (prevalence) and for beta (abundance), each by refitting the constrained
#     model from the full-model solution; multiple starts.
# Comparators: naive logistic, logistic + logN, occ_ztnb (two-stage, step 2b),
#              pscl ZINB (n >= 50), and for the abundance arm lm(log r ~ x + logN).
# =============================================================================
suppressMessages({ library(parallel); library(pscl) })
args <- commandArgs(trailingOnly = TRUE)
out_csv <- if (length(args) >= 1) args[1] else "step4_results.csv"
quick   <- length(args) >= 2 && args[2] == "quick"
set.seed(20260912)

m       <- if (quick) 120 else 300
n_grid  <- if (quick) c(40) else c(20, 50, 100)
c_depth <- if (quick) c(1, 4) else c(1, 4)
disp    <- if (quick) c("poisson", "nb_pure") else c("poisson", "nb", "nb_pure")
n_rep   <- if (quick) 2 else 8
n_cores <- max(1L, detectCores()); alpha <- 0.05
d0_prior <- 10          # prior weight (pseudo-observations) for dispersion shrinkage
tau_ridge <- 5          # sd of the Gaussian ridge on the zero-part coefficients

# ---------------------------------------------------------------- generator ---
simulate_cell <- function(n, cdep, kind) {
  x <- rep(c(0L, 1L), each = n / 2)
  N <- round(exp(rnorm(n, log(1e4), 0.5)) * ifelse(x == 1, cdep, 1))
  a  <- rnorm(m, 0.5, 1.0); mu <- rnorm(m, log(2e-4), 1.5)
  sg <- if (kind == "nb_pure") rep(0, m) else runif(m, 0.5, 1.5)
  typ <- sample(c("null", "prev_only", "abund_only", "both"), m, TRUE, prob = c(0.80, 0.07, 0.07, 0.06))
  b    <- ifelse(typ %in% c("prev_only", "both"),  sample(c(-1.5, 1.5), m, TRUE), 0)
  beta <- ifelse(typ %in% c("abund_only", "both"), sample(c(-1, 1), m, TRUE), 0)
  X <- matrix(0L, n, m)
  for (j in seq_len(m)) {
    S <- rbinom(n, 1, plogis(a[j] + b[j] * x))
    lam <- exp(mu[j] + beta[j] * x + rnorm(n, 0, sg[j])); mu_c <- N * lam
    cnt <- if (kind == "poisson") rpois(n, mu_c) else rnbinom(n, size = 2, mu = mu_c)
    X[, j] <- S * cnt
  }
  list(X = X, x = x, N = N, typ = typ)
}

safe <- function(expr) tryCatch(expr, error = function(e) NA_real_)
glm_p_x <- function(D, df) { cf <- summary(suppressWarnings(glm(D ~ ., data = df, family = binomial)))$coefficients
  if ("x" %in% rownames(cf)) cf["x", 4] else NA_real_ }

# --------------------------------------------- stage 1: per-taxon ZTNB theta ---
ztnb_theta <- function(Xj, x, N) {
  det <- Xj > 0; if (sum(det) < 6 || length(unique(x[det])) < 2) return(c(NA, NA))
  xd <- Xj[det]; lN <- log(N[det]); xx <- x[det]
  nll <- function(par) { mu_i <- exp(lN + par[1] + par[2] * xx); th <- exp(par[3])
    lp0 <- th * (log(th) - log(th + mu_i))
    -sum(dnbinom(xd, size = th, mu = mu_i, log = TRUE) - log1p(-pmin(exp(lp0), 1 - 1e-12))) }
  f <- nlminb(c(mean(log(xd / N[det])), 0, 0), nll, lower = c(-25, -10, -5), upper = c(25, 10, 8))
  c(log_theta = f$par[3], log_mean = mean(log(xd / N[det])))
}

# ------------------------------------------ stage 2: joint fit, theta fixed ---
jzi_fit <- function(Xj, x, N, theta) {
  D <- Xj > 0; lN <- log(N); n <- length(Xj)
  nll <- function(par, free) {          # par = (a, b, c, beta); free = logical mask
    p <- numeric(4); p[free] <- par
    psi <- plogis(p[1] + p[2] * x); mu <- exp(lN + p[3] + p[4] * x)
    l0 <- dnbinom(0, size = theta, mu = mu)
    P0 <- (1 - psi) + psi * l0
    ll <- sum(log(pmax(P0[!D], 1e-300))) +
          sum(log(pmax(psi[D], 1e-300)) + dnbinom(Xj[D], size = theta, mu = mu[D], log = TRUE))
    -ll + (p[1]^2 + p[2]^2) / (2 * tau_ridge^2)
  }
  fit <- function(free, start) {
    lo <- c(-10, -10, -25, -10)[free]; hi <- -lo
    best <- nlminb(start[free], nll, free = free, lower = lo, upper = hi)
    for (s in list(c(0, 1, 0, 0), c(0, -1, 0, 0), c(0, 0, 0, 1), c(0, 0, 0, -1))) {
      st <- start; st <- st + s; f2 <- nlminb(st[free], nll, free = free, lower = lo, upper = hi)
      if (is.finite(f2$objective) && f2$objective < best$objective) best <- f2 }
    best
  }
  det <- D; c0 <- if (any(det)) mean(log(Xj[det] / N[det])) else -10
  st_full <- c(qlogis(min(max(mean(D), 0.05), 0.95)), 0, c0, 0)
  full <- fit(c(TRUE, TRUE, TRUE, TRUE), st_full)
  pf <- numeric(4); pf[] <- full$par
  null_b    <- fit(c(TRUE, FALSE, TRUE, TRUE), pf)     # b = 0
  null_beta <- fit(c(TRUE, TRUE, TRUE, FALSE), pf)     # beta = 0
  lrt <- function(f0) { s <- 2 * (f0$objective - full$objective); if (is.finite(s) && s < 0 && s > -1e-6) s <- 0
    if (!is.finite(s) || s < 0) NA_real_ else pchisq(s, 1, lower.tail = FALSE) }
  c(p_prev = lrt(null_b), p_abund = lrt(null_beta), b = pf[2], beta = pf[4])
}

# ------------------------------------------------------------- run one rep ---
run_rep <- function(n, cdep, kind, rep_id) {
  sim <- simulate_cell(n, cdep, kind); X <- sim$X; x <- sim$x; N <- sim$N; logN <- log(N)
  # stage 1: thetas + trend + shrinkage
  th <- t(sapply(seq_len(m), function(j) ztnb_theta(X[, j], x, N)))
  ok_th <- is.finite(th[, 1]) & is.finite(th[, 2])
  n_det <- colSums(X > 0)
  trend <- if (sum(ok_th) >= 10) { lo <- suppressWarnings(loess(th[ok_th, 1] ~ th[ok_th, 2], span = 0.75, degree = 1))
             pr <- predict(lo, newdata = th[, 2]); pr[!is.finite(pr)] <- median(th[ok_th, 1]); pr } else rep(median(th[ok_th, 1], na.rm = TRUE), m)
  w <- n_det / (n_det + d0_prior)
  log_theta_star <- ifelse(ok_th, w * th[, 1] + (1 - w) * trend, trend)
  log_theta_star[!is.finite(log_theta_star)] <- 0
  theta_star <- exp(pmin(pmax(log_theta_star, -4), 6))

  res <- lapply(seq_len(m), function(j) {
    Xj <- X[, j]; D <- as.integer(Xj > 0)
    if (sum(D) < 3 || sum(D) > n - 3) return(NULL)
    out <- c()
    out["naive"]    <- safe(glm_p_x(D, data.frame(x = x)))
    out["depthcov"] <- safe(glm_p_x(D, data.frame(x = x, logN = logN)))
    jz <- tryCatch(jzi_fit(Xj, x, N, theta_star[j]), error = function(e) c(p_prev = NA, p_abund = NA, b = NA, beta = NA))
    out["jzi_prev"]  <- jz["p_prev"]; out["jzi_abund"] <- jz["p_abund"]
    out["zinb"] <- if (n >= 50) safe({ zf <- suppressWarnings(zeroinfl(Xj ~ x + offset(logN) | x, dist = "negbin"))
      summary(zf)$coefficients$zero["x", 4] }) else NA_real_
    det <- D == 1
    out["abund_lm_depth"] <- if (length(unique(x[det])) == 2 && sum(det) >= 6)
      safe(summary(lm(log(Xj[det] / N[det]) ~ x[det] + logN[det]))$coefficients[2, 4]) else NA_real_
    out
  })
  keep <- !vapply(res, is.null, logical(1)); P <- do.call(rbind, res[keep]); typ <- sim$typ[keep]
  prev_true <- typ %in% c("prev_only", "both"); abund_true <- typ %in% c("abund_only", "both")
  rows <- lapply(colnames(P), function(meth) {
    pv <- P[, meth]; ok <- is.finite(pv)
    truth <- if (meth %in% c("jzi_abund", "abund_lm_depth")) abund_true else prev_true
    q <- p.adjust(pv[ok], "BH"); rej <- q <= alpha
    tp <- sum(rej & truth[ok]); fp <- sum(rej & !truth[ok]); rr <- tapply(pv[ok] < 0.05, typ[ok], mean)
    data.frame(n = n, c_depth = cdep, disp = kind, rep = rep_id, method = meth,
               n_taxa = sum(ok), n_fail = sum(!ok),
               fdr = if (tp + fp > 0) fp / (tp + fp) else 0, tpr = if (sum(truth[ok]) > 0) tp / sum(truth[ok]) else NA,
               rej_null = rr["null"], rej_prev_only = rr["prev_only"], rej_abund_only = rr["abund_only"], rej_both = rr["both"],
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

grid <- expand.grid(n = n_grid, c_depth = c_depth, disp = disp, rep = seq_len(n_rep), stringsAsFactors = FALSE)
cat(sprintf("m=%d, %d tasks, %d cores\n", m, nrow(grid), n_cores))
t0 <- Sys.time()
res <- mclapply(seq_len(nrow(grid)), function(i) { g <- grid[i, ]; set.seed(4e6 + i)
  tryCatch(run_rep(g$n, g$c_depth, g$disp, g$rep), error = function(e) { message("task ", i, ": ", conditionMessage(e)); NULL }) },
  mc.cores = n_cores)
res <- do.call(rbind, res); write.csv(res, out_csv, row.names = FALSE)
cat(sprintf("done in %.1f min -> %s (%d rows)\n", as.numeric(Sys.time() - t0, units = "mins"), out_csv, nrow(res)))
