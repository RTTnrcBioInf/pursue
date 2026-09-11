#!/usr/bin/env Rscript
# =============================================================================
# Step 3 -- compositional correction: where does the reference go?
#
# Generative model: absolute abundances A_ij = exp(mu_j + beta_j x_i + eps_ij),
# relative r_i = A_i / sum_j A_ij, counts X_i ~ Multinomial(N_i, r_i).
# Truth is beta_j on the ABSOLUTE scale. The relative-scale effect of a null
# taxon is -shift_i, where shift = log(sum_j A_ij | x=1) - log(sum_j A_ij | x=0).
#
# Every method starts from the same per-taxon fit: lm(log((X+0.5)/N) ~ x) (the
# field baseline; step 2 handles zeros/depth separately). They differ only in
# how the effect estimates are centred:
#
#   none      raw relative-abundance effect (no correction)
#   clr       per-sample CLR transform, then lm  (ALDEx2-style centering)
#   median    LinDA: beta - median(beta), SE unchanged
#   mode      LinDA alternative: beta - mode(beta) (density mode), SE unchanged
#   huber     smoothed-median (Huber M-estimate, precision-weighted) centre with the
#             centre's uncertainty propagated into every SE  [the radEmu idea, linear]
#   enull     empirical-null centre (Efron): 2-component mixture on beta_j with known
#             SE_j -- null N(delta, SE_j^2), alt N(eta, omega^2 + SE_j^2) -- by EM;
#             centre = delta, SE(delta) from the null component's precision
#   dacomp    outcome-blind reference: 20% most stable taxa by median pairwise
#             log-ratio SD (DACOMP / PURSUE score); effect = coef of log(X_j / sum_ref X)
#   adapt     reference = half of taxa nearest median(beta) (ADAPT), same log-ratio refit
#   oracle    beta - true shift (ceiling)
#
# Scenarios: DA fraction pi1, direction balance, community evenness (dominant
# taxon = low-complexity, vaginal-like), and a bloom (one taxon x20).
# =============================================================================
suppressMessages(library(parallel))
args <- commandArgs(trailingOnly = TRUE)
out_csv <- if (length(args) >= 1) args[1] else "step3_results.csv"
quick   <- length(args) >= 2 && args[2] == "quick"
set.seed(20260911)

m <- if (quick) 150 else 300; n <- 100; n_rep <- if (quick) 2 else 30
pi1_grid <- c(0.05, 0.20, 0.40); balance_grid <- c("balanced", "all_up")
comm_grid <- c("even", "dominant"); bloom_grid <- c(FALSE, TRUE)
effect_grid <- c("fixed", "graded")     # fixed: |beta| = 1 ; graded: |beta| ~ U(0.2, 1.5) (stress for mixture centring)
delta <- 1.0; alpha <- 0.05
n_cores <- max(1L, detectCores())

simulate_cell <- function(pi1, balance, comm, bloom, effect) {
  x <- rep(c(0L, 1L), each = n / 2)
  N <- round(exp(rnorm(n, log(2e4), 0.4)))
  mu <- rnorm(m, 0, 1.5)
  if (comm == "dominant") mu[1] <- mu[1] + 6            # one taxon ~ 60-90% of the community
  sg <- runif(m, 0.6, 1.2)
  is_da <- rep(FALSE, m); is_da[sample(seq(2, m), round(pi1 * m))] <- TRUE
  sgn <- if (balance == "balanced") sample(c(-1, 1), m, TRUE) else rep(1, m)
  mag <- if (effect == "fixed") rep(delta, m) else runif(m, 0.2, 1.5)
  beta <- ifelse(is_da, mag * sgn, 0)
  bloom_idx <- NA
  if (bloom) { bloom_idx <- 2; is_da[2] <- TRUE; beta[2] <- 3 }   # x20 bloom
  logA <- matrix(mu, n, m, byrow = TRUE) + outer(x, beta) + matrix(rnorm(n * m, 0, rep(sg, each = n)), n, m)
  A <- exp(logA)
  shift <- mean(log(rowSums(A)[x == 1])) - mean(log(rowSums(A)[x == 0]))   # realised community shift
  R <- A / rowSums(A)
  X <- t(sapply(seq_len(n), function(i) rmultinom(1, N[i], R[i, ])))
  list(X = X, x = x, N = N, beta = beta, is_da = is_da, shift = shift, bloom_idx = bloom_idx)
}

fit_lm_all <- function(Y, x) {          # Y: n x m response; returns beta, se, df per column
  Xd <- cbind(1, x); qrX <- qr(Xd); df <- n - 2
  B <- qr.coef(qrX, Y); E <- qr.resid(qrX, Y)
  s2 <- colSums(E^2) / df; vxx <- chol2inv(qr.R(qrX))[2, 2]
  list(beta = B[2, ], se = sqrt(s2 * vxx), df = df)
}

huber_center <- function(b, se, k = 1.345) {
  w0 <- 1 / se^2; mu <- median(b)
  for (it in 1:50) {
    r <- (b - mu) / se; w <- w0 * pmin(1, k / pmax(abs(r), 1e-12))
    mu_new <- sum(w * b) / sum(w); if (abs(mu_new - mu) < 1e-8) break; mu <- mu_new
  }
  psi <- pmax(pmin((b - mu) / se, k), -k)                     # Huber psi
  # sandwich variance of the M-estimate: sum(psi^2/se^2) / (sum(dpsi/se^2))^2
  dpsi <- as.numeric(abs((b - mu) / se) <= k)
  v <- sum((psi / se)^2) / (sum(dpsi / se^2))^2
  list(center = mu, se = sqrt(v))
}

dens_mode <- function(b) { d <- density(b, bw = "SJ"); d$x[which.max(d$y)] }

enull_center <- function(b, se, iters = 200) {
  delta <- dens_mode(b); eta <- mean(b); omega2 <- var(b); pi0 <- 0.8
  for (it in seq_len(iters)) {
    f0 <- dnorm(b, delta, se); f1 <- dnorm(b, eta, sqrt(omega2 + se^2))
    w <- pi0 * f0 / (pi0 * f0 + (1 - pi0) * f1 + 1e-300)
    pi0_new <- mean(w)
    delta_new <- sum(w * b / se^2) / sum(w / se^2)
    v1 <- omega2 + se^2
    eta_new <- sum((1 - w) * b / v1) / sum((1 - w) / v1)
    omega2_new <- max(sum((1 - w) * ((b - eta_new)^2 - se^2) / v1) / sum((1 - w) / v1), 1e-4)
    conv <- abs(delta_new - delta) < 1e-7
    delta <- delta_new; eta <- eta_new; omega2 <- omega2_new; pi0 <- min(max(pi0_new, 0.05), 0.99)
    if (conv) break
  }
  list(center = delta, se = sqrt(1 / sum(w / se^2)), pi0 = pi0)
}

dacomp_reference <- function(X, frac = 0.20) {
  L <- log(X + 0.5); L <- L - rowMeans(L)                      # any per-sample constant cancels in ratios
  G <- crossprod(L) ; d <- diag(G); D2 <- outer(d, d, "+") - 2 * G; D2[D2 < 0] <- 0
  S <- sqrt(D2 / (nrow(X) - 1)); diag(S) <- NA
  score <- apply(S, 1, median, na.rm = TRUE)
  which(rank(score) <= max(5, round(frac * ncol(X))))
}

ratio_refit <- function(X, x, ref) {
  denom <- rowSums(X[, ref, drop = FALSE]) + 0.5
  Y <- log((X + 0.5) / denom)
  f <- fit_lm_all(Y, x); f$beta[ref] <- NA; f$se[ref] <- NA; f
}

run_rep <- function(pi1, balance, comm, bloom, effect, rep_id) {
  sim <- simulate_cell(pi1, balance, comm, bloom, effect)
  X <- sim$X; x <- sim$x; N <- sim$N
  base <- fit_lm_all(log((X + 0.5) / N), x)
  b <- base$beta; se <- base$se; df <- base$df
  est <- list(); ses <- list()
  est$none   <- b;                     ses$none   <- se
  clr <- log(X + 0.5); clr <- clr - rowMeans(clr); fc <- fit_lm_all(clr, x)
  est$clr    <- fc$beta;               ses$clr    <- fc$se
  est$median <- b - median(b);         ses$median <- se
  est$mode   <- b - dens_mode(b);      ses$mode   <- se
  hc <- huber_center(b, se)
  est$huber  <- b - hc$center;         ses$huber  <- sqrt(se^2 + hc$se^2)
  en <- enull_center(b, se)
  est$enull  <- b - en$center;         ses$enull  <- sqrt(se^2 + en$se^2)
  ref_d <- dacomp_reference(X); fd <- ratio_refit(X, x, ref_d)
  est$dacomp <- fd$beta;               ses$dacomp <- fd$se
  ref_a <- order(abs(b - median(b)))[seq_len(round(m / 2))]; fa <- ratio_refit(X, x, ref_a)
  est$adapt  <- fa$beta;               ses$adapt  <- fa$se
  est$oracle <- b + sim$shift;         ses$oracle <- se      # relative = absolute - shift  => absolute = relative + shift

  truth <- sim$beta; is_da <- sim$is_da
  rows <- lapply(names(est), function(meth) {
    e <- est[[meth]]; s <- ses[[meth]]; ok <- is.finite(e) & is.finite(s)
    tstat <- e / s; p <- 2 * pt(-abs(tstat), df); q <- p.adjust(p[ok], "BH"); rej <- q <= alpha
    tp <- sum(rej & is_da[ok]); fp <- sum(rej & !is_da[ok])
    nul <- ok & !is_da; ci_cov <- mean(abs(e[ok] - truth[ok]) <= qt(0.975, df) * s[ok])
    data.frame(pi1 = pi1, balance = balance, comm = comm, bloom = bloom, effect = effect, rep = rep_id, method = meth,
               n_eval = sum(ok), shift = sim$shift,
               bias_null = mean(e[nul]), rmse_all = sqrt(mean((e[ok] - truth[ok])^2)),
               type1 = mean(p[nul] < 0.05), power = mean(p[ok & is_da] < 0.05),
               fdr = if (tp + fp > 0) fp / (tp + fp) else 0, tpr = tp / sum(is_da[ok]),
               ci_cov = ci_cov,
               ref_purity = if (meth == "dacomp") mean(!is_da[ref_d]) else if (meth == "adapt") mean(!is_da[ref_a]) else NA,
               stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

grid <- expand.grid(pi1 = pi1_grid, balance = balance_grid, comm = comm_grid, bloom = bloom_grid,
                    effect = effect_grid, rep = seq_len(n_rep), stringsAsFactors = FALSE)
cat(sprintf("m=%d, n=%d, %d tasks, %d cores\n", m, n, nrow(grid), n_cores))
t0 <- Sys.time()
res <- mclapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]; set.seed(3e6 + i)
  tryCatch(run_rep(g$pi1, g$balance, g$comm, g$bloom, g$effect, g$rep),
           error = function(e) { message("task ", i, " failed: ", conditionMessage(e)); NULL })
}, mc.cores = n_cores)
res <- do.call(rbind, res)
write.csv(res, out_csv, row.names = FALSE)
cat(sprintf("done in %.1f min -> %s (%d rows)\n", as.numeric(Sys.time() - t0, units = "mins"), out_csv, nrow(res)))
