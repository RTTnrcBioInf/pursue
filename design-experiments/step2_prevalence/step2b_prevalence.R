#!/usr/bin/env Rscript
# =============================================================================
# Step 2 -- which prevalence arm?
#
# Generative model (per taxon j, sample i):
#   structural presence  S_ij ~ Bern(psi_ij),  logit psi_ij = a_j + b_j x_i
#   abundance if present lambda_ij = exp(mu_j + beta_j x_i + eps_ij)
#   counts               X_ij = S_ij * Poisson(N_i * lambda_ij)   (or NB)
#   depth                N_i ~ LogNormal, multiplied by c_depth in group x=1
# Observed detection D_ij = 1[X_ij > 0] conflates S (structural) with sampling
# non-detection, which depends on N_i and lambda_ij.
#
# Truth types: null | prev_only (b!=0) | abund_only (beta!=0) | both.
#
# Prevalence-arm candidates (all test "b_j = 0"):
#   naive     glm(D ~ x)
#   depthcov  glm(D ~ x + log N)
#   er_lm     PURSUE expected-rarefied presence P^ER (rarefy to d = min N), lm(P^ER ~ x)
#   occ_depth occupancy likelihood P(D)=psi_i * p_i, logit psi = a + b x, with
#             p_i = E_eps[1 - exp(-N_i lam_j e^eps)] plugged in from the abundance arm
#             WITHOUT x (depth-only detection). LRT on b.
#   occ_full  same, lam_ij from the abundance arm INCLUDING x -> absorbs the
#             abundance shadow.
#   occ_cal   occ_full plus a free detection-scale theta: p_i = 1-exp(-theta N lam), so
#             truncation bias in lam_hat is absorbed by theta. LRT on b.
#   zip       pscl::zeroinfl(X ~ x + offset(log N) | x): zero-part coefficient on x
#             (the "proper" zero-inflated reference; slow; known to be fragile)
# Abundance-arm cross-checks (test beta_j = 0 on detected cells):
#   abund_lm    lm(log r ~ x) on detected cells  (MaAsLin 3 / PURSUE style)
#   abund_trunc truncation-corrected lognormal: each detected cell's likelihood is
#               divided by P(detect | mu_i, sigma, N_i), so depth-dependent detection
#               of low-abundance instances no longer biases beta. LRT on beta.
#   occ_trunc   occupancy arm with lam_ij from abund_trunc
#   abund_ztnb  zero-truncated negative binomial on the COUNTS of detected cells,
#               log mu_i = log N_i + mu0 + beta x_i, closed-form likelihood (count part
#               of a hurdle-NB). Handles count discreteness (X=1 cells carry log(1/N),
#               i.e. depth) and overdispersion. LRT on beta.
#   occ_ztnb    occupancy arm with detection p_i = 1 - (theta/(theta+mu_i))^theta from
#               abund_ztnb -- the coherent two-stage detection-aware hurdle.
#   lm_pseudo   lm(log((X+0.5)/N) ~ x) on ALL cells (the LinDA / limma-on-log-TSS
#               baseline): zeros enter as log(0.5/N).
# =============================================================================
suppressMessages({ library(parallel); library(pscl) })

args <- commandArgs(trailingOnly = TRUE)
out_csv <- if (length(args) >= 1) args[1] else "step2b_results.csv"
quick   <- length(args) >= 2 && args[2] == "quick"
set.seed(20260910)

# probabilists' Gauss-Hermite nodes/weights via Golub-Welsch (integrate f(eps) phi(eps) d eps)
gh_quad <- function(K) {
  J <- matrix(0, K, K); od <- sqrt(seq_len(K - 1)); J[cbind(1:(K-1), 2:K)] <- od; J[cbind(2:K, 1:(K-1))] <- od
  e <- eigen(J, symmetric = TRUE); list(x = e$values, w = e$vectors[1, ]^2)
}
GH <- gh_quad(9); gh <- GH$x; gw <- GH$w

m       <- if (quick) 120 else 300
n_grid  <- if (quick) c(40) else c(20, 50, 100)
c_depth <- c(1, 2, 4)
disp    <- c("poisson", "nb")
n_rep   <- if (quick) 2 else 8
run_zip <- TRUE
n_cores <- max(1L, detectCores())
alpha   <- 0.05

simulate_cell <- function(n, cdep, kind) {
  x <- rep(c(0L, 1L), each = n / 2)
  N <- round(exp(rnorm(n, log(1e4), 0.5)) * ifelse(x == 1, cdep, 1))
  a  <- rnorm(m, 0.5, 1.0)                      # baseline logit presence
  mu <- rnorm(m, log(2e-4), 1.5)                # log rel. abundance when present (median expected count ~2 at N=1e4)
  sg <- runif(m, 0.5, 1.5)
  typ <- sample(c("null", "prev_only", "abund_only", "both"), m, replace = TRUE,
                prob = c(0.80, 0.07, 0.07, 0.06))
  b    <- ifelse(typ %in% c("prev_only", "both"),  sample(c(-1.5, 1.5), m, TRUE), 0)
  beta <- ifelse(typ %in% c("abund_only", "both"), sample(c(-1, 1), m, TRUE), 0)
  X <- matrix(0L, n, m)
  for (j in seq_len(m)) {
    S   <- rbinom(n, 1, plogis(a[j] + b[j] * x))
    lam <- exp(mu[j] + beta[j] * x + rnorm(n, 0, sg[j]))
    mu_c <- N * lam
    cnt <- if (kind == "poisson") rpois(n, mu_c) else rnbinom(n, size = 2, mu = mu_c)
    X[, j] <- S * cnt
  }
  list(X = X, x = x, N = N, typ = typ, b = b, beta = beta)
}

# expected-rarefied presence (PURSUE closed form), rarefy to d = min N
er_presence <- function(Xj, N, d) {
  out <- numeric(length(Xj))
  pos <- Xj > 0
  out[pos] <- 1 - exp(lchoose(N[pos] - Xj[pos], d) - lchoose(N[pos], d))
  pmin(pmax(out, 1e-12), 1 - 1e-12)
}

safe_p <- function(expr) tryCatch(expr, error = function(e) NA_real_, warning = function(w) suppressWarnings(expr))

glm_p_x <- function(D, design_df, offset = NULL) {
  fit <- suppressWarnings(glm(D ~ ., data = design_df, family = binomial, offset = offset))
  cf <- summary(fit)$coefficients
  if (!"x" %in% rownames(cf)) return(NA_real_)
  cf["x", "Pr(>|z|)"]
}

test_taxon <- function(Xj, x, N, d) {
  D <- as.integer(Xj > 0); n <- length(Xj)
  if (sum(D) < 3 || sum(D) > n - 3) return(NULL)
  logN <- log(N); lNc <- logN - mean(logN)
  out <- c()
  out["naive"]    <- safe_p(glm_p_x(D, data.frame(x = x)))
  out["depthcov"] <- safe_p(glm_p_x(D, data.frame(x = x, logN = logN)))
  out["er_lm"]    <- safe_p(summary(lm(er_presence(Xj, N, d) ~ x))$coefficients["x", 4])

  det <- D == 1; r <- Xj[det] / N[det]; xd <- x[det]; N_det <- N[det]; lN_det <- log(N_det)
  mu_nox <- mean(log(r)); sg_nox <- if (sum(det) > 2) sd(log(r)) else 1
  enough <- length(unique(xd)) == 2 && sum(det) >= 6

  # --- abundance arms ---
  out["abund_lm"] <- out["abund_lm_depth"] <- out["abund_ztnb"] <- out["abund_ztnb_t"] <- NA_real_
  mu_lm <- rep(mu_nox, n); sg_lm <- sg_nox; mu_lmd <- mu_lm; sg_lmd <- sg_nox
  if (enough) {
    af <- lm(log(r) ~ xd); out["abund_lm"] <- summary(af)$coefficients["xd", 4]
    mu_lm <- coef(af)[1] + coef(af)[2] * x; sg_lm <- summary(af)$sigma
    ad <- lm(log(r) ~ xd + lN_det); cf <- summary(ad)$coefficients
    out["abund_lm_depth"] <- cf["xd", 4]
    mu_lmd <- coef(ad)[1] + coef(ad)[2] * x + coef(ad)[3] * logN; sg_lmd <- summary(ad)$sigma
  }
  # detection prob given presence under Poisson sampling, integrating lognormal spread
  det_prob <- function(mu_vec, sg, theta = 1) {
    p <- numeric(n)
    for (k in seq_along(gh)) p <- p + gw[k] * (1 - exp(-theta * N * exp(mu_vec + sg * gh[k])))
    pmin(pmax(p, 1e-6), 1 - 1e-6)
  }
  # zero-truncated NB (closed form) + NB detection probability; LRT and Wald-t versions
  x_det <- Xj[det]; p_ztnb <- NULL
  ztnb_nll <- function(par, with_x) {
    mu0 <- par[1]; bb <- if (with_x) par[2] else 0; th <- exp(par[length(par)])
    mu_i <- exp(lN_det + mu0 + bb * xd); lp0 <- th * (log(th) - log(th + mu_i))
    -sum(dnbinom(x_det, size = th, mu = mu_i, log = TRUE) - log1p(-pmin(exp(lp0), 1 - 1e-12)))
  }
  if (enough) {
    q0 <- nlminb(c(mu_nox, 0), ztnb_nll, with_x = FALSE, lower = c(-25, -5), upper = c(25, 8))
    q1 <- nlminb(c(q0$par[1], 0, q0$par[2]), ztnb_nll, with_x = TRUE, lower = c(-25, -10, -5), upper = c(25, 10, 8))
    for (b0 in c(-1, 1)) { q2 <- nlminb(c(q0$par[1], b0, q0$par[2]), ztnb_nll, with_x = TRUE, lower = c(-25, -10, -5), upper = c(25, 10, 8))
      if (is.finite(q2$objective) && q2$objective < q1$objective) q1 <- q2 }
    st <- 2 * (q0$objective - q1$objective); if (is.finite(st) && st < 0 && st > -1e-6) st <- 0
    out["abund_ztnb"] <- if (is.finite(st) && st >= 0) pchisq(st, 1, lower.tail = FALSE) else NA_real_
    H <- tryCatch(optimHess(q1$par, ztnb_nll, with_x = TRUE), error = function(e) NULL)
    if (!is.null(H)) { V <- tryCatch(solve(H), error = function(e) NULL)
      if (!is.null(V) && is.finite(V[2, 2]) && V[2, 2] > 0)
        out["abund_ztnb_t"] <- 2 * pt(-abs(q1$par[2] / sqrt(V[2, 2])), df = sum(det) - 3) }
    th1 <- exp(q1$par[3]); mu_all <- exp(logN + q1$par[1] + q1$par[2] * x)
    p_ztnb <- pmin(pmax(1 - exp(th1 * (log(th1) - log(th1 + mu_all))), 1e-6), 1 - 1e-6)
  }

  # --- occupancy arm: P(D=1) = psi_i * p_i ; logit psi = a + b x (+ g lNc) ; LRT on b ---
  occ_lrt_p <- function(pd, with_lN = FALSE) {
    nll <- function(par, with_x) {
      a <- par[1]; b <- if (with_x) par[2] else 0; g <- if (with_lN) par[length(par)] else 0
      P <- pmin(pmax(plogis(a + b * x + g * lNc) * pd, 1e-9), 1 - 1e-9)
      -sum(D * log(P) + (1 - D) * log(1 - P))
    }
    k0 <- 1 + with_lN; lo0 <- rep(-8, k0); hi0 <- rep(8, k0)
    f0 <- nlminb(rep(0, k0), nll, with_x = FALSE, lower = lo0, upper = hi0)
    st1 <- if (with_lN) c(f0$par[1], 0, f0$par[2]) else c(f0$par[1], 0)
    lo1 <- rep(-8, k0 + 1); hi1 <- rep(8, k0 + 1)
    f1 <- nlminb(st1, nll, with_x = TRUE, lower = lo1, upper = hi1)
    for (b0 in c(-1, 1)) { st2 <- st1; st2[2] <- b0
      f2 <- nlminb(st2, nll, with_x = TRUE, lower = lo1, upper = hi1)
      if (is.finite(f2$objective) && f2$objective < f1$objective) f1 <- f2 }
    stat <- 2 * (f0$objective - f1$objective); if (is.finite(stat) && stat < 0 && stat > -1e-6) stat <- 0
    if (!is.finite(stat) || stat < 0) return(NA_real_)
    pchisq(stat, 1, lower.tail = FALSE)
  }
  out["occ_depth"]   <- safe_p(occ_lrt_p(det_prob(rep(mu_nox, n), sg_nox)))
  out["occ_full"]    <- safe_p(occ_lrt_p(det_prob(mu_lm, sg_lm)))
  out["occ_lmd"]     <- safe_p(occ_lrt_p(det_prob(mu_lmd, sg_lmd)))
  out["occ_lmd_lN"]  <- safe_p(occ_lrt_p(det_prob(mu_lmd, sg_lmd), with_lN = TRUE))
  out["occ_ztnb"]    <- if (is.null(p_ztnb)) NA_real_ else safe_p(occ_lrt_p(p_ztnb))
  out["occ_ztnb_lN"] <- if (is.null(p_ztnb)) NA_real_ else safe_p(occ_lrt_p(p_ztnb, with_lN = TRUE))

  # --- joint references ---
  out["zip"] <- safe_p({ zf <- suppressWarnings(zeroinfl(Xj ~ x + offset(logN) | x, dist = "poisson"))
    summary(zf)$coefficients$zero["x", "Pr(>|z|)"] })
  out["zinb"] <- if (n >= 50) safe_p({ zf <- suppressWarnings(zeroinfl(Xj ~ x + offset(logN) | x, dist = "negbin"))
    summary(zf)$coefficients$zero["x", "Pr(>|z|)"] }) else NA_real_
  out
}

run_rep <- function(n, cdep, kind, rep_id) {
  sim <- simulate_cell(n, cdep, kind)
  d <- min(sim$N)
  res <- lapply(seq_len(m), function(j) test_taxon(sim$X[, j], sim$x, sim$N, d))
  keep <- !vapply(res, is.null, logical(1))
  P <- do.call(rbind, res[keep]); typ <- sim$typ[keep]
  prev_true  <- typ %in% c("prev_only", "both")
  abund_true <- typ %in% c("abund_only", "both")
  rows <- list()
  for (meth in colnames(P)) {
    pv <- P[, meth]; ok <- is.finite(pv)
    truth <- if (meth %in% c("abund_lm", "abund_lm_depth", "abund_ztnb", "abund_ztnb_t")) abund_true else prev_true
    q <- p.adjust(pv[ok], "BH"); rej <- q <= alpha
    tp <- sum(rej & truth[ok]); fp <- sum(rej & !truth[ok])
    # per-truth-type rejection rate at p < 0.05
    rr <- tapply(pv[ok] < 0.05, typ[ok], mean)
    rows[[meth]] <- data.frame(
      n = n, c_depth = cdep, disp = kind, rep = rep_id, method = meth,
      n_taxa = sum(ok), n_fail = sum(!ok),
      fdr = if (tp + fp > 0) fp / (tp + fp) else 0,
      tpr = if (sum(truth[ok]) > 0) tp / sum(truth[ok]) else NA,
      rej_null = rr["null"], rej_prev_only = rr["prev_only"],
      rej_abund_only = rr["abund_only"], rej_both = rr["both"],
      stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

grid <- expand.grid(n = n_grid, c_depth = c_depth, disp = disp, rep = seq_len(n_rep),
                    stringsAsFactors = FALSE)
cat(sprintf("m=%d taxa, %d tasks, %d cores\n", m, nrow(grid), n_cores))
t0 <- Sys.time()
res <- mclapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]; set.seed(2e6 + i)
  tryCatch(run_rep(g$n, g$c_depth, g$disp, g$rep),
           error = function(e) { message("task ", i, " failed: ", conditionMessage(e)); NULL })
}, mc.cores = n_cores)
res <- do.call(rbind, res)
write.csv(res, out_csv, row.names = FALSE)
cat(sprintf("done in %.1f min -> %s (%d rows)\n", as.numeric(Sys.time() - t0, units = "mins"), out_csv, nrow(res)))
