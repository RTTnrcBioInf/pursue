# --- it42: the permutation path for every design (Phase 2, target 1) ------------------------------------------
# eperm_cs1 leaves the exact permutation path, for the asymptotic efull_b3w path, whenever the design is not two
# exchangeable groups: designed depth confounding (R17-R19, rho < 0.51), covariates (R20, R21, B_conf04/07),
# continuous exposure (R22, B_cont) and repeated measures (R23). p13 puts every one of those regimes at 0.49-0.87
# of the best comparator, and the one miscalibrated cell (mid R19, FDR 0.123) is on that path.
#
# The pairwise statistic generalises to all four. With the label-free pair graph E and antisymmetric kernel A_ij,
#   sum over pairs of (x_i - x_j) A_ij = 2 sum_i x_i s_i,  s_i = sum_{j ~ i} A_ij,
# so for ANY exposure x the U-statistic is linear in fixed sample scores s, and the test is a regression of s on x:
# - depth confounding (binary, rho < 0.51): the common-depth thinning makes every score mean-zero under H0 whatever
#   the depths; only the score variances differ by depth, which the studentised (Welch) permutation tolerates
#   (Janssen 1997; Chung & Romano 2013). So: eperm_cs1 itself with rho_min = 0.
# - continuous exposure, no covariates: permuting x is exact; the statistic is the HC0-studentised slope of s on x.
# - covariates Z: Freedman-Lane -- residualise s on (1, Z), permute the residuals, refit s ~ 1 + Z + x and take
#   the HC0-studentised x coefficient (asymptotically valid; the studentised form is the robust choice, Winkler et
#   al. 2014). Every permuted quantity is a matrix product of the residual scores with a permuted design vector,
#   so 2p + 2 products per batch (p = columns of the full design).
# - repeated measures with subject-level exposure: sum the scores within subject (within-subject pairs cancel) and
#   permute subjects -- exact.
# The compositional centre is re-estimated inside every permutation exactly as in it39 (linearised in gamma), and
# the scale is chosen by leave-one-out BH counts with margin 1 as in eperm_cs1. Two-group exchangeable designs run
# eperm_cs1 unchanged. New names only: 60_erd.R is frozen while the conf1 confirmation runs.
.epg_design <- function(meta, formula, tested_term) {
  tl <- attr(stats::terms(formula), "term.labels"); x <- meta[[tested_term]]
  if (is.null(x) || !tested_term %in% tl || anyNA(x)) return(NULL)
  ux <- unique(x); binary <- length(ux) == 2L
  if (binary) { g1 <- x == sort(ux)[2]; xv <- as.numeric(g1) } else if (is.numeric(x)) { g1 <- NULL; xv <- as.numeric(x) } else return(NULL)
  X <- stats::model.matrix(formula, meta); asg <- attr(X, "assign"); zc <- which(asg != 0 & asg != which(tl == tested_term))
  cl <- .cluster_of(meta)
  cl_ok <- is.null(cl) || all(tapply(xv, cl, function(v) length(unique(v))) == 1L)       # exposure constant within subject
  list(x = xv, g1 = g1, binary = binary, Z = if (length(zc)) X[, zc, drop = FALSE] else NULL, cl = cl, cl_ok = cl_ok)
}
.epg_rho <- function(depth, d) {                                   # design thinning factor, as .rho_design; continuous x via r -> d
  if (d$binary) return(.rho_design(depth, d$g1))
  r <- suppressWarnings(stats::cor(log(depth), d$x)); if (!is.finite(r)) return(1)
  dd <- 2 * abs(r) / sqrt(max(1 - r^2, 1e-8)); 1 - 0.5 * min(max((dd - 0.25) / 0.25, 0), 1)
}
.fl_prep <- function(x, Z) {                                       # Freedman-Lane pieces for s ~ 1 + Z + x
  X0 <- if (is.null(Z)) matrix(1, length(x), 1) else cbind(1, Z); X <- cbind(X0, x); q0 <- qr(X0); m0x <- qr.resid(q0, x); w <- m0x / sum(m0x * x)   # w: x row of (X'X)^-1 X'
  list(n = length(x), X = X, q0 = q0, w = w, XtXi = solve(crossprod(X)), W2X = X * w^2, Q = crossprod(X * w))
}
.fl_resid <- function(S, pr) t(qr.resid(pr$q0, t(S)))            # scores residualised on (1, Z), taxa x samples
.fl_t <- function(R, R2, pr, P) {                                  # HC0 t of x for permutations P (n x B indices; identity = observed)
  B <- ncol(P); pm <- function(a) matrix(a[P], nrow(P), B); p <- ncol(pr$X)
  beta <- R %*% pm(pr$w); V <- R2 %*% pm(pr$w^2)
  Ck <- lapply(seq_len(p), function(k) R %*% pm(pr$X[, k])); Ak <- lapply(seq_len(p), function(k) R %*% pm(pr$W2X[, k]))
  cc <- lapply(seq_len(p), function(l) { z <- 0; for (k in seq_len(p)) z <- z + pr$XtXi[l, k] * Ck[[k]]; z })
  for (l in seq_len(p)) V <- V - 2 * cc[[l]] * Ak[[l]]
  for (k in seq_len(p)) for (l in seq_len(p)) V <- V + pr$Q[k, l] * cc[[k]] * cc[[l]]
  beta / sqrt(pmax(V, 1e-300))
}
.epg_perm_fl <- function(S0, pr, combos, B = 4000L, B2 = 40000L, cmin = 10L, seed = 17L, chunk = 1000L, tobs, kap, med, ck = "log") {
  n <- pr$n; m <- nrow(S0[[1]]); R <- lapply(S0, .fl_resid, pr = pr); R2 <- lapply(R, function(r) r^2)
  stat <- function(tm, cmb) { a <- abs(tm[[cmb[1]]]); for (h in cmb[-1]) a <- pmax(a, abs(tm[[h]])); a }
  obs <- matrix(sapply(combos, function(cmb) { a <- abs(tobs[, cmb[1]]); for (h in cmb[-1]) a <- pmax(a, abs(tobs[, h])); a }), m, length(combos))
  cnt <- matrix(0, m, length(combos), dimnames = list(NULL, names(combos))); tot <- rep(0, m)
  run <- function(rows, Bn, sd0) .with_seed(sd0, { done <- 0L; rr <- union(rows, med); at <- match(rows, rr)
    while (done < Bn) { bb <- min(chunk, Bn - done); P <- vapply(seq_len(bb), function(k) sample.int(n), integer(n))
      tt <- function(h, r) { x <- .fl_t(R[[h]][r, , drop = FALSE], R2[[h]][r, , drop = FALSE], pr, P); x[!is.finite(x)] <- 0; x }
      tc <- tt(ck, rr); dl <- .root_shift(tc[match(med, rr), , drop = FALSE], kap[med, ck])
      tm <- lapply(names(S0), function(h) (if (h == ck) tc[at, , drop = FALSE] else tt(h, rows)) + outer(kap[rows, h], dl)); names(tm) <- names(S0)
      for (cc in seq_along(combos)) cnt[rows, cc] <<- cnt[rows, cc] + rowSums(stat(tm, combos[[cc]]) >= obs[rows, cc] - 1e-10)
      tot[rows] <<- tot[rows] + bb; done <- done + bb } })
  run(seq_len(m), B, seed)
  low <- which(apply(cnt, 1, min) < cmin); if (length(low) && B2 > 0) run(low, B2, seed + 1L)
  list(p = (cnt + 1) / (tot + 1), tot = tot)
}
.epg_memo <- new.env()
.epg_fit <- function(counts, meta, formula, tested_term, thr = 1e-3, margin = 1L, K = 60L, B = 4000L, B2 = 40000L, winsor = 0.03,
                     paths = c("rho", "cont", "cov", "clus"), adet = FALSE, nodet_conf = FALSE, pk = character(0)) {
  key <- list(counts, meta, formula, tested_term, thr, margin, K, B, B2, winsor, paths, adet, nodet_conf, pk)
  if (!is.null(.epg_memo$key) && identical(.epg_memo$key, key)) return(.epg_memo$val)
  d <- .epg_design(meta, formula, tested_term); depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts)
  base <- function(rmin) { v <- .eperm_fit(counts, meta, formula, tested_term, thr = thr, margin = margin, K = K, B = B, B2 = B2, winsor = winsor,
                                           rho_min = rmin, recentre = TRUE, sel = "bh"); v$path <- if (isTRUE(v$perm)) "two" else "efull"; v }
  type <- if (is.null(d)) NA else if (!is.null(d$cl)) { if (d$binary && is.null(d$Z) && d$cl_ok) "clus" else NA } else
    if (!is.null(d$Z)) "cov" else if (d$binary) "two" else "cont"
  if (is.na(type) || (type != "two" && !type %in% paths)) { val <- base(0.51)
  } else if (type == "two") { rmin <- if ("rho" %in% paths) 0 else 0.51
    val <- if ((adet || length(pk)) && isTRUE(round(.rho_design(depth, d$g1), 2) >= rmin)) .eps_fit(counts, meta, formula, tested_term, thr = thr, margin = margin, K = K, B = B,
      B2 = B2, winsor = winsor, ps = character(0), adet = adet, rho_min = rmin, paths = paths, pk = pk) else base(rmin)
    if (nodet_conf && isFALSE(val$perm) && isTRUE(val$rho < 0.51)) {         # it46: designed depth confounding -- no detection kernel
      val$p_full <- .sel_bh(cbind(log = val$p_log, sqrt = val$p_sqrt), margin = margin); val$path <- "efull_nodet" }
  } else {
    x <- d$x; n <- ncol(counts); rho <- round(.epg_rho(depth, d), 2)
    sf <- if (rho == 1) .size_factors(counts, depth) else rep(0, n)
    if (winsor > 0 && rho == 1) counts <- .winsorize_counts(counts, depth, winsor)
    rows <- .with_seed(7L, if (nrow(counts) > 300L) sort(sample.int(nrow(counts), 300L)) else seq_len(nrow(counts)))
    ct <- counts[rows, , drop = FALSE]; ct <- ct[rowSums(ct > 0) >= 3, , drop = FALSE]
    E <- .pair_graph(n, K); ok <- rowSums(counts > 0) >= 3
    hs <- c(det = "det", log = "log", sqrt = "sqrt"); combos <- list(lead = c("det", "log"), det = "det", log = "log", sqrt = "sqrt"); alts <- c(if (adet) "det", "log", "sqrt")
    for (h in pk) { hs[[h]] <- h; combos[[h]] <- h; alts <- c(alts, h) }             # it47: extra pairwise kernels (pos)
    if (type == "clus") {                                          # subject-level units: scores summed within subject
      Mem <- stats::model.matrix(~ d$cl - 1); agg <- function(S) S %*% Mem
      g1u <- as.numeric(tapply(x, d$cl, `[`, 1)) == 1; nu <- length(g1u); G1 <- matrix(as.numeric(g1u), nu, 1)
      tstat <- function(SS) { v <- sapply(hs, function(h) { a <- agg(SS[[h]]); a <- a - rowMeans(a); .perm_tstat(a, a^2, G1, sum(g1u), nu - sum(g1u))[, 1] })
        v[!is.finite(v)] <- 0; matrix(v, ncol = length(hs), dimnames = list(NULL, hs)) }
    } else {
      pr <- .fl_prep(x, d$Z); I1 <- matrix(seq_len(n), n, 1)
      tstat <- function(SS) { v <- sapply(hs, function(h) { r <- .fl_resid(SS[[h]], pr); .fl_t(r, r^2, pr, I1)[, 1] })
        v[!is.finite(v)] <- 0; matrix(v, ncol = length(hs), dimnames = list(NULL, hs)) }
    }
    sc <- function(cnt, gm, h) if (h == "pos") .perm_scores_pos(cnt, depth, exp(gm * x + sf), E, rho = rho) else .perm_scores(cnt, depth, exp(gm * x + sf), h, E, rho = rho)
    medt <- function(gm) { SS <- list(det = NULL, log = sc(ct, gm, "log"), sqrt = NULL)
      t <- if (type == "clus") { a <- agg(SS$log); a <- a - rowMeans(a); .perm_tstat(a, a^2, G1, sum(g1u), nu - sum(g1u))[, 1] } else { r <- .fl_resid(SS$log, pr); .fl_t(r, r^2, pr, I1)[, 1] }
      stats::median(t[is.finite(t)]) }
    flo <- medt(-2); fhi <- medt(2)
    gm <- if (is.finite(flo) && is.finite(fhi) && sign(flo) != sign(fhi)) stats::uniroot(medt, c(-2, 2), f.lower = flo, f.upper = fhi, tol = 2e-3)$root else 0
    cok <- counts[ok, , drop = FALSE]
    S <- lapply(hs, function(h) sc(cok, gm, h)); S0 <- lapply(hs, function(h) sc(cok, 0, h))
    t1 <- tstat(S); t0 <- tstat(S0)
    kap <- if (abs(gm) >= 0.02) (t1 - t0) / gm else (tstat(lapply(hs, function(h) sc(cok, 0.05, h))) - t0) / 0.05
    kap[!is.finite(kap)] <- 0
    med <- which(rownames(cok) %in% rownames(ct))
    tobs <- t1 + kap * .root_shift(t1[med, "log", drop = FALSE], kap[med, "log"])
    pm <- if (type == "clus") .perm_multi(lapply(S0, agg), g1u, combos, B = B, B2 = B2, tobs = tobs, kap = kap, med = med) else
      .epg_perm_fl(S0, pr, combos, B = B, B2 = B2, tobs = tobs, kap = kap, med = med)
    P <- matrix(NA_real_, nrow(counts), ncol(pm$p), dimnames = list(rownames(counts), colnames(pm$p))); P[ok, ] <- pm$p
    pl <- P[, "lead"]; Q <- P[, c("lead", alts), drop = FALSE]
    hit <- apply(Q, 2, function(p) { q <- rep(1, length(p)); f <- is.finite(p); q[f] <- stats::p.adjust(p[f], "BH"); q <= 0.05 })
    loo <- sweep(-hit, 2, colSums(hit), "+"); la <- loo[, alts, drop = FALSE]; ia <- max.col(la, ties.method = "last"); altn <- la[cbind(seq_len(nrow(la)), ia)]
    use <- altn >= loo[, "lead"] + margin; pfull <- ifelse(use, Q[cbind(seq_len(nrow(Q)), 1L + ia)], pl)
    xc <- x - mean(x); deg <- 2 * nrow(E) / n; est <- rep(NA_real_, nrow(counts))
    est[ok] <- if (type == "clus") rowSums(S$log[, x == 1, drop = FALSE]) / sum(x[E[, 1]] != x[E[, 2]]) / log(2) else drop(.fl_resid(S$log, pr) %*% pr$w) / deg / log(2)
    val <- list(feature = rownames(counts), p_full = pfull, p_ms = pl, p_det = P[, "det"], p_log = P[, "log"], p_sqrt = P[, "sqrt"],
                est = est, gam = gm, rho = rho, perm = TRUE, fallback = FALSE, path = type, choice = ifelse(use, alts[ia], "lead"))
  }
  .epg_memo$key <- key; .epg_memo$val <- val; val
}
register_candidate("eperm_g", function(counts, meta, formula, tested_term) {
  v <- .epg_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_full, estimate = v$est)
}, notes = "it42: eperm_cs1 with the permutation path for every design -- depth confounding, continuous exposure, covariates (Freedman-Lane), repeated measures (subject permutation)")
# it44: detection alone as a selection option. The lead max(det, log) pays a combination penalty that the selection
# cannot undo, because det alone was never offered: one twinsuk B_prev cell, BH discoveries det 71 / lead 67 / log 24
# (true positives 68 / 64 / 24); logistic presence 77 (74). Offering det lets presence-driven data use it whole.
register_candidate("eperm_gd", function(counts, meta, formula, tested_term) {
  v <- .epg_fit(counts, meta, formula, tested_term, adet = TRUE); data.frame(feature = v$feature, p = v$p_full, estimate = v$est)
}, notes = "it42 + it44: eperm_g with detection alone among the selection options")

# --- srv18 verdicts (2026-10-09): keep the continuous-exposure path (house R22 22.85 -> 25.50, B_cont 10.60 -> 11.75;
# nulls to be confirmed at more reps) and detection-as-option in exchangeable two-group designs (B_prev 33.2 -> 35.3);
# drop the depth-confounded path (power down on msq / house R19, P(any) up: house R19.null 0.125, msq R19.null 0.175,
# mid R19.null 0.70 -- heteroscedastic scores break the extreme tail), the subject-permutation path (R23 5.35 -> 3.25:
# 20 units) and Freedman-Lane for covariates (R21 +0.2, B_conf07 -1.45). eperm_c2 is the combination.
register_candidate("eperm_c2", function(counts, meta, formula, tested_term) {
  v <- .epg_fit(counts, meta, formula, tested_term, paths = "cont", adet = TRUE); data.frame(feature = v$feature, p = v$p_full, estimate = v$est)
}, notes = "it42/44 kept parts: eperm_cs1 + continuous-exposure permutation path + detection as a selection option (exchangeable two-group)")
# it45: Storey's adaptive BH as PURSUE's own FDR procedure. BH controls FDR at pi0 q, so with 10-40% DA a fixed q leaves
# 5-30% of the error budget unused (dev suite FDR 0.007-0.06 at q 0.05). q = pi0 x BH, pi0 = min(1, (1 + #{p > 0.5}) /
# (m / 2)): Storey, Taylor & Siegmund 2004, finite-sample FDR control under independence; the +1 keeps pi0 = 1 under a
# global null, so P(any) there is unchanged. A method's own FDR procedure is part of the method (elementary.R).
.storey_q <- function(p) { f <- is.finite(p); q <- rep(NA_real_, length(p)); m <- sum(f); if (!m) return(q)
  pi0 <- min(1, (1 + sum(p[f] > 0.5)) / (m * 0.5)); q[f] <- pmin(1, pi0 * stats::p.adjust(p[f], "BH")); q }
register_candidate("eperm_c2q", function(counts, meta, formula, tested_term) {
  v <- .epg_fit(counts, meta, formula, tested_term, paths = "cont", adet = TRUE); data.frame(feature = v$feature, p = v$p_full, q = .storey_q(v$p_full), estimate = v$est)
}, notes = "it45: eperm_c2 with Storey's adaptive BH (pi0 from p > 0.5, +1 correction)")

# it46: no detection kernel under designed depth confounding. Local mid R19 nulls (hmp_tongue, 8 cells, efull path):
# BH false discoveries in 8 / 8 cells with det or the lead, 2 / 8 with log or sqrt alone; null tail at p < 0.01: det 3.3x
# nominal, lead 3.0x, log 1.5x, sqrt 1.8x. MIDASim draws presence from a depth-dependent model on top of the sampling at
# that depth, so deep samples carry presences that no rarefaction can thin away -- detection cannot be made
# depth-invariant there, abundance given the sampling can. Cost on house (4 cells): R17 lead 8.75 vs log 7.25, R19 10.75
# vs 10.75. Used only where the design is depth-confounded beyond chance (rho < 0.51, the efull path).
.sel_bh <- function(Q, margin = 1L) {                              # first column the default; LOO BH-count selection (it40)
  hit <- apply(Q, 2, function(p) { q <- rep(1, length(p)); f <- is.finite(p); q[f] <- stats::p.adjust(p[f], "BH"); q <= 0.05 })
  loo <- sweep(-hit, 2, colSums(hit), "+"); ifelse(loo[, 2] >= loo[, 1] + margin, Q[, 2], Q[, 1]) }
register_candidate("eperm_c3", function(counts, meta, formula, tested_term) {
  v <- .epg_fit(counts, meta, formula, tested_term, paths = "cont", adet = TRUE, nodet_conf = TRUE); data.frame(feature = v$feature, p = v$p_full, estimate = v$est)
}, notes = "it46: eperm_c2 without the detection kernel under designed depth confounding (rho < 0.51)")
