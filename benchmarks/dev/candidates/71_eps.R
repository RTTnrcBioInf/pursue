# --- it43: per-sample scales in the exchangeable two-group path (Phase 2, target 2: implants) -------------------
# p13, axis B: logistic presence -- presence ~ group, no depth term at all -- leads 10 of 13 implant rows (B_prev
# 8.15 vs 6.99 TP per cell, B_small 1.20 vs 0.84). Our detection kernel compares each pair at the pair's common
# depth, so a deep sample's detections are thinned to its partner's depth; that keeps depth confounding out, but in
# a randomised design depth is exchangeable and the alignment only throws reads away.
# What validity needs there is narrower: the compositional shift between groups (gamma) must not leak into null
# taxa. Binomial thinning removes it exactly (a thinned multinomial is a multinomial at the thinned depth), so:
# thin only the group whose null taxa are inflated, by exp(-|gamma|), keep each sample at its own depth, and take
#   dps: expected detection at that depth,  lps: expected log(1 + count) at that depth,
# each residualised per taxon on log(thinned depth) across all samples (label-free, removes depth noise). Unlike
# it41's lin / sqp (a multiplicative correction, which matches means but not sampling distributions and failed the
# x4 bloom), thinning matches the whole sampling distribution. Size factors stay in the pairwise kernels only:
# per-sample compositional noise is exchangeable here, so it costs precision, not validity.
# The two scales join the leave-one-out BH selection (margin 1) as options next to log and sqrt; gamma is
# re-estimated in every permutation as in it39. Only where the permutation path is exact for them -- two groups,
# no covariates, depth imbalance at chance level (rho >= 0.51); every other design runs it42 (eperm_g).
.ps_scores <- function(counts, depth, x, gm, h) {
  e <- gm * x; D <- pmax(floor(depth * exp(e - max(e))), 1)
  Y <- as.matrix(counts); storage.mode(Y) <- "double"
  F <- if (h == "dps") .erd_f_var(Y, depth, D) else .erl_f_var(Y, depth, D)
  q <- qr(cbind(1, log(D))); F - t(qr.fitted(q, t(F)))
}
.eps_memo <- new.env()
.eps_fit <- function(counts, meta, formula, tested_term, thr = 1e-3, margin = 1L, K = 60L, B = 4000L, B2 = 40000L, winsor = 0.03,
                     ps = c("dps", "lps"), adet = FALSE, rho_min = 0.51, paths = c("rho", "cont", "cov", "clus"), pk = character(0)) {
  key <- list(counts, meta, formula, tested_term, thr, margin, K, B, B2, winsor, ps, adet, rho_min, paths, pk)
  if (!is.null(.eps_memo$key) && identical(.eps_memo$key, key)) return(.eps_memo$val)
  d <- .epg_design(meta, formula, tested_term); depth <- if (!is.null(meta$depth)) meta$depth else colSums(counts)
  rho <- if (!is.null(d) && d$binary) round(.rho_design(depth, d$g1), 2) else NA
  if (length(ps) && !isTRUE(rho >= 0.51)) ps <- character(0)       # per-sample scales need exchangeable depth
  if (is.null(d) || !d$binary || !is.null(d$Z) || !is.null(d$cl) || !isTRUE(rho >= rho_min) || (!length(ps) && !adet && !length(pk))) {
    val <- .epg_fit(counts, meta, formula, tested_term, thr = thr, margin = margin, K = K, B = B, B2 = B2, winsor = winsor, adet = adet, paths = paths, pk = pk)
  } else {
    g1 <- d$g1; x <- d$x; n <- ncol(counts); raw <- counts
    sf <- if (rho == 1) .size_factors(counts, depth) else rep(0, n)
    if (winsor > 0 && rho == 1) counts <- .winsorize_counts(counts, depth, winsor)
    rows <- .with_seed(7L, if (nrow(counts) > 300L) sort(sample.int(nrow(counts), 300L)) else seq_len(nrow(counts)))
    ct <- counts[rows, , drop = FALSE]; ct <- ct[rowSums(ct > 0) >= 3, , drop = FALSE]
    medf <- function(gm) { u <- .u_stat_cov(ct, depth, g1, NULL, gm, "log", rho = rho, sf = if (rho == 1) sf else NULL); stats::median(u$U / sqrt(u$v), na.rm = TRUE) }
    flo <- medf(-2); fhi <- medf(2)
    gm <- if (is.finite(flo) && is.finite(fhi) && sign(flo) != sign(fhi)) stats::uniroot(medf, c(-2, 2), f.lower = flo, f.upper = fhi, tol = 2e-3)$root else 0
    E <- .pair_graph(n, K); ok <- rowSums(counts > 0) >= 3; cok <- counts[ok, , drop = FALSE]; rok <- raw[ok, , drop = FALSE]
    hs <- c(det = "det", log = "log", sqrt = "sqrt"); combos <- list(lead = c("det", "log"), det = "det", log = "log", sqrt = "sqrt")
    for (h in c(pk, ps)) combos[[h]] <- h
    hs <- c(hs, setNames(pk, pk))
    alts <- setdiff(names(combos), c("lead", if (!adet) "det"))
    scores <- function(gm, sfb) { S <- lapply(hs, function(h) if (h == "pos") .perm_scores_pos(cok, depth, exp(gm * x + sfb), E, rho = rho) else .perm_scores(cok, depth, exp(gm * x + sfb), h, E, rho = rho))
      for (h in ps) S[[h]] <- .ps_scores(rok, depth, x, gm, h); S }                 # per-sample scales: raw counts (no winsorising)
    S <- scores(gm, sf); S0 <- scores(0, sf)
    G1 <- matrix(as.numeric(g1), n, 1); hh <- names(S)
    tst <- function(SS) { v <- sapply(hh, function(h) { a <- SS[[h]] - rowMeans(SS[[h]]); .perm_tstat(a, a^2, G1, sum(g1), n - sum(g1))[, 1] })
      v[!is.finite(v)] <- 0; matrix(v, ncol = length(hh), dimnames = list(NULL, hh)) }
    t1 <- tst(S); t0 <- tst(S0)
    kap <- if (abs(gm) >= 0.02) (t1 - t0) / gm else (tst(scores(0.05, sf)) - t0) / 0.05
    kap[!is.finite(kap)] <- 0
    med <- which(rownames(cok) %in% rownames(ct))
    tobs <- t1 + kap * .root_shift(t1[med, "log", drop = FALSE], kap[med, "log"])
    pm <- .perm_multi(S0, g1, combos, B = B, B2 = B2, tobs = tobs, kap = kap, med = med)
    P <- matrix(NA_real_, nrow(counts), ncol(pm$p), dimnames = list(rownames(counts), colnames(pm$p))); P[ok, ] <- pm$p
    pl <- P[, "lead"]; Q <- P[, c("lead", alts), drop = FALSE]
    hit <- apply(Q, 2, function(p) { q <- rep(1, length(p)); f <- is.finite(p); q[f] <- stats::p.adjust(p[f], "BH"); q <= 0.05 })
    loo <- sweep(-hit, 2, colSums(hit), "+"); la <- loo[, alts, drop = FALSE]; ia <- max.col(la, ties.method = "last"); altn <- la[cbind(seq_len(nrow(la)), ia)]
    use <- altn >= loo[, "lead"] + margin; pfull <- ifelse(use, Q[cbind(seq_len(nrow(Q)), 1L + ia)], pl)
    ncross <- sum(g1[E[, 1]] != g1[E[, 2]]); est <- rep(NA_real_, nrow(counts)); est[ok] <- rowSums(S$log[, g1, drop = FALSE]) / ncross / log(2)
    val <- list(feature = rownames(counts), p_full = pfull, p_ms = pl, est = est, gam = gm, rho = rho, perm = TRUE, fallback = FALSE, path = if (length(ps)) "two_ps" else "two",
                choice = ifelse(use, alts[ia], "lead"), P = P)
  }
  .eps_memo$key <- key; .eps_memo$val <- val; val
}
register_candidate("eperm_ps", function(counts, meta, formula, tested_term) {
  v <- .eps_fit(counts, meta, formula, tested_term); data.frame(feature = v$feature, p = v$p_full, estimate = v$est)
}, notes = "it43: eperm_g + per-sample-thinned detection and log scales (dps, lps) as selection options in exchangeable two-group designs")
