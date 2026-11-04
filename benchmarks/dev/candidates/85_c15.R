# --- it59 / it60: the two gaps chosen after conf2 (2026-11-04) -------------------------------------------------------------
# conf2 read without house: about three-quarters of c14's remaining power deficit is SimulateMSeq against LDM (largest at
# n = 100 / 200 per group), and c14's only failing cells are one-directional designs (mid R10 FDR 0.101, B_allup 0.108),
# where the Vandeputte direction gap most likely also comes from (srv34).
#
# it59 -- the compositional centre on a reference set, weighted by sensitivity. Every path centres on the median taxon: the
# shift d with median over the centring rows of (t + kap d) = 0 (.root_shift; observed statistic and every permutation alike).
# When most changes go one way the changed taxa pull that median, so null taxa inherit a shift -- the larger kap, the larger:
# abundant taxa first. "lad" (below) weights each taxon by its sensitivity and re-solves on the reference rows; "trim" is the
# unweighted ablation. The same rule in the observed statistic and in every permutation, so exactness is as before.
# options(pursue.centre).
#
# it60 -- the plain-proportion kernel with its own centre. probe1: LDM's msq edge is its frequency-scale test under
# permutation; it52's raw linear kernel (lin0) matched it on banked msq cells but failed the x4 bloom (FDR 0.74) because it
# borrowed the log kernel's centre, which on the linear scale leaves every null taxon a residual group difference. it56's
# own-centre rule (a kernel centred on its own median taxon, observed and per permutation) is that missing piece. lin0 is
# offered only in exchangeable two-group designs (the path where own-centring runs); continuous / covariate designs keep
# c14's kernels.
# it60 read locally (2026-11-04) -- REJECTED. On bank1's 20 msq cells it did what it52's lin0 did (TP R00 107 -> 115,
# R11 100 -> 121), but its own centre does not save the bloom: bloom x4 tongue / stool 35-77 / 1-58 false calls per cell
# (c14 0-3), house R25 stool 72-211 (c14 0-18); the selection takes lin0 for 260-380 taxa per failing cell. A bloom scales
# the other taxa by a factor that varies by sample; on the log scale with size factors that is one shift per sample, on
# the plain proportion scale each taxon's group difference depends on how it co-varies with the bloomer, so no single
# centre removes it. LDM's frequency scale is the same statistic: on the dev suite's bloom x4 cells (srv9, 10 cells) LDM
# has 198.6 false calls per cell (FDR 0.90; ZicoSeq 0.22). c14's msq gap to LDM is, in this part, the price of bloom robustness.
#
# Default options reproduce every earlier candidate exactly: .root_shift_med is 60_erd.R's .root_shift, copied.
.root_shift_med <- function(T, k, iters = 8L, tol = 1e-3) {
  f <- function(d) .colmed(T + outer(k, d)); s <- stats::median(k); if (!is.finite(s) || abs(s) < 1e-8) return(rep(0, ncol(T)))
  d0 <- rep(0, ncol(T)); f0 <- f(d0); d1 <- -f0 / s; f1 <- f(d1)
  for (it in seq_len(iters)) { den <- f1 - f0; dn <- ifelse(abs(den) > 1e-12, d1 - f1 * (d1 - d0) / den, d1 - f1 / s)
    dn <- pmin(pmax(dn, -3), 3); d0 <- d1; f0 <- f1; d1 <- dn; f1 <- f(d1); if (max(abs(f1)) < tol) break }
  d1
}
.colmed_na <- function(M) { v <- if (requireNamespace("matrixStats", quietly = TRUE)) matrixStats::colMedians(M, na.rm = TRUE) else
  apply(M, 2, stats::median, na.rm = TRUE); v[!is.finite(v)] <- 0; v }
.root_shift_trim <- function(T, k, c0 = 2, rounds = 2L, iters = 8L, tol = 1e-3) {   # median root on the reference rows (|centred t| <= c0)
  d <- .root_shift_med(T, k, iters, tol); s <- stats::median(k); if (!is.finite(s) || abs(s) < 1e-8) return(d)
  for (r in seq_len(rounds)) {
    keep <- abs(T + outer(k, d)) <= c0
    f <- function(dd) { M <- T + outer(k, dd); M[!keep] <- NA; .colmed_na(M) }
    d0 <- d; f0 <- f(d0); d1 <- d0 - f0 / s; f1 <- f(d1)
    for (it in seq_len(iters)) { den <- f1 - f0; dn <- ifelse(abs(den) > 1e-12, d1 - f1 * (d1 - d0) / den, d1 - f1 / s)
      dn <- pmin(pmax(dn, -3), 3); d0 <- d1; f0 <- f1; d1 <- dn; f1 <- f(d1); if (max(abs(f1)) < tol) break }
    d <- d1 }
  d
}
# "lad": the shift minimising sum_j |k_j| |T_j / k_j + d| = sum_j |T_j + k_j d| over the reference rows -- a weighted median of
# -T / k with weights |k|, re-solved twice on the rows with |centred t| <= pursue.trim. The median rule weights every taxon
# alike, but a taxon's sensitivity to the centre (|kap|) is largest for abundant taxa, so a centre error shows there most:
# local bias check (73 banked cells, 2026-11-04): median centred t among null taxa in the most abundant third, mean |.| over
# ten settings 0.119 (median rule), 0.090 (trimmed median), 0.060 (lad); all-up house -0.20 / -0.10 / -0.08, all-up msq -0.22 /
# -0.14 / -0.07, bloom x4 -0.09 / -0.06 / 0.00, balanced house R00 +0.03 / +0.08 / +0.08.
# Calls on the same banked cells (c14 -> lad, TP / FP summed; tongue, 3 reps, and bank1's 20 SimulateMSeq cells): all-up house
# R10 21/0 -> 27/1, R11 16/0 -> 18/0, B_allup 21/5 -> 23/4, msq R11 100/10 -> 106/13, bloom x4 15/4 -> 15/3; balanced house
# R00 24/5 -> 23/5, B_ref 24/2 -> 22/2, msq R00 107/8 -> 108/8; global nulls (house R06, B_null) no rejection either way. The
# price: on balanced designs the weighting leaves null taxa a small positive offset (+0.04 to +0.09), worth a few discoveries.
.wmed_cols <- function(V, W) {                                     # weighted median of each column (V, W: rows x columns)
  n <- nrow(V); B <- ncol(V); o <- apply(V, 2, order); ix <- cbind(as.vector(o), rep(seq_len(B), each = n))
  Vs <- matrix(V[ix], n, B); cw <- apply(matrix(W[ix], n, B), 2, cumsum); if (n == 1L) cw <- matrix(cw, 1L)
  pos <- colSums(cw < rep(cw[n, ] / 2, each = n)) + 1L; Vs[cbind(pmin(pos, n), seq_len(B))]
}
.root_shift_lad <- function(T, k, c0 = 2, rounds = 2L) {
  ok <- is.finite(k) & abs(k) > 1e-6; if (sum(ok) < 3L) return(.root_shift_med(T, k))
  T <- T[ok, , drop = FALSE]; k <- k[ok]; V <- -T / k; W <- matrix(abs(k), nrow(T), ncol(T))
  d <- .wmed_cols(V, W)
  for (r in seq_len(rounds)) { keep <- abs(T + outer(k, d)) <= c0; Wk <- W * keep; none <- colSums(Wk) <= 0; if (any(none)) Wk[, none] <- W[, none]
    d <- .wmed_cols(V, Wk) }
  pmin(pmax(d, -3), 3)
}
.root_shift <- function(T, k, iters = 8L, tol = 1e-3) {
  ce <- getOption("pursue.centre", "median")
  if (identical(ce, "lad")) return(.root_shift_lad(T, k, c0 = getOption("pursue.trim", 2)))
  if (identical(ce, "trim")) return(.root_shift_trim(T, k, c0 = getOption("pursue.trim", 2), iters = iters, tol = tol))
  .root_shift_med(T, k, iters, tol)
}
.c15 <- function(lin0 = FALSE, centre = "median", trim = 2) function(counts, meta, formula, tested_term) {
  d <- .epg_design(meta, formula, tested_term); two <- !is.null(d) && d$binary && is.null(d$Z) && is.null(d$cl)
  pk <- if (lin0 && two) c("pos", "lin0") else "pos"
  .with_opts({
    v <- if (.is_paired(meta, formula, tested_term)) .epp_fit(counts, meta, formula, tested_term, lin = TRUE, kapfree = 20L, rho1 = TRUE) else
      .epg_fit(counts, meta, formula, tested_term, paths = "cont", adet = TRUE, nodet_conf = TRUE, pk = pk, minsw = 2L)
    data.frame(feature = v$feature, p = v$p_full, q = .storey_gq(v$p_full), estimate = v$est) },
    pursue.flstat = "cor", pursue.own = pk, pursue.centre = centre, pursue.trim = trim)
}
register_candidate("eperm_c15c", .c15(centre = "lad"), notes = "it59: c14 with the compositional centre as a |kap|-weighted median on the reference rows (lad)")
register_candidate("eperm_c15t", .c15(centre = "trim"), notes = "it59 ablation: c14 with the plain median re-solved on the reference rows (|centred t| <= 2, twice)")
register_candidate("eperm_c15l", .c15(lin0 = TRUE), notes = "it60 (rejected: bloom / R25): c14 + the raw plain-proportion kernel (lin0) centred on its own median taxon")
register_candidate("eperm_c15", .c15(lin0 = TRUE, centre = "lad"), notes = "it59 + it60 (rejected with it60)")
