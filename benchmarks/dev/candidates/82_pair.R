# --- it57: exact within-subject path for paired designs (repeated measures, exposure changing within subject) -------
# srv30 (hmp_tongue with HMP subjects, 40 subjects x 2 visits, one visit exposed at random): c12 is valid but far too
# conservative here -- null-taxon FPR at p < 0.05 0.014 (nominal 0.05), no false call in 40 null cells -- and finds 4.35
# TP per cell against 5.2-5.4 for LDM (cluster.id), ZicoSeq (strata) and LinDA (random intercept). A within-subject
# exposure sends c12 to the large-sample path (subject-clustered CR1 t) with the label-free pair graph: each sample is
# compared with ~60 random partners, almost all of them other people. The subject-clustered variance keeps that valid,
# but the partners' between-person variation stays in every score, so the pairing is not used: c12's null FPR (0.014)
# sits with the methods that ignore subjects (0.015-0.025).
#
# The pairwise statistic has a natural paired form. Take the pair graph to be the within-subject pairs only (all of a
# subject's visit pairs: label-free, so the within-subject permutation is exact; same-label pairs cancel), compare each pair at its common depth as everywhere else (same kernels, size factors, thinning by
# exp(gamma x) for the compositional centre), and permute the exposure WITHIN subject -- for two visits, a sign flip of
# each subject's kernel difference. Under H0 (a subject's visits exchangeable with respect to the exposure, which a
# randomised within-subject assignment guarantees) that is exact; the Welch t of the scores over the flipped labels is
# the paired t up to a constant. The centre is re-estimated in every permutation; kernels, LOO BH selection (margin 1,
# minsw 2), the pos kernel on its own centre and gated Storey q are c12's. Subjects without both exposure levels carry
# no within-subject contrast and are left out. Every other design runs c12 unchanged.
.epp_memo <- new.env()
.epp_fit <- function(counts, meta, formula, tested_term, thr = 1e-3, margin = 1L, B = 4000L, B2 = 40000L, winsor = 0.03,
                     adet = TRUE, pk = "pos", minsw = 2L) {
  key <- list(counts, meta, formula, tested_term, thr, margin, B, B2, winsor, adet, pk, minsw, getOption("pursue.own"))
  if (!is.null(.epp_memo$key) && identical(.epp_memo$key, key)) return(.epp_memo$val)
  d <- .epg_design(meta, formula, tested_term); depth0 <- if (!is.null(meta$depth)) meta$depth else colSums(counts)
  cl0 <- d$cl; both <- tapply(d$x, cl0, function(v) length(unique(v)) == 2L); use <- cl0 %in% names(both)[both]
  feats <- rownames(counts); counts <- counts[, use, drop = FALSE]; depth <- depth0[use]; x <- d$x[use]; g1 <- d$g1[use]
  cl <- droplevels(cl0[use]); n <- ncol(counts)
  blk <- lapply(split(seq_len(n), cl), function(ii) list(i = ii, k = sum(x[ii] == 1)))
  # every within-subject pair, whatever its labels: the graph must not depend on the exposure, or permuting the labels
  # within subject is not exact (same-label pairs cancel in the group sum; with two visits it is the one pair)
  E <- do.call(rbind, lapply(blk, function(b) t(utils::combn(b$i, 2L))))
  two <- all(vapply(blk, function(b) length(b$i) == 2L && b$k == 1L, logical(1)))       # two visits, one exposed: a sign flip
  gen <- if (two) { A <- vapply(blk, function(b) b$i[1], 1L); Z <- vapply(blk, function(b) b$i[2], 1L); nb <- length(blk)
    function() { f <- stats::runif(nb) < 0.5; A * f + Z * (!f) } } else
    function() unlist(lapply(blk, function(b) b$i[sample.int(length(b$i), b$k)]), use.names = FALSE)
  rho <- round(.rho_design(depth, g1), 2)
  sf <- if (rho == 1) .size_factors(counts, depth) else rep(0, n)
  raw <- counts; if (winsor > 0 && rho == 1) counts <- .winsorize_counts(counts, depth, winsor)
  rows <- .with_seed(7L, if (nrow(counts) > 300L) sort(sample.int(nrow(counts), 300L)) else seq_len(nrow(counts)))
  ct <- counts[rows, , drop = FALSE]; ct <- ct[rowSums(ct > 0) >= 3, , drop = FALSE]
  ok <- rowSums(counts > 0) >= 3
  hs <- c(det = "det", log = "log", sqrt = "sqrt"); combos <- list(lead = c("det", "log"), det = "det", log = "log", sqrt = "sqrt"); alts <- c(if (adet) "det", "log", "sqrt")
  for (h in pk) { hs[[h]] <- h; combos[[h]] <- h; alts <- c(alts, h) }
  sc <- function(cnt, gm, h) { bb <- exp(gm * x + sf)
    if (h == "pos") .perm_scores_pos(cnt, depth, bb, E, rho = rho) else .perm_scores(cnt, depth, bb, h, E, rho = rho) }
  G1 <- matrix(as.numeric(g1), n, 1); n1 <- sum(g1); n0 <- n - n1
  tw <- function(a) { a <- a - rowMeans(a); v <- .perm_tstat(a, a^2, G1, n1, n0)[, 1]; v[!is.finite(v)] <- 0; v }
  tstat <- function(SS) matrix(sapply(hs, function(h) tw(SS[[h]])), ncol = length(hs), dimnames = list(NULL, hs))
  medt <- function(gm) stats::median(tw(sc(ct, gm, "log")))
  flo <- medt(-2); fhi <- medt(2)
  gm <- if (is.finite(flo) && is.finite(fhi) && sign(flo) != sign(fhi)) stats::uniroot(medt, c(-2, 2), f.lower = flo, f.upper = fhi, tol = 2e-3)$root else 0
  cok <- counts[ok, , drop = FALSE]
  S <- lapply(hs, function(h) sc(cok, gm, h)); S0 <- lapply(hs, function(h) sc(cok, 0, h))
  t1 <- tstat(S); t0 <- tstat(S0)
  kap <- if (abs(gm) >= 0.02) (t1 - t0) / gm else (tstat(lapply(hs, function(h) sc(cok, 0.05, h))) - t0) / 0.05
  kap[!is.finite(kap)] <- 0
  med <- which(rownames(cok) %in% rownames(ct))
  tobs <- t1 + kap * .root_shift(t1[med, "log", drop = FALSE], kap[med, "log"])
  own <- intersect(getOption("pursue.own", character(0)), setdiff(colnames(t1), "log"))
  for (h in own) tobs[, h] <- t1[, h] + kap[, h] * .root_shift(t1[med, h, drop = FALSE], kap[med, h])
  pm <- .perm_multi(S0, g1, combos, B = B, B2 = B2, tobs = tobs, kap = kap, med = med, own = own, gen = gen)
  P <- matrix(NA_real_, length(feats), ncol(pm$p), dimnames = list(feats, colnames(pm$p))); P[ok, ] <- pm$p
  pl <- P[, "lead"]; Q <- P[, c("lead", alts), drop = FALSE]
  hit <- apply(Q, 2, function(p) { q <- rep(1, length(p)); f <- is.finite(p); q[f] <- stats::p.adjust(p[f], "BH"); q <= 0.05 })
  loo <- sweep(-hit, 2, colSums(hit), "+"); la <- loo[, alts, drop = FALSE]; ia <- max.col(la, ties.method = "last"); altn <- la[cbind(seq_len(nrow(la)), ia)]
  usek <- altn >= loo[, "lead"] + margin & altn >= minsw; pfull <- ifelse(usek, Q[cbind(seq_len(nrow(Q)), 1L + ia)], pl)
  est <- rep(NA_real_, length(feats)); est[ok] <- rowSums(S$log[, g1, drop = FALSE]) / sum(x[E[, 1]] != x[E[, 2]]) / log(2)
  val <- list(feature = feats, p_full = pfull, p_ms = pl, est = est, gam = gm, rho = rho, perm = TRUE, fallback = FALSE, path = "pair",
              n_subj = length(blk), choice = ifelse(usek, alts[ia], "lead"), P = P)
  .epp_memo$key <- key; .epp_memo$val <- val; val
}
# srv31 (house, > 2 visits): the within-subject path is valid with any visit pattern (null FPR 0.050), but it uses only
# within-subject contrasts. Where people differ in how much of their time is exposed (HW5: 1-4 of 5 visits) or some are
# never exposed (HWv), the exposure also varies BETWEEN people, and that information is thrown away (c13 5.80 / 2.40 TP
# vs c12 6.25 / 5.25). Dispatch on the design alone (labels and subjects, never the counts, so each path keeps its own
# validity): the within-subject path only when >= 90% of the exposure's variance lies within subjects -- every person
# exposed for about the same share of visits (paired, crossover, balanced pre/post); otherwise c12's clustered path.
.within_share <- function(x, cl) { v <- sum((x - mean(x))^2); if (v <= 0) return(0); sum((x - stats::ave(x, cl))^2) / v }
.is_paired <- function(meta, formula, tested_term) {
  d <- .epg_design(meta, formula, tested_term)
  !is.null(d) && d$binary && is.null(d$Z) && !is.null(d$cl) && !d$cl_ok &&
    sum(tapply(d$x, d$cl, function(v) length(unique(v)) == 2L)) >= 6L && .within_share(d$x, d$cl) >= 0.9
}
register_candidate("eperm_c13", function(counts, meta, formula, tested_term) {
  .with_opts({
    v <- if (.is_paired(meta, formula, tested_term)) .epp_fit(counts, meta, formula, tested_term) else
      .epg_fit(counts, meta, formula, tested_term, paths = "cont", adet = TRUE, nodet_conf = TRUE, pk = "pos", minsw = 2L)
    data.frame(feature = v$feature, p = v$p_full, q = .storey_gq(v$p_full), estimate = v$est) }, pursue.flstat = "cor", pursue.own = "pos")
}, notes = "it57: eperm_c12 plus an exact within-subject permutation path (within-subject pairs, labels permuted within subject) for paired designs")
