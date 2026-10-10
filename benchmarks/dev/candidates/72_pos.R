# --- it47: abundance given presence as a selection option (Phase 2, target 3: sd2) --------------------------------
# Where it came from. The confirmation run (fresh seed) and p13 agree that on sd2 PURSUE 0.2's abundance arm has
# 3-5x eperm_cs1's calibrated power (R03 0.178 vs 0.033, R13 0.150 vs 0.030; LDM 0.057 / 0.055) -- and 0.2 models
# log abundance among the non-zero counts. SparseDOSSA2 puts DA in the log-normal part of a zero-inflated model and
# leaves presence alone, so every kernel that mixes the (unchanged) zeros in -- detection, log(1 + Y), sqrt --
# dilutes the signal. On the banked sd2 cells (tuning template hmp_tongue, never the evaluation data) a Welch test
# of log relative abundance among positives found 13.2 / 17.6 TP per cell at sd2 R08 / R13 against ours 3.2 / 5.8 and
# LDM 6.8 / 9.8; on msq it is weaker than our scales (6.2 vs 12.6 on twinsuk R00), so it is an option, not a default.
# The kernel keeps the exact pairwise frame (it30's): at the pair's common effective depth,
#   A_ij = f_j L_i - f_i L_j,   f = P(Y_D > 0),   L = E[log(Y_D / (D b)); Y_D > 0]   (expected rarefaction, exact)
# = f_i f_j (mean log abundance given presence of i - of j). Both members are thinned to the same effective depth,
# so A_ij is antisymmetric in distribution under H0 and the linear permutation statistic stays exact. Pairwise on
# the same sd2 cells (Welch normal p): 11.4 / 15.4 TP; msq R00 tongue 6.4 vs sqrt 9.6. Joins the leave-one-out BH
# selection; the centre is re-estimated in every permutation as for every other kernel.
.perm_scores_pos <- function(counts, depth, b, E, rho = 1, chunk = 4000L) {
  Y <- as.matrix(counts); storage.mode(Y) <- "double"; m <- nrow(Y); n <- ncol(Y); S <- matrix(0, m, n)
  for (st in seq(1L, nrow(E), by = chunk)) { e <- E[st:min(nrow(E), st + chunk - 1L), , drop = FALSE]; i <- e[, 1]; j <- e[, 2]
    De <- rho * pmin(depth[i] * b[i], depth[j] * b[j]); Di <- pmax(floor(De / b[i]), 1); Dj <- pmax(floor(De / b[j]), 1)
    fi <- .erd_f_var(Y[, i, drop = FALSE], depth[i], Di); fj <- .erd_f_var(Y[, j, drop = FALSE], depth[j], Dj)
    Li <- .erl_f_var(Y[, i, drop = FALSE], depth[i], Di, off = 0) - fi * rep(log(Di * b[i]), each = m)
    Lj <- .erl_f_var(Y[, j, drop = FALSE], depth[j], Dj, off = 0) - fj * rep(log(Dj * b[j]), each = m)
    A <- fj * Li - fi * Lj
    Inc <- matrix(0, nrow(e), n); Inc[cbind(seq_len(nrow(e)), i)] <- 1; Inc[cbind(seq_len(nrow(e)), j)] <- -1; S <- S + A %*% Inc }
  S
}
register_candidate("eperm_c4", function(counts, meta, formula, tested_term) {
  v <- .epg_fit(counts, meta, formula, tested_term, paths = "cont", adet = TRUE, nodet_conf = TRUE, pk = "pos"); data.frame(feature = v$feature, p = v$p_full, estimate = v$est)
}, notes = "it47: eperm_c3 + abundance-given-presence kernel (pos) as a selection option")
register_candidate("eperm_c4q", function(counts, meta, formula, tested_term) {
  v <- .epg_fit(counts, meta, formula, tested_term, paths = "cont", adet = TRUE, nodet_conf = TRUE, pk = "pos")
  data.frame(feature = v$feature, p = v$p_full, q = .storey_q(v$p_full), estimate = v$est)
}, notes = "it47 + it45: eperm_c4 with Storey's adaptive BH")
