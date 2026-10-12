# --- it51: abundance given presence at each sample's own depth (Phase 2, target 3: low-signal sd2 / sps) -----------------
# p14: PURSUE 0.2 still beats c5q in 11 sd2 and 11 sps regimes -- the low-signal ones (R00, R08, R09, R14, R15 ...),
# while c5q leads where the signal is strong (R02, R03, R13). 0.2's abundance arm fits log abundance among positives at
# full depth; our pos kernel compares each pair at the pair's common depth, which throws reads away. On the banked sd2
# cells a direct Welch test of log RA among positives (full depth, no compositional correction) found 13.2 / 17.6 TP per
# cell at R08 / R13 against the pairwise kernel's 10.0 / 17.0 under permutation.
# The per-sample version uses it43's construction for exchangeable designs (two groups, no covariates, depth imbalance
# at chance level): thin only the group whose null taxa are inflated, by exp(-|gamma|) -- binomial thinning matches the
# whole sampling distribution, so validity needs nothing else -- keep each sample at its own depth D_i and score
#   s_i = E[log(Y/D_i); Y > 0] - c P(Y > 0),   c = pooled mean log RA among positives (label-free)
# = P(present) x (log RA given presence - c), residualised on log D_i across samples. A linear statistic in fixed
# sample scores, so the permutation stays exact; gamma is re-estimated in every permutation as for the other options.
# Offered to the selection beside the pairwise pos kernel.
# Result (banked sd2 R08 / R13, 10 cells): pps alone 9.2 / 13.8 TP per cell against the pairwise pos kernel's 10.8 / 14.8;
# c7p 12.6 / 13.6 vs c5q 12.4 / 13.4, with one cell of 5 false calls where pps was chosen. Not promoted. The 13.2 / 17.6
# that motivated it did not reproduce: a positives-only Welch test on these cells finds 7.8 / 12.2 (t reference) and
# 10.4 / 14.8 (permutation, with false calls) -- the pairwise kernel already matches a direct test.
register_candidate("eperm_c7p", function(counts, meta, formula, tested_term) {
  v <- .epg_fit(counts, meta, formula, tested_term, paths = "cont", adet = TRUE, nodet_conf = TRUE, pk = "pos", minsw = 2L, ps = "pps")
  data.frame(feature = v$feature, p = v$p_full, q = .storey_gq(v$p_full), estimate = v$est)
}, notes = "it51: eperm_c5q + per-sample abundance-given-presence scale (pps) in exchangeable two-group designs")
