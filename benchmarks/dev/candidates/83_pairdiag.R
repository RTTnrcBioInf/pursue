# --- srv32: what makes c13's within-subject path non-exact? -----------------------------------------------------------
# conf2 (held-out, read 2026-10-30): on hmp_skin_ear's within-subject nulls c13 had P(any false call) 0.45 and 3.8 false
# calls per cell (per-taxon FPR 0.059); on the other five templates <= 0.05. An exact sign-flip test cannot do that, so a
# step of the path depends on the observed labels in a way the permutations do not reproduce. Suspects, one switch each
# (82_pair.R, .epp_fit; every other design runs c12 unchanged):
#   lin  -- the observed statistic is re-thinned at the fitted centre gamma (data level, nonlinear) while every permutation
#           is re-centred by the linearisation t(0) + kap x shift: with one partner per sample the thinning is pair-specific
#           (min of the two depths), so the linearisation can be poor and the two disagree. lin computes the observed
#           statistic by the permutations' own formula.
#   kapfree -- kap (each taxon's sensitivity to gamma) is taken at the observed labels and reused for every permutation;
#           with one partner per sample it depends on which member is exposed. kapfree averages it over 20 flips.
#   rho1 -- rho (thinning factor) is chosen from the observed depth-exposure imbalance and kept for every permutation.
#   nosel -- per-taxon kernel choice from the other taxa's discoveries (LOO), a data-driven choice the p-value ignores.
#   g0  -- no compositional centre at all: fixed scores, so exact by construction (the reference; loses the centring).
# lin + kapfree makes the centring a fixed function of the data applied identically to every labelling; with rho1 the
# path has no label-dependent step left except the kernel choice (c13fx).
.pdiag <- function(...) { a <- list(...); function(counts, meta, formula, tested_term) .with_opts({
  v <- if (.is_paired(meta, formula, tested_term)) do.call(.epp_fit, c(list(counts, meta, formula, tested_term), a)) else
    .epg_fit(counts, meta, formula, tested_term, paths = "cont", adet = TRUE, nodet_conf = TRUE, pk = "pos", minsw = 2L)
  data.frame(feature = v$feature, p = v$p_full, q = .storey_gq(v$p_full), estimate = v$est) }, pursue.flstat = "cor", pursue.own = "pos") }
register_candidate("eperm_c13lin", .pdiag(lin = TRUE), notes = "srv32: c13, observed statistic by the permutations' linearised re-centring")
register_candidate("eperm_c13lf", .pdiag(lin = TRUE, kapfree = 20L), notes = "srv32: c13lin + label-free kap (mean over 20 within-subject flips)")
register_candidate("eperm_c13r1", .pdiag(rho1 = TRUE), notes = "srv32: c13 with rho fixed at 1 on the paired path")
register_candidate("eperm_c13ns", .pdiag(nosel = TRUE), notes = "srv32: c13, lead kernel only (no per-taxon kernel choice) on the paired path")
register_candidate("eperm_c13g0", .pdiag(g0 = TRUE), notes = "srv32: c13 without the compositional centre on the paired path (exact reference)")
register_candidate("eperm_c13fx", .pdiag(lin = TRUE, kapfree = 20L, rho1 = TRUE), notes = "srv32: lin + kapfree + rho1 (label-free centring and thinning)")
