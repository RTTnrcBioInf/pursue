# --- it58: c14 = c12 + the corrected within-subject path (srv32) -------------------------------------------------------
# srv32 (development data: tongue, buccal, throat, baboons, mixed-site and thinned pairs) found two faults in c13's exact
# within-subject path: the observed statistic was re-thinned at the fitted centre while every permutation was re-centred
# by a linearisation (baboons, a scrambled second site), and rho was chosen from the observed labels, so a chance depth
# imbalance switched off size factors and winsorisation and the selected kernel picked up compositional spill-over
# (throat). c14 computes the observed statistic by the permutations' own formula, with kap averaged over 20 within-
# subject flips and rho fixed at 1: every step of the centring and thinning is a fixed function of the data, applied
# identically to each labelling (per-taxon kernel choice as in c12). srv32: null P(any) 0.056 (c13 0.084; exact
# reference without centring 0.056), signal FDR 0.076 (c13 0.169) at the same power (TP 6.94 vs 6.89). Designs that are
# not within-subject run c12 unchanged, so A-C and D/E are c13's.
register_candidate("eperm_c14", .pdiag(lin = TRUE, kapfree = 20L, rho1 = TRUE),
                   notes = "it58: c12 + within-subject path with label-free centring (linearised observed statistic, kap over 20 flips) and rho 1")
