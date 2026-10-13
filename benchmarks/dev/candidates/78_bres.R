# --- it53: permutation p-values fine enough for the BH threshold at large m (Phase 2, targets 2 / 3: implants, B_cont) ---
# srv24 / bank3: on B_cont c5q finds 13.0 TP per cell against ZicoSeq's 22.05 (axis B: 1.12 vs 2.17). The kernels are not
# the problem: on banked twinsuk B_cont cells our sqrt kernel's observed t correlates 0.98 with a plain regression of
# sqrt(RA) on exposure and is nearly as large (median |t| on DA taxa 1.8 vs 2.0), and BH on its normal p finds 15-22 TP
# in a cell where the permutation p finds 2. The permutation null is fine (per-taxon sd of permuted t 1.03). It is the
# resolution: implant cells keep the template's whole feature set (m ~ 2200 on twinsuk), where BH's first threshold is
# 0.05 / m = 2.3e-5 -- below the floor of B + B2 = 44 000 permutations (1 / 44 001). The strongest taxa all sit at the
# floor, BH never starts. (ZicoSeq pools the permutations of all taxa, which escapes the floor; asymptotic tests have
# none -- one reason logistic presence leads several large-m implant rows.)
# Fix: a third stage for the taxa still below cmin exceedances, up to bres x m permutations in total (bres 200 = ten
# exceedances at p = alpha / m). Only the extreme taxa get it, so the cost is seconds; m <= 220 changes nothing.
# Read from options("pursue.bres") inside both permutation engines (.perm_multi, .epg_perm_fl); default 0 = off.
# Result (bank3, twinsuk B_cont, 4 cells): TP 24 / 29 / 26 / 0 -> 24 / 27 / 23 / 5, at 3-4 extra minutes per cell. The
# floor was not what held B_cont back -- the studentisation was (it54, 79_flcor.R). Not promoted; the stage stays
# available (default off).
.with_bres <- function(expr, bres = 200) { old <- options(pursue.bres = bres); .clear_fit_memos(); on.exit({ options(old); .clear_fit_memos() }); expr }
register_candidate("eperm_c9", function(counts, meta, formula, tested_term) {
  .with_bres({ v <- .epg_fit(counts, meta, formula, tested_term, paths = "cont", adet = TRUE, nodet_conf = TRUE, pk = "pos", minsw = 2L)
    data.frame(feature = v$feature, p = v$p_full, q = .storey_gq(v$p_full), estimate = v$est) })
}, notes = "it53: eperm_c5q with permutation p resolved down to ~0.05 / m (bres 200)")
