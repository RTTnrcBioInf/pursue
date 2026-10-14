# --- it56: the abundance-given-presence kernel centred on its own median taxon (Phase 2, spread calibration) --------
# srv27 (spread-only implants on real samples; every call false): c10 fails H_all -- every taxon 3x more variable in the
# 30% group, R25's design -- on twinsuk_stool only (P(any) 0.79, 54 calls per cell; tongue 0.05). Every comparator but
# logistic presence fails it too. Local diagnosis on stool reps 18 and 7: the selection moves every taxon to the pos
# kernel (it47), whose p < 0.05 rate is 0.23 / 0.18 against 0.05-0.07 for lead, det, log and sqrt. The pos bulk is
# shifted, not over-dispersed: observed t median -1.16 (sd 1.13) under the run's centring, which is the log kernel's
# compositional centre (gamma 0.39 here) carried to every kernel by its sensitivity kap. A spread change moves each
# kernel's typical value differently -- det +0.28, sqrt +0.06, pos -1.16 in t units when log is held at 0 -- so one
# centre cannot serve them all. Centred on its own median taxon, pos's bulk sits at -0.12 and its BH calls (normal p)
# go 88 -> 0; on four banked sd2 cells its true positives are unchanged (8 / 8, 9 / 8, 20 / 21, 13 / 13).
# The rule is the same one the log kernel already uses -- the median taxon does not change -- applied on pos's own
# scale, in the observed statistic and in every permutation (option pursue.own, read by .eps_fit and .perm_multi;
# default none, so every earlier candidate is bit-identical). Two-group exchangeable path only: that is where pos runs
# and where H_all failed; the continuous setting passed.
# Robert (2026-10-14): a taxon more variable in one group but unchanged in typical abundance is not to be called.
register_candidate("eperm_c12", function(counts, meta, formula, tested_term) {
  .with_opts({ v <- .epg_fit(counts, meta, formula, tested_term, paths = "cont", adet = TRUE, nodet_conf = TRUE, pk = "pos", minsw = 2L)
    data.frame(feature = v$feature, p = v$p_full, q = .storey_gq(v$p_full), estimate = v$est) }, pursue.flstat = "cor", pursue.own = "pos")
}, notes = "it56: eperm_c10 with the pos kernel centred on its own median taxon (not the log kernel's centre)")
