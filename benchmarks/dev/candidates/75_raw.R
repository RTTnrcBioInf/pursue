# --- it50: unadjusted scales as selection options -- raw counts, no size factors (Phase 2, target 3: msq) ---------------
# p14: ZicoSeq leads msq in 19 of 21 regimes by a steady 10-30% (strict TPR 0.105 vs 0.092). On the banked msq tuning
# cells (bank1) c5q trails ZicoSeq on hmp_tongue (R00 8.4 vs 9.4, R11 6.2 vs 6.8 TP per cell; level on twinsuk). The
# taxa only ZicoSeq calls are prevalent (median prevalence 0.74) and mid-abundant (log10 RA -2.5), with our p ~10^-2.4
# against its 10^-3; the ones only we call are rarer (0.43, -3.5). On the same cells, eperm_cs1's selection without size
# factors and winsorising reaches 9.8 (sqrt kernel alone 9.2; LDM 10.0) -- SimulateMSeq resamples real compositions,
# where a per-sample median-ratio factor estimated from a few dozen prevalent taxa adds noise instead of removing a
# common scaling. srv21 showed the opposite elsewhere (dropping them everywhere cost the continuous / covariate paths and
# raised the x4 bloom FDR), so the data should choose: "<kernel>0" is the same pairwise kernel on raw (unwinsorised)
# counts with the compositional centre only (no size factors), offered to the leave-one-out selection beside the
# adjusted one. Label-free, so the permutation stays exact; the floor (minsw 2) guards the extra options.
# Result (banked msq hmp_tongue, 10 cells): c7 8.6 / 6.0 TP per cell at R00 / R11 vs c5q 8.4 / 6.2 -- sqrt0 alone is no
# better than sqrt with size factors. Not sent to the server. The raw-count kernel variants stay: lin0 (77_lin0.R) is one.
register_candidate("eperm_c7", function(counts, meta, formula, tested_term) {
  v <- .epg_fit(counts, meta, formula, tested_term, paths = "cont", adet = TRUE, nodet_conf = TRUE, pk = c("pos", "sqrt0", "log0"), minsw = 2L)
  data.frame(feature = v$feature, p = v$p_full, q = .storey_gq(v$p_full), estimate = v$est)
}, notes = "it50: eperm_c5q + sqrt and log kernels on raw counts without size factors as selection options")
