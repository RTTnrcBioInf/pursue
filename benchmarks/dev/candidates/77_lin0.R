# --- it52: the frequency scale, unadjusted, and a sqrt / log omnibus as selection options (Phase 2, target 3: msq) --------
# probe1: LDM's msq edge is its frequency-scale (relative abundance) test under permutation. Our linear kernel looked
# hopeless (5.6 TP per cell on banked msq R00 tongue vs LDM's 9.6) until bank tests showed why: size factors (and the
# 3% winsorising) wreck the linear scale on SimulateMSeq's resampled real compositions. "lin0" (it50's raw-count,
# no-size-factor variant of the linear kernel; compositional centre kept, so it still tests absolute change) reaches
# 10.5 alone there. Second, "sl" = max(|t_sqrt|, |t_log|) under the joint permutation -- LDM's omnibus construction
# on our kernels -- for taxa that sit between the two scales. Post hoc on the banked msq / sd2 cells (selection rule
# as in c5, gated Storey; TP per cell; ZicoSeq / LDM in brackets):
#              msq R00 tongue   R00 twinsuk   R11 tongue   R11 twinsuk   sd2 R08 / R13
#   c5q             8.75           13.2           6.2          14.6       12.4 / 14.75
#   + lin0         10.8            12.4           6.8          14.2       12.4 / 14.8
#   + lin0, sl     10.8            13.8           6.8          15.4       12.4 / 14.8
#   (ZicoSeq)       9.5            12.8           6.8          14.6        4.8 / 6.25
#   (LDM)          10.25           13.8           7.4          13.2        6.8 / 8.75
# False positives equal or lower. 4-5 cells per set, so the dev suite decides.
register_candidate("eperm_c8", function(counts, meta, formula, tested_term) {
  v <- .epg_fit(counts, meta, formula, tested_term, paths = "cont", adet = TRUE, nodet_conf = TRUE, pk = c("pos", "lin0"), minsw = 2L)
  data.frame(feature = v$feature, p = v$p_full, q = .storey_gq(v$p_full), estimate = v$est)
}, notes = "it52: eperm_c5q + the unadjusted linear (frequency) kernel as a selection option")
register_candidate("eperm_c8s", function(counts, meta, formula, tested_term) {
  v <- .epg_fit(counts, meta, formula, tested_term, paths = "cont", adet = TRUE, nodet_conf = TRUE, pk = c("pos", "lin0"), minsw = 2L, cmb = list(sl = c("sqrt", "log")))
  data.frame(feature = v$feature, p = v$p_full, q = .storey_gq(v$p_full), estimate = v$est)
}, notes = "it52: eperm_c8 + max(|t_sqrt|, |t_log|) omnibus as a selection option")
