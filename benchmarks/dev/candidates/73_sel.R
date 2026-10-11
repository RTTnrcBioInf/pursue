# --- it48: an option must earn its switch; adaptivity only once BH has found something (Phase 2, target 4) ------------
# srv20n: eperm_c4's extra null false discoveries (house / mid R06, B_null, sd2 R06; 4 of 440 cells, c3 0 in each)
# are all exactly two per cell -- the signature of the leave-one-out selection. When an alternative option has
# exactly two BH discoveries and the lead none, each of the two sees one discovery among the others, switches
# (1 >= 0 + margin 1) and is called; one spurious pair supports itself. Every option added to the selection is one
# more chance of such a pair (c3 null P(any) 0.052 -> c4 0.061), and the same mechanism adds false positives in
# low-signal cells (srv20, c4 - c3: mid R00 +4, house R15 +2, B_abund +2, sps R00 +2). Rule: an alternative needs at
# least `minsw` discoveries among the other taxa (minsw = 2: an option is used only once it finds three), on top of
# the margin over the lead. minsw = 1 is the previous behaviour (margin 1 already implies one).
# Adaptive BH (it45) multiplies BH's q by pi0 <= 1. Under a global null the Storey estimate dips below 1 about half the
# time, which is where c2q / c4q add null cells (srv19n, srv20n: +1 to +3 in 440). Gate: use pi0 only when plain BH
# already rejects at q <= 0.05. Under a global null any BH rejection is already a false discovery, so P(any false
# discovery) is exactly BH's; on srv20's signal cells the gate keeps 295 of c4q's 296 extra true positives.
.storey_gq <- function(p, alpha = 0.05) { f <- is.finite(p); q <- rep(NA_real_, length(p)); if (!any(f)) return(q)
  q[f] <- stats::p.adjust(p[f], "BH"); if (!any(q[f] <= alpha)) return(q)
  .storey_q(p) }
register_candidate("eperm_c5", function(counts, meta, formula, tested_term) {
  v <- .epg_fit(counts, meta, formula, tested_term, paths = "cont", adet = TRUE, nodet_conf = TRUE, pk = "pos", minsw = 2L); data.frame(feature = v$feature, p = v$p_full, estimate = v$est)
}, notes = "it48: eperm_c4 with a minimum of two leave-one-out discoveries before an option is used")
register_candidate("eperm_c5q", function(counts, meta, formula, tested_term) {
  v <- .epg_fit(counts, meta, formula, tested_term, paths = "cont", adet = TRUE, nodet_conf = TRUE, pk = "pos", minsw = 2L)
  data.frame(feature = v$feature, p = v$p_full, q = .storey_gq(v$p_full), estimate = v$est)
}, notes = "it48: eperm_c5 with Storey's adaptive BH, applied only when plain BH rejects")
# it48b: no winsorising. it36 capped each taxon's top 3% of proportions (ZicoSeq's default), which helped msq under the
# asymptotic t reference (srv13 +1.6 TP) by taming skew. Under exact permutation skew costs no validity, and the cap
# only removes signal: banked msq R00 / R11 hmp_tongue cells, eperm_cs1 selection, winsor 0.03 -> 0: 8.8 -> 9.4 and
# 6.6 -> 7.0 TP per cell, sd2 R08 3.4 -> 3.6, false positives 0.4-0.6 either way. (It does not explain why our linear
# kernel stays far below LDM's frequency-scale test: 5.8 vs 9.6 on msq R00.)
register_candidate("eperm_c5w", function(counts, meta, formula, tested_term) {
  v <- .epg_fit(counts, meta, formula, tested_term, paths = "cont", adet = TRUE, nodet_conf = TRUE, pk = "pos", minsw = 2L, winsor = 0)
  data.frame(feature = v$feature, p = v$p_full, q = .storey_gq(v$p_full), estimate = v$est)
}, notes = "it48b: eperm_c5q without winsorising")
# it48c: no size factors. With winsorising off, the same banked cells without size factors (balanced-depth designs
# are the only ones that use them): msq R00 hmp_tongue p_full 9.4 -> 10.2 TP per cell (LDM's own omnibus 10.0), the
# linear kernel 5.8 -> 8.6 (LDM frequency scale 9.6), sqrt 8.6 -> 9.2; sd2 R08 p_full 3.6 -> 2.8 (pos-dominated in
# c4 anyway). Recentring off changed nothing (lin 5.6, sqrt 9.2). Size factors came in with it35 for the asymptotic
# path; whether they still pay under the permutation path is for the dev suite to say. The switch is an option read
# by a wrapper around .size_factors, so every path (permutation, continuous, efull) drops them together; the memo
# caches are cleared around the call because their keys do not include the option.
if (!exists(".size_factors_it35")) .size_factors_it35 <- .size_factors
.size_factors <- function(counts, depth, ...) if (isTRUE(getOption("pursue.nosf"))) rep(0, ncol(counts)) else .size_factors_it35(counts, depth, ...)
.clear_fit_memos <- function() for (e in c(".efull_memo", ".eperm_memo", ".epg_memo", ".eps_memo")) if (exists(e)) assign("key", NULL, envir = get(e))
register_candidate("eperm_c5ws", function(counts, meta, formula, tested_term) {
  old <- options(pursue.nosf = TRUE); .clear_fit_memos(); on.exit({ options(old); .clear_fit_memos() })
  v <- .epg_fit(counts, meta, formula, tested_term, paths = "cont", adet = TRUE, nodet_conf = TRUE, pk = "pos", minsw = 2L, winsor = 0)
  data.frame(feature = v$feature, p = v$p_full, q = .storey_gq(v$p_full), estimate = v$est)
}, notes = "it48c: eperm_c5w without size factors")
