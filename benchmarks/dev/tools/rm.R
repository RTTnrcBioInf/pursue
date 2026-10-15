# -----------------------------------------------------------------------------
# rm.R -- repeated-measures cells on real data (2026-10-15). Axis C's "cluster" scheme uses random groups of four that
# carry no real correlation and are not passed to the methods; R23 exists only on `house`. These cells keep the real
# within-subject structure: every sample is a real sample, two visits per subject.
#
#   P_subj  exposure fixed per subject (half the subjects are cases; both of a subject's visits share the label)
#   P_pair  exposure varies within subject (one of each subject's two visits, chosen at random, is the exposed one)
#   *.null  the same with no implanted taxa: every call is false, and the label assignment makes the null exact
#
# Effects are implanted exactly as implant() does it (benchmarks/R/simulators/implant.R: abundance by count scaling,
# prevalence by moving non-zero cells, 1/3 each of abundance / prevalence / both, medium effects), then every sample is
# redrawn multinomially at its own depth. meta carries `subject` (what c12 and the subject-aware comparators use) and
# `visit`. The formula is ~ group; how a method uses `subject` is its own business.
#
# Template: a dev template with subject labels attached by rm_template() from a sample map (SRS id -> RSID, VISITNO;
# benchmarks/devdata/hmp16s_subjects.csv, written on the server by dev/tools/hmp_subjects.R from HMP16SData). Technical
# replicates (one SRS sequenced twice) keep the deeper run; a subject's two earliest distinct visits are used.
# -----------------------------------------------------------------------------

RM_SPECS <- list(
  P_subj = list(design = "subject", n_subj = 40L, da_frac = 0.10, effect = "medium", signal_type = "mixed", balance = "balanced"),
  P_pair = list(design = "paired",  n_subj = 40L, da_frac = 0.10, effect = "medium", signal_type = "mixed", balance = "balanced"))
RM_SPECS[["P_subj.null"]] <- utils::modifyList(RM_SPECS$P_subj, list(da_frac = 0))
RM_SPECS[["P_pair.null"]] <- utils::modifyList(RM_SPECS$P_pair, list(da_frac = 0))

rm_template <- function(tpl, map) {
  srs <- sub("^[^.]+[.]([^.]+)[.].*$", "\\1", tpl$meta$X.SampleID)
  i <- match(srs, map$SRS_SAMPLE_ID); ok <- !is.na(i)
  ct <- tpl$counts[, ok, drop = FALSE]; srs <- srs[ok]; mp <- map[i[ok], , drop = FALSE]
  N <- colSums(ct); o <- order(srs, -N); keep <- o[!duplicated(srs[o])]                     # technical replicates: deeper run
  ct <- ct[, keep, drop = FALSE]; mp <- mp[keep, , drop = FALSE]
  key <- paste(mp$RSID, mp$VISITNO); N <- colSums(ct); o <- order(key, -N); keep <- o[!duplicated(key[o])]
  ct <- ct[, keep, drop = FALSE]; mp <- mp[keep, , drop = FALSE]
  meta <- data.frame(row.names = colnames(ct), subject = as.character(mp$RSID), visit = as.integer(mp$VISITNO),
                     sex = mp$SEX, run_center = mp$RUN_CENTER, stringsAsFactors = FALSE)
  ct <- ct[rowSums(ct > 0) > 0, , drop = FALSE]
  list(counts = ct, meta = meta, id = paste0(tpl$id, "_long"),
       n_matched = sum(ok), n_samples = ncol(ct), n_subjects = length(unique(meta$subject)),
       n_subj2 = sum(table(meta$subject) >= 2L))
}

# implant()'s effect block, applied to a given exposure vector x (0/1)
.rm_effects <- function(X, x, spec) {
  m <- nrow(X); es <- effect_sizes(spec$effect); n_da <- round(spec$da_frac * m)
  typ <- rep("none", m); lfc2 <- numeric(m); shift <- numeric(m); Xf <- X * 1.0
  if (n_da > 0) {
    da <- sample.int(m, n_da)
    kind <- switch(spec$signal_type, abundance = rep("abundance", n_da), prevalence = rep("prevalence", n_da),
                   mixed = sample(c("abundance", "prevalence", "both"), n_da, TRUE))
    p_up <- switch(spec$balance, balanced = 0.5, "80_20" = 0.8, "100_0" = 1.0)
    sgn <- ifelse(stats::runif(n_da) < p_up, 1, -1)
    mag_a <- if (es$graded) stats::runif(n_da, es$lfc2[1], es$lfc2[2]) else rep(es$lfc2, n_da)
    mag_p <- if (es$graded) stats::runif(n_da, 0.2, 0.5) else rep(switch(spec$effect, small = 0.2, medium = 0.35, large = 0.5), n_da)
    typ[da] <- kind
    for (k in seq_len(n_da)) {
      j <- da[k]
      if (kind[k] %in% c("abundance", "both")) { lfc2[j] <- sgn[k] * mag_a[k]; Xf[j, ] <- Xf[j, ] * ifelse(x == 1L, 2^lfc2[j], 1) }
      if (kind[k] %in% c("prevalence", "both")) {
        shift[j] <- sgn[k] * mag_p[k]; up <- if (sgn[k] > 0) x == 1L else x == 0L
        nz_down <- which(!up & Xf[j, ] > 0); z_up <- which(up & Xf[j, ] == 0)
        k_move <- min(length(nz_down), length(z_up), round(mag_p[k] * length(nz_down)))
        if (k_move > 0) { med_nz <- stats::median(Xf[j, Xf[j, ] > 0])
          Xf[j, nz_down[sample.int(length(nz_down), k_move)]] <- 0; Xf[j, z_up[sample.int(length(z_up), k_move)]] <- med_nz }
      }
    }
  }
  list(Xf = Xf, typ = typ, lfc2 = lfc2, shift = shift)
}

rm_cell <- function(tpl, spec, seed = 1L) {
  set.seed(seed)
  v <- tpl$meta; subj2 <- names(which(table(v$subject) >= 2L))
  if (length(subj2) < spec$n_subj) return(list(unsupported = sprintf("subjects_with_2_visits[%d<%d]", length(subj2), spec$n_subj)))
  S <- sort(sample(subj2, spec$n_subj))
  idx <- unlist(lapply(S, function(s) { w <- which(v$subject == s); w[order(v$visit[w])][1:2] }))
  X <- tpl$counts[, idx, drop = FALSE]; sj <- v$subject[idx]; vis <- v$visit[idx]; n <- length(idx)
  x <- if (spec$design == "subject") { cs <- sample(S, spec$n_subj %/% 2L); as.integer(sj %in% cs) } else
    as.integer(unlist(lapply(seq_along(S), function(k) sample(c(0L, 1L)))))                    # pairs are adjacent: 2k-1, 2k
  ef <- .rm_effects(X, x, spec)
  N <- colSums(X); m <- nrow(X)
  Xn <- sapply(seq_len(n), function(i) { p <- ef$Xf[, i]; if (sum(p) <= 0) rep(0L, m) else stats::rmultinom(1L, N[i], p / sum(p))[, 1] })
  dimnames(Xn) <- list(rownames(X), paste0("S", seq_len(n))); storage.mode(Xn) <- "integer"
  truth <- data.frame(feature = rownames(X), truth_type = ef$typ, truth_abs = as.integer(ef$typ != "none"), truth_lfc2 = ef$lfc2,
                      truth_prev_shift = ef$shift, truth_rel = as.integer(ef$typ != "none"), stringsAsFactors = FALSE)
  meta <- data.frame(row.names = colnames(Xn), group = factor(ifelse(x == 1L, "case", "control"), levels = c("control", "case")),
                     exposure = x, subject = factor(sj), visit = vis, depth = N, stringsAsFactors = FALSE)
  list(counts = Xn, meta = meta, truth = truth, tested_term = "group", formula = ~ group, design = spec$design)
}
