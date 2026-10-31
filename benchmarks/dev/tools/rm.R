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

# srv32 stress templates for the within-subject path (development data only). Both are fixed before any exposure is
# drawn, so a random within-subject exposure stays exchangeable and an exact test must hold on them.
#  rm_mixed: one sample from site A and one from site B per subject (earliest visit at each; HMP RSIDs shared), taxa
#            unioned (absent = 0) -- pairs whose two members differ systematically and strongly.
#  rm_thinned: one of each subject's first two samples (chosen at random, fixed seed) binomially thinned to `frac` of its
#            depth -- pairs whose two members differ in depth.
rm_mixed <- function(tA, tB, id = "mixed") {
  pick1 <- function(t) { m <- t$meta; o <- order(m$subject, m$visit); k <- o[!duplicated(m$subject[o])]; list(ct = t$counts[, k, drop = FALSE], s = m$subject[k]) }
  a <- pick1(tA); b <- pick1(tB); S <- intersect(a$s, b$s)
  ia <- match(S, a$s); ib <- match(S, b$s); f <- union(rownames(a$ct), rownames(b$ct))
  fill <- function(ct, k) { z <- matrix(0L, length(f), length(k), dimnames = list(f, colnames(ct)[k])); z[rownames(ct), ] <- ct[, k, drop = FALSE]; z }
  ct <- cbind(fill(a$ct, ia), fill(b$ct, ib)); ct <- ct[rowSums(ct > 0) > 0, , drop = FALSE]
  meta <- data.frame(row.names = colnames(ct), subject = c(S, S), visit = rep(1:2, each = length(S)), stringsAsFactors = FALSE)
  list(counts = ct, meta = meta, id = id, n_matched = ncol(ct), n_samples = ncol(ct), n_subjects = length(S), n_subj2 = length(S))
}
rm_thinned <- function(t, frac = 0.1, seed = 3201L, id = "thinned") {
  m <- t$meta; ct <- t$counts; s2 <- names(which(table(m$subject) >= 2L))
  .with_seed <- function(sd, e) { old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()); set.seed(sd)
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv())); e }
  .with_seed(seed, for (s in s2) { w <- which(m$subject == s); w <- w[order(m$visit[w])][1:2]; j <- w[sample.int(2L, 1L)]
    ct[, j] <- stats::rbinom(nrow(ct), ct[, j], frac) })
  storage.mode(ct) <- "integer"; t$counts <- ct; t$id <- id; t
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
  rnd <- identical(tpl$pick, "random")                                     # conf2 baboons: two samples per animal at random
  idx <- unlist(lapply(S, function(s) { w <- which(v$subject == s); w <- w[order(v$visit[w])]; if (rnd) w[sort(sample.int(length(w), 2L))] else w[1:2] }))
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

# ---- house with a within-subject exposure and more than two visits (2026-10-15, it57 check) -------------------------
# A code check for c13's paired path where the development data cannot reach (hmp_tongue has two subjects with a third
# visit): simulate_house()'s generator (R/simulators/sim_house.R, unchanged) with the subject structure of R23 -- a
# per-taxon subject intercept on log abundance, ICC 0.3 -- but the exposure assigned WITHIN subject. Copied here so the
# benchmark generator stays as it is. Per subject: `vis` visits (or a random number in vis_range), and for a fraction
# `mixed` of the subjects a random 1 .. v-1 of them exposed; the other subjects are all unexposed (no within-subject
# contrast: c13 must leave them out). Truth as in house: abundance / prevalence / both, absolute scale.
RM_SPECS[["HW5"]] <- list(gen = "house", design = "within", vis = 5L, vis_range = NULL, mixed = 1, n_per_group = 50L, m = 500L,
                          da_frac = 0.10, effect = "medium", signal_type = "mixed", balance = "balanced")
RM_SPECS[["HWv"]] <- utils::modifyList(RM_SPECS$HW5, list(vis_range = 2:8, mixed = 0.8))
RM_SPECS[["HW5.null"]] <- utils::modifyList(RM_SPECS$HW5, list(da_frac = 0))
RM_SPECS[["HWv.null"]] <- utils::modifyList(RM_SPECS$HWv, list(da_frac = 0))

simulate_house_within <- function(template, spec, seed = 1L) {
  set.seed(seed)
  prof <- template_profile(template); m <- spec$m; n_subj <- max(6L, round(2L * spec$n_per_group / 5L))
  v <- if (is.null(spec$vis_range)) rep(spec$vis, n_subj) else spec$vis_range[sample.int(length(spec$vis_range), n_subj, TRUE)]
  subject <- rep(seq_len(n_subj), v); n <- length(subject)
  mix <- seq_len(n_subj) %in% sample.int(n_subj, round(spec$mixed * n_subj))
  x <- integer(n)
  for (s in which(mix)) { ii <- which(subject == s); x[ii[sample.int(length(ii), sample.int(length(ii) - 1L, 1L))]] <- 1L }
  N <- pmax(sample(prof$depth, n, replace = TRUE), 200L)
  idx <- sample(seq_len(prof$n_features), m, replace = m > prof$n_features)
  psi <- pmin(pmax(prof$prevalence[idx], 0.05), 0.995)
  mu  <- prof$mean_log_relab[idx]; mu[!is.finite(mu)] <- stats::median(mu, na.rm = TRUE)
  sg  <- prof$sd_log_relab[idx]; sg[!is.finite(sg) | sg < 0.3] <- stats::median(sg[is.finite(sg) & sg >= 0.3], na.rm = TRUE); sg <- pmin(sg, 2.5)
  es <- effect_sizes(spec$effect); n_da <- round(spec$da_frac * m)
  typ <- rep("none", m); lfc2 <- numeric(m); logor <- numeric(m)
  if (n_da > 0) {
    da <- sample.int(m, n_da)
    kind <- switch(spec$signal_type, abundance = rep("abundance", n_da), prevalence = rep("prevalence", n_da),
                   mixed = sample(c("abundance", "prevalence", "both"), n_da, replace = TRUE))
    p_up <- switch(spec$balance, balanced = 0.5, "80_20" = 0.8, "100_0" = 1.0)
    sgn <- ifelse(stats::runif(n_da) < p_up, 1, -1)
    typ[da] <- kind
    lfc2[da] <- ifelse(kind %in% c("abundance", "both"), sgn * es$lfc2, 0)
    logor[da] <- ifelse(kind %in% c("prevalence", "both"), sgn * es$logor, 0)
  }
  beta <- lfc2 * log(2); icc <- 0.3
  A <- matrix(0, n, m)
  for (j in seq_len(m)) {
    S <- stats::rbinom(n, 1L, stats::plogis(stats::qlogis(psi[j]) + logor[j] * x))
    u_subj <- stats::rnorm(n_subj, 0, sqrt(icc / (1 - icc)) * sg[j])[subject]
    A[, j] <- exp(mu[j] + beta[j] * x + u_subj + stats::rnorm(n, 0, sg[j])) * S
  }
  tot <- rowSums(A); tot[tot <= 0] <- 1; R <- A / tot
  X <- matrix(0L, n, m); for (i in seq_len(n)) if (sum(R[i, ]) > 0) X[i, ] <- stats::rpois(m, N[i] * R[i, ])
  feat <- sprintf("F%04d", seq_len(m)); colnames(X) <- feat; rownames(X) <- sprintf("S%04d", seq_len(n))
  truth <- data.frame(feature = feat, truth_type = typ, truth_abs = as.integer(typ != "none"), truth_lfc2 = lfc2, truth_logor = logor,
                      stringsAsFactors = FALSE)
  meta <- data.frame(row.names = rownames(X), group = factor(ifelse(x == 1L, "case", "control"), levels = c("control", "case")),
                     exposure = x, subject = factor(subject), depth = N, stringsAsFactors = FALSE)
  list(counts = t(X), meta = meta, truth = truth, tested_term = "group", formula = ~ group, design = spec$design)
}
