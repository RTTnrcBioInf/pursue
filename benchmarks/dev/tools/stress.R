# tools/stress.R -- stress transforms applied to simulated cells (dev suite "ext" settings)
#
# apply_bloom(sim, fold): a hard compositional shock. A new taxon appears in every case sample and
# takes (fold-1)/fold of its reads; every other taxon keeps its ABSOLUTE abundance, so its relative
# abundance in cases falls by `fold`. Sequencing depth does not grow with a bloom, so each case
# sample is then resampled (multinomially) back to its original depth. Absolute truth: the bloom
# taxon is DA; every other taxon keeps its label. A relative-scale test calls most taxa "down";
# PURSUE 0.2's centring does not survive fold 4 either (FPR 0.157 on hmp_tongue).
apply_bloom <- function(sim, fold = 4, seed = 1L) {
  set.seed(seed); ct <- sim$counts; g <- as.integer(sim$meta[[sim$tested_term]])
  g <- if (is.factor(sim$meta[[sim$tested_term]])) g - 1L else as.integer(g > 0)
  N <- colSums(ct); bl <- round(N * (fold - 1)) * g
  ct2 <- rbind(ct, BLOOM = bl)
  for (i in which(g == 1L)) ct2[, i] <- stats::rmultinom(1L, N[i], ct2[, i] / sum(ct2[, i]))
  tr <- sim$truth; add <- tr[1, , drop = FALSE]; add[1, ] <- NA; add$feature <- "BLOOM"; add$truth_abs <- 1L
  if ("truth_type" %in% names(add)) add$truth_type <- "bloom"
  sim$counts <- ct2; sim$truth <- rbind(tr, add); sim$meta$depth <- colSums(ct2); sim
}
