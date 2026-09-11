# -----------------------------------------------------------------------------
# nullperm.R -- Axis C: null calibration on real data by label permutation.
# scheme in {full, stratified, depth_quintile, cluster}. The template must have a
# binary grouping variable (template$group_var) or one is created by a random split.
# Returns the real counts with a permuted binary `group` and, for `stratified` /
# `cluster`, the stratum/cluster variable in meta.
# -----------------------------------------------------------------------------

null_permute <- function(template, scheme = c("full", "stratified", "depth_quintile", "cluster"),
                         seed = 1L, strat_var = NULL, cluster_var = NULL, n_max = 200L) {
  scheme <- match.arg(scheme); set.seed(seed)
  ct <- template$counts; meta <- template$meta
  if (ncol(ct) > n_max) { s <- sort(sample.int(ncol(ct), n_max)); ct <- ct[, s, drop = FALSE]; meta <- meta[s, , drop = FALSE] }
  n <- ncol(ct)
  g <- if (!is.na(template$group_var) && template$group_var %in% colnames(meta)) {
    v <- meta[[template$group_var]]; as.integer(v == sort(unique(v))[1]) } else stats::rbinom(n, 1L, 0.5)
  perm <- switch(scheme,
    full = sample(g),
    stratified = { sv <- if (!is.null(strat_var) && strat_var %in% colnames(meta)) meta[[strat_var]] else sample(rep(1:2, length.out = n))
                   ave(g, sv, FUN = sample) },
    depth_quintile = { q <- cut(colSums(ct), stats::quantile(colSums(ct), 0:5 / 5), include.lowest = TRUE); ave(g, q, FUN = sample) },
    cluster = { cv <- if (!is.null(cluster_var) && cluster_var %in% colnames(meta)) meta[[cluster_var]] else sample(rep(seq_len(ceiling(n / 4)), length.out = n))
                lab <- tapply(g, cv, function(z) z[1]); lab <- lab[sample(length(lab))]; names(lab) <- unique(cv); as.integer(lab[as.character(cv)]) })
  meta_out <- data.frame(row.names = colnames(ct), group = factor(ifelse(perm == 1L, "case", "control"), levels = c("control", "case")),
                         depth = colSums(ct), stringsAsFactors = FALSE)
  if (scheme == "stratified") meta_out$stratum <- factor(if (!is.null(strat_var) && strat_var %in% colnames(meta)) meta[[strat_var]] else sample(rep(1:2, length.out = n)))
  truth <- data.frame(feature = rownames(ct), truth_type = "none", truth_abs = 0L, truth_rel = 0L, stringsAsFactors = FALSE)
  list(counts = ct, meta = meta_out, truth = truth, tested_term = "group",
       formula = if (scheme == "stratified") ~ group + stratum else ~ group, scheme = scheme)
}
