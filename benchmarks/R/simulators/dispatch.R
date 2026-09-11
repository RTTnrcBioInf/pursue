# -----------------------------------------------------------------------------
# dispatch.R -- one entry point for every generator, plus the external simulator
# wrappers. Every generator returns list(counts [features x samples], meta, truth,
# tested_term, formula). Regime factors a simulator cannot express are recorded in
# `unsupported` and the cell is skipped by run_cell (status "regime_unsupported").
#
# External wrappers are [UNVERIFIED] in the development container (no CRAN /
# Bioconductor access); hpc/smoke_test_simulators.R exercises them on the cluster.
# -----------------------------------------------------------------------------

simulator_registry <- function() data.frame(stringsAsFactors = FALSE,
  id      = c("house", "msq", "mid", "sd2", "sps"),
  label   = c("in-house ZI lognormal-Poisson", "GUniFrac::SimulateMSeq", "MIDASim (parametric)", "sparseDOSSA 2", "SPsimSeq"),
  package = c("", "GUniFrac", "MIDASim", "SparseDOSSA2", "SPsimSeq"),
  truth_scale = c("absolute", "absolute", "relative", "absolute", "absolute"))

simulator_available <- function(id) { pk <- simulator_registry()$package[simulator_registry()$id == id]; !nzchar(pk) || requireNamespace(pk, quietly = TRUE) }

simulate_cell_data <- function(simulator, template, regime, seed, cache_dir = NULL) {
  switch(simulator,
    house = simulate_house(template, regime, seed),
    msq   = simulate_msq(template, regime, seed),
    mid   = simulate_mid(template, regime, seed, cache_dir),
    sd2   = simulate_sd2(template, regime, seed, cache_dir),
    sps   = simulate_sps(template, regime, seed),
    stop("unknown simulator: ", simulator))
}

.unsupported <- function(regime, ok_factors) {
  bad <- c(); ref <- reference_regime()
  for (f in names(ref)) if (!(f %in% ok_factors) && !identical(as.character(regime[[f]]), as.character(ref[[f]]))) bad <- c(bad, f)
  bad
}
.effects_for <- function(regime) { es <- effect_sizes(regime$effect); list(lfc2 = if (es$graded) mean(es$lfc2) else es$lfc2, graded = es$graded) }

# ---- GUniFrac::SimulateMSeq [verified API from the original engine] ----
# Supports: n, m, da_frac, balance (balanced/unbalanced), effect (covariate.eff.mean), depth
# confounding (depth.conf.factor = log(c)), covariate confounding (confounder.type =
# "continuous", conf.cov.cor = phi), exposure continuous (covariate.type = "continuous").
# Not supported: prevalence-only signals, repeated measures, bloom, hetero.
simulate_msq <- function(template, regime, seed) {
  un <- .unsupported(regime, c("n_per_group", "m", "da_frac", "balance", "effect", "depth_conf", "conf_phi", "exposure", "signal_type"))
  if (identical(regime$signal_type, "prevalence")) un <- c(un, "signal_type")
  if (length(un)) return(list(unsupported = un))
  set.seed(seed)
  ref <- template$counts; storage.mode(ref) <- "numeric"
  paras <- template$msq_paras; if (is.null(paras)) paras <- GUniFrac:::EstPara(ref)
  es <- .effects_for(regime); ef <- es$lfc2 * log(2)
  sim <- GUniFrac::SimulateMSeq(ref.otu.tab = ref, model.paras = paras, nSam = 2L * regime$n_per_group, nOTU = regime$m,
    diff.otu.pct = regime$da_frac, diff.otu.direct = if (regime$balance == "balanced") "balanced" else "unbalanced",
    diff.otu.mode = "mix", covariate.type = if (regime$exposure == "continuous") "continuous" else "binary", grp.ratio = 1,
    covariate.eff.mean = ef, covariate.eff.sd = if (es$graded) ef / 2 else 0, error.sd = 0,
    depth.mu = 10000, depth.theta = 5, depth.conf.factor = log(regime$depth_conf),
    confounder.type = if (regime$conf_phi > 0) "continuous" else "none", conf.cov.cor = regime$conf_phi,
    conf.diff.otu.pct = 0, conf.nondiff.otu.pct = if (regime$conf_phi > 0) 0.1 else 0, confounder.eff.mean = if (regime$conf_phi > 0) log(2) else 0, confounder.eff.sd = 0)
  X <- sim$otu.tab.sim; x <- as.numeric(sim$covariate); if (is.matrix(x)) x <- x[, 1]
  da <- rep(FALSE, nrow(X)); da[sim$diff.otu.ind] <- TRUE
  N <- colSums(X); R <- sweep(X, 2L, N, "/")
  binary <- regime$exposure != "continuous"
  rel <- if (binary) log2((rowMeans(R[, x == 1, drop = FALSE]) + 1e-12) / (rowMeans(R[, x == 0, drop = FALSE]) + 1e-12)) else NA
  truth <- data.frame(feature = rownames(X), truth_type = ifelse(da, "abundance", "none"), truth_abs = as.integer(da),
                      truth_lfc2 = ifelse(da, NA, 0), truth_rel_lfc2 = rel, truth_rel = as.integer(if (binary) abs(rel) > 0.25 else da), stringsAsFactors = FALSE)
  meta <- data.frame(row.names = colnames(X), group = if (binary) factor(ifelse(x == 1, "case", "control"), levels = c("control", "case")) else NA,
                     exposure = x, confounder = if (!is.null(sim$confounder)) as.numeric(sim$confounder)[seq_len(ncol(X))] else stats::rnorm(ncol(X)), depth = N)
  list(counts = X, meta = meta, truth = truth, tested_term = if (binary) "group" else "exposure",
       formula = stats::as.formula(paste("~", if (binary) "group" else "exposure", if (regime$conf_phi > 0) "+ confounder" else "")))
}

# ---- MIDASim [UNVERIFIED] -- parametric mode; effects via MIDASim.modify on rel.abund / lib.size ----
simulate_mid <- function(template, regime, seed, cache_dir = NULL) {
  un <- .unsupported(regime, c("n_per_group", "m", "da_frac", "balance", "effect", "depth_conf", "signal_type"))
  if (identical(regime$signal_type, "prevalence")) un <- c(un, "signal_type")
  if (length(un)) return(list(unsupported = un))
  set.seed(seed)
  ct <- t(template$counts); m_avail <- ncol(ct)
  keep <- order(colMeans(ct > 0), decreasing = TRUE)[seq_len(min(regime$m, m_avail))]
  ct <- ct[, keep, drop = FALSE]
  setup <- template$mid_setup
  if (is.null(setup)) { f <- if (!is.null(cache_dir)) file.path(cache_dir, paste0("midasim_", template$id, "_", regime$m, ".rds")) else NULL
    if (!is.null(f) && file.exists(f)) setup <- readRDS(f) else { setup <- MIDASim::MIDASim.setup(ct, mode = "parametric", n.break.ties = 10); if (!is.null(f)) saveRDS(setup, f) } }
  n_per <- regime$n_per_group; n <- 2L * n_per; m <- ncol(ct)
  es <- .effects_for(regime); n_da <- round(regime$da_frac * m); da <- sample.int(m, n_da)
  p_up <- switch(regime$balance, balanced = 0.5, "80_20" = 0.8, "100_0" = 1.0)
  fold <- rep(1, m); if (n_da > 0) fold[da] <- 2^(es$lfc2 * ifelse(stats::runif(n_da) < p_up, 1, -1))
  base_rel <- setup$mean.rel.abund
  gen_group <- function(ngrp, mult, depth_mult) {
    rel <- base_rel * mult; rel <- rel / sum(rel)
    mod <- MIDASim::MIDASim.modify(setup, lib.size = round(sample(colSums(template$counts), ngrp, replace = TRUE) * depth_mult), mean.rel.abund = rel)
    MIDASim::MIDASim(mod)$sim_count }
  X0 <- gen_group(n_per, rep(1, m), 1); X1 <- gen_group(n_per, fold, regime$depth_conf)
  X <- rbind(X0, X1); colnames(X) <- colnames(ct); rownames(X) <- paste0("S", seq_len(n)); X <- t(X)
  x <- rep(c(0L, 1L), each = n_per); N <- colSums(X); R <- sweep(X, 2L, N, "/")
  rel <- log2((rowMeans(R[, x == 1, drop = FALSE]) + 1e-12) / (rowMeans(R[, x == 0, drop = FALSE]) + 1e-12))
  isda <- fold != 1
  truth <- data.frame(feature = rownames(X), truth_type = ifelse(isda, "abundance", "none"), truth_abs = as.integer(isda),
                      truth_lfc2 = log2(fold), truth_rel_lfc2 = rel, truth_rel = as.integer(abs(rel) > 0.25), stringsAsFactors = FALSE)
  meta <- data.frame(row.names = colnames(X), group = factor(ifelse(x == 1, "case", "control"), levels = c("control", "case")), exposure = x, depth = N)
  list(counts = X, meta = meta, truth = truth, tested_term = "group", formula = ~ group)
}

# ---- sparseDOSSA 2 [UNVERIFIED] -- fit_SparseDOSSA2 on the template (cached), spike-in via metadata_effects ----
simulate_sd2 <- function(template, regime, seed, cache_dir = NULL) {
  un <- .unsupported(regime, c("n_per_group", "m", "da_frac", "balance", "effect", "signal_type", "conf_phi"))
  if (identical(regime$signal_type, "prevalence")) un <- c(un, "signal_type")
  if (length(un)) return(list(unsupported = un))
  set.seed(seed)
  fit <- template$sd2_fit
  if (is.null(fit)) { f <- if (!is.null(cache_dir)) file.path(cache_dir, paste0("sd2_", template$id, ".rds")) else NULL
    if (!is.null(f) && file.exists(f)) fit <- readRDS(f) else { fit <- SparseDOSSA2::fit_SparseDOSSA2(data = template$counts, control = list(verbose = FALSE)); if (!is.null(f)) saveRDS(fit, f) } }
  n_per <- regime$n_per_group; n <- 2L * n_per; m <- min(regime$m, nrow(template$counts))
  es <- .effects_for(regime); n_da <- round(regime$da_frac * m); da <- sample.int(m, n_da)
  p_up <- switch(regime$balance, balanced = 0.5, "80_20" = 0.8, "100_0" = 1.0)
  x <- rep(c(0, 1), each = n_per); conf <- if (regime$conf_phi > 0) regime$conf_phi * as.numeric(scale(x)) + sqrt(1 - regime$conf_phi^2) * stats::rnorm(n) else stats::rnorm(n)
  md <- cbind(group = x, confounder = conf)
  spike <- if (n_da > 0) data.frame(metadata_datum = 1, feature_spiked = rownames(template$counts)[da], associated_property = "abundance",
                                    effect_size = es$lfc2 * log(2) * ifelse(stats::runif(n_da) < p_up, 1, -1)) else NULL
  sim <- SparseDOSSA2::SparseDOSSA2(template = fit, n_sample = n, n_feature = m, spike_metadata = if (is.null(spike)) "none" else spike,
                                    metadata_matrix = md, verbose = FALSE)
  X <- sim$simulated_data; N <- colSums(X); R <- sweep(X, 2L, N, "/")
  isda <- rownames(X) %in% (if (n_da > 0) spike$feature_spiked else character(0))
  rel <- log2((rowMeans(R[, x == 1, drop = FALSE]) + 1e-12) / (rowMeans(R[, x == 0, drop = FALSE]) + 1e-12))
  truth <- data.frame(feature = rownames(X), truth_type = ifelse(isda, "abundance", "none"), truth_abs = as.integer(isda),
                      truth_lfc2 = ifelse(isda, NA, 0), truth_rel_lfc2 = rel, truth_rel = as.integer(abs(rel) > 0.25), stringsAsFactors = FALSE)
  meta <- data.frame(row.names = colnames(X), group = factor(ifelse(x == 1, "case", "control"), levels = c("control", "case")), exposure = x, confounder = conf, depth = N)
  list(counts = X, meta = meta, truth = truth, tested_term = "group", formula = if (regime$conf_phi > 0) ~ group + confounder else ~ group)
}

# ---- SPsimSeq [UNVERIFIED] -- semi-parametric resampling with SpiecEasi-style correlation ----
simulate_sps <- function(template, regime, seed) {
  un <- .unsupported(regime, c("n_per_group", "m", "da_frac", "effect", "signal_type"))
  if (identical(regime$signal_type, "prevalence")) un <- c(un, "signal_type")
  if (length(un)) return(list(unsupported = un))
  set.seed(seed)
  ct <- template$counts; m <- min(regime$m, nrow(ct)); n_per <- regime$n_per_group
  g <- template$meta[[template$group_var]]
  grp <- if (!is.na(template$group_var) && !is.null(g)) as.integer(factor(g)) else sample(rep(1:2, length.out = ncol(ct)))
  es <- .effects_for(regime)
  sim <- SPsimSeq::SPsimSeq(n.sim = 1, s.data = ct, group = grp, n.genes = m, batch.config = 1, group.config = c(0.5, 0.5),
                            tot.samples = 2L * n_per, pDE = regime$da_frac, lfc.thrld = es$lfc2, t.thrld = 2.5, llStat.thrld = 5,
                            model.zero.prob = TRUE, genewiseCor = TRUE, result.format = "list", return.details = TRUE, verbose = FALSE)
  d <- sim$sim.data.list[[1]]; X <- as.matrix(d$counts); md <- d$colData; rd <- d$rowData
  x <- as.integer(md$Group == levels(factor(md$Group))[2]); N <- colSums(X); R <- sweep(X, 2L, N, "/")
  isda <- as.logical(rd$DE.ind)
  rel <- log2((rowMeans(R[, x == 1, drop = FALSE]) + 1e-12) / (rowMeans(R[, x == 0, drop = FALSE]) + 1e-12))
  truth <- data.frame(feature = rownames(X), truth_type = ifelse(isda, "abundance", "none"), truth_abs = as.integer(isda),
                      truth_lfc2 = ifelse(isda, NA, 0), truth_rel_lfc2 = rel, truth_rel = as.integer(abs(rel) > 0.25), stringsAsFactors = FALSE)
  meta <- data.frame(row.names = colnames(X), group = factor(ifelse(x == 1, "case", "control"), levels = c("control", "case")), exposure = x, depth = N)
  list(counts = X, meta = meta, truth = truth, tested_term = "group", formula = ~ group)
}
