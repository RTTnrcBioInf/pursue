# -----------------------------------------------------------------------------
# realism.R -- the realism gate (protocol section 4.4) for a (simulator, template)
# pair: a classifier's ability to tell simulated samples from real ones, plus
# marginal and correlation-structure distances. Uses ranger or randomForest when
# installed, otherwise a regularised logistic regression on principal components.
# -----------------------------------------------------------------------------

realism_gate <- function(template, simulator, n_sim = 200L, seed = 1L, cache_dir = NULL) {
  set.seed(seed)
  reg <- reference_regime(); reg$da_frac <- 0; reg$n_per_group <- as.integer(ceiling(n_sim / 2)); reg$m <- nrow(template$counts)
  sim <- simulate_cell_data(simulator, template, reg, seed, cache_dir)
  if (!is.null(sim$unsupported)) return(data.frame(simulator = simulator, template = template$id, auroc = NA, ks_median = NA, frob = NA, note = "unsupported"))
  real <- template$counts; n_real <- min(ncol(real), n_sim); real <- real[, sample.int(ncol(real), n_real), drop = FALSE]
  common <- intersect(rownames(real), rownames(sim$counts))
  if (length(common) < 20L) { # simulators that rename features: match by rank of mean abundance
    o_r <- order(rowMeans(real), decreasing = TRUE); o_s <- order(rowMeans(sim$counts), decreasing = TRUE)
    k <- min(nrow(real), nrow(sim$counts)); real <- real[o_r[seq_len(k)], ]; simc <- sim$counts[o_s[seq_len(k)], ]; rownames(simc) <- rownames(real)
  } else { real <- real[common, ]; simc <- sim$counts[common, ] }
  feat <- function(ct) { N <- colSums(ct); R <- sweep(ct, 2L, N, "/")
    cbind(t(log(R + 1e-6)), richness = colSums(ct > 0), log_depth = log(N),
          shannon = apply(R, 2L, function(p) { p <- p[p > 0]; -sum(p * log(p)) })) }
  Fr <- feat(real); Fs <- feat(simc); Z <- rbind(Fr, Fs); y <- factor(c(rep("real", nrow(Fr)), rep("sim", nrow(Fs))))
  auc <- tryCatch({
    if (requireNamespace("ranger", quietly = TRUE)) { fit <- ranger::ranger(x = Z, y = y, probability = TRUE, num.trees = 500, seed = seed)
      pr <- fit$predictions[, "sim"] } else if (requireNamespace("randomForest", quietly = TRUE)) {
      fit <- randomForest::randomForest(Z, y, ntree = 500); pr <- fit$votes[, "sim"] } else {
      pc <- stats::prcomp(Z, rank. = min(20L, ncol(Z))); d <- data.frame(y = y, pc$x); idx <- sample(rep(1:5, length.out = nrow(d)))
      pr <- numeric(nrow(d)); for (k in 1:5) { f <- suppressWarnings(stats::glm(y ~ ., data = d[idx != k, ], family = stats::binomial()))
        pr[idx == k] <- stats::predict(f, d[idx == k, ], type = "response") } }
    auroc(pr, y == "sim") }, error = function(e) NA_real_)
  ks <- stats::median(sapply(seq_len(nrow(real)), function(j) suppressWarnings(stats::ks.test(Fr[, j], Fs[, j])$statistic)))
  top <- order(rowMeans(real > 0), decreasing = TRUE)[seq_len(min(100L, nrow(real)))]
  Cr <- stats::cor(Fr[, top], method = "spearman"); Cs <- stats::cor(Fs[, top], method = "spearman"); Cr[is.na(Cr)] <- 0; Cs[is.na(Cs)] <- 0
  frob <- sqrt(sum((Cr - Cs)^2)) / length(top)
  data.frame(simulator = simulator, template = template$id, n_real = nrow(Fr), n_sim = nrow(Fs), auroc = auc, ks_median = ks, frob = frob,
             flag_low_realism = is.finite(auc) && auc > 0.9, note = "", stringsAsFactors = FALSE)
}

auroc <- function(score, truth) { r <- rank(score); n1 <- sum(truth); n0 <- sum(!truth); (sum(r[truth]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
