# -----------------------------------------------------------------------------
# Permutation design and index generation
# -----------------------------------------------------------------------------

# Shared exchangeability design for both components.
# Called from run_component_permutations() once per pipeline run; the result
# is reused by both the prevalence and abundance null kernels so the two
# components see the same permutation maps.
build_permutation_design <- function(
    metadata,
    tested_vars = NULL,
    tested_term = NULL,
    cluster_id = NULL,
    strata = NULL,
    cluster_exposure_mode = c("sample", "cluster"),
    verbose = FALSE
) {
  metadata <- as_plain_data_frame(metadata)

  # resolve tested_vars and the cluster/strata vectors.
  tested_vars <- check_tested_vars(metadata, tested_vars, tested_term)
  cluster_exposure_mode <- match.arg(cluster_exposure_mode)

  cluster_vector <- extract_cluster_vector(metadata, cluster_id)
  strata_vector  <- extract_strata_vector(metadata, strata)

  unit_level <- cluster_exposure_mode

  # cluster mode preconditions: cluster_id supplied AND tested exposure
  # constant within each cluster
  if (identical(unit_level, "cluster")) {
    if (is.null(cluster_vector)) {
      stop("cluster_exposure_mode = 'cluster' requires cluster_id.")
    }
    if (!detect_cluster_level_test(metadata, tested_vars, cluster_vector)) {
      stop("cluster_exposure_mode = 'cluster' requires each tested variable to be constant within cluster.")
    }
  }

  if (unit_level == "sample") {
    # sample level: rows are the permutation units.
    # cluster_id supplied -> exchange only within cluster
    # no cluster_id -> single all in one block (free shuffle)
    sample_blocks <- if (!is.null(cluster_vector)) {
      factor(cluster_vector)
    } else {
      factor(rep("all", nrow(metadata)))
    }

    # shuffle only inside the stratum and
    # cluster combination. The interaction builds those combined block ids,
    if (!is.null(strata_vector)) {
      sample_blocks <- interaction(strata_vector, sample_blocks, drop = TRUE, lex.order = TRUE)
    }

    out <- list(
      unit_level = "sample",
      tested_vars = tested_vars,
      sample_blocks = factor(sample_blocks),
      cluster_id = cluster_vector,
      strata = strata_vector,
      n = nrow(metadata)
    )
  } else {
    # cluster level: clusters are the permutation units. the kernel uses a target to donor cluster map.
    cluster_factor <- factor(cluster_vector)
    cluster_levels <- levels(cluster_factor)

    # for cluster permutation under strata, each cluster must lie in exactly one
    # stratum.
    cluster_strata <- NULL
    if (!is.null(strata_vector)) {
      strata_by_cluster <- split(strata_vector, cluster_factor)
      cluster_has_one_stratum <- vapply(
        strata_by_cluster,
        function(x) length(unique(x)) <= 1L,
        logical(1)
      )

      if (!all(cluster_has_one_stratum)) {
        stop("For cluster-level permutation, each cluster must belong to one stratum.")
      }

      # collapse stratum per row to stratum per cluster (one label per cluster)
      cluster_strata <- vapply(strata_by_cluster, function(x) as.character(x[1L]), character(1))
      names(cluster_strata) <- names(strata_by_cluster)
    }

    out <- list(
      unit_level = "cluster",
      tested_vars = tested_vars,
      cluster_factor = cluster_factor,
      cluster_levels = cluster_levels,
      cluster_strata = cluster_strata,
      cluster_id = cluster_vector,
      strata = strata_vector,
      n = nrow(metadata)
    )
  }

  class(out) <- c("pursue_permutation_design", "list")    # tag so downstream kernels can sanity check inputs

  if (verbose) {
    message("Permutation design built.")
    message("  Unit level: ", out$unit_level)
    message("  Tested vars: ", paste(out$tested_vars, collapse = ", "))
    if (out$unit_level == "sample") {
      message("  Sample blocks: ", nlevels(out$sample_blocks))
    } else {
      message("  Clusters: ", length(out$cluster_levels))
    }
  }

  out
}

# Draw n_perm permutation objects according to the design.
# Returned list is consumed by both null kernels in lockstep.
generate_permutation_indices <- function(
    perm_design,
    n_perm = 199,
    seed = 1,
    verbose = FALSE
) {

  if (!is.null(seed)) set.seed(seed)

  permutations <- vector("list", n_perm)

  if (perm_design$unit_level == "sample") {
    # sample level: ask permute for n_perm row shuffles inside sample_blocks.
    # Each row of shuf_mat is a permutation of 1..n that respects the blocking.
    shuf_mat <- permute::shuffleSet(
      n = perm_design$n,
      nset = n_perm,
      control = permute::how(blocks = perm_design$sample_blocks),
      check = FALSE                                                # tolerate duplicate draws
    )

    # wrap each row into the { type, index } shape so kernels can dispatch on $type
    for (i in seq_len(n_perm)) {
      permutations[[i]] <- list(
        type = "sample",
        index = as.integer(shuf_mat[i, ])
      )
    }
  } else {
    # cluster level: shuffle cluster indices 1..K, then map indices to labels.
    # K typically much smaller than n, so the perm matrix is cheaper to draw.
    cluster_levels <- perm_design$cluster_levels
    n_clusters <- length(cluster_levels)

    # strata factor over cluster indices 1..K (single block if no strata supplied)
    strata_blocks <- if (is.null(perm_design$cluster_strata)) {
      factor(rep("all", n_clusters))
    } else {
      factor(perm_design$cluster_strata[cluster_levels])
    }

    shuf_mat <- permute::shuffleSet(
      n = n_clusters,
      nset = n_perm,
      control = permute::how(blocks = strata_blocks),
      check = FALSE
    )

    # each draw becomes a target to donor label map.
    # target = original cluster (names(cluster_map)), donor = shuffled cluster (values).
    for (i in seq_len(n_perm)) {
      permuted_clusters <- cluster_levels[as.integer(shuf_mat[i, ])]
      names(permuted_clusters) <- cluster_levels
      permutations[[i]] <- list(
        type = "cluster",
        cluster_map = permuted_clusters
      )
    }
  }

  if (verbose) {
    message("Generated ", length(permutations), " permutation objects.")
  }

  permutations
}


# -----------------------------------------------------------------------------
# Fixed effect and statistic helpers
# -----------------------------------------------------------------------------

# Wald stat from precomputed QR and (Z^T Z)^-1. This is the only function that runs
# in the inner null loop.
linear_block_wald_from_precomp <- function(
    y,
    qr_full,
    xtx_inverse,
    n_nuisance_coef,
    n_test_coef
) {
  coef_estimates <- qr.coef(qr_full, y)                    # full coefficient vector
  residuals <- qr.resid(qr_full, y)                        # r_hat = y - Z coef

  n_obs <- length(y)                                       # n
  n_coef <- n_nuisance_coef + n_test_coef                  # p
  df_resid <- n_obs - n_coef                               # n - p

  # guard: not enough samples to estimate residual variance
  if (!is.finite(df_resid) || df_resid <= 0L) {
    return(list(
      stat = 0,
      df = NA_real_, status = "df_resid_nonpositive",
      estimate = NA_real_, se = NA_real_
    ))
  }

  rss <- sum(residuals^2)                                  # sum_i r_hat_i^2
  sigma2 <- rss / df_resid                                 # sigma_hat^2 = sum_i r_hat_i^2 / (n - p)

  # guard: sigma_hat^2 NaN or negative under degenerate numerics
  if (!is.finite(sigma2) || sigma2 < 0) {
    return(list(
      stat = 0,
      df = NA_real_, status = "sigma2_invalid",
      estimate = NA_real_, se = NA_real_
    ))
  }

  # prepare_component_fl_cache() builds full_model_matrix = cbind(nuisance, tested),
  # so the tested indices are always (p_nuis + 1) through (p_nuis + p_test).
  # The kernels depend on this layout when they index into xtx_inverse.
  tested_coef_index <- seq.int(n_nuisance_coef + 1L, n_nuisance_coef + n_test_coef)

  coef_tested <- coef_estimates[tested_coef_index]                                         # theta_hat
  vcov_tested <- sigma2 * xtx_inverse[tested_coef_index, tested_coef_index, drop = FALSE]  # sigma_hat^2 [(Z^T Z)^-1]_{theta theta} (scalar) or Sigma_hat (multi df)

  # delegate the scalar vs multi df formula choice to wald_stat()
  wald_result <- wald_stat(coef_estimates = coef_tested, vcov_matrix = vcov_tested)

  # only the 1 df case has a meaningful scalar estimate and SE; for multi df
  # tested terms the estimate is a vector and the SE is per coefficient, so
  # we return NA at this level and let downstream code rely on the stat alone
  estimate_out <- if (n_test_coef == 1L) unname(coef_tested[1L]) else NA_real_
  se_out       <- if (n_test_coef == 1L) wald_result$se          else NA_real_

  list(
    stat = unname(wald_result$stat),
    df = if (n_test_coef == 1L) df_resid else wald_result$df,
    status = wald_result$status,
    estimate = estimate_out,
    se = se_out
  )
}

# Median centered residual winsorization
# Centers r_red on med(r), clips at q_0.05 / q_0.95 of the centered values,
# then adds med(r) back to produce r_red_wins. Returns input unchanged when
# finite residuals are fewer than min_n.
winsorize_residuals <- function(e, quantile = 0.05, min_n = 20L) {
  residuals <- as.numeric(e)
  finite_rows <- is.finite(residuals)
  residuals_out <- residuals

  # too few finite residuals
  if (sum(finite_rows) < as.integer(min_n)) {
    return(list(
      e_wins = residuals_out, lower = NA_real_, upper = NA_real_,
      n_wins_low = 0L, n_wins_high = 0L,
      applied = FALSE, status = "too_few_residuals"
    ))
  }

  finite_residuals <- residuals[finite_rows]
  residual_center <- stats::median(finite_residuals)            # med(r)
  centered_residuals <- finite_residuals - residual_center      # r_red - med(r)

  lower_cut <- as.numeric(stats::quantile(centered_residuals, probs = quantile,     names = FALSE, na.rm = TRUE))
  upper_cut <- as.numeric(stats::quantile(centered_residuals, probs = 1 - quantile, names = FALSE, na.rm = TRUE))

  clipped <- pmin(pmax(centered_residuals, lower_cut), upper_cut)   # clip
  winsorized_residuals <- residual_center + clipped                  # r_red_wins = med(r) + clipped

  n_low  <- sum(centered_residuals < lower_cut, na.rm = TRUE)
  n_high <- sum(centered_residuals > upper_cut, na.rm = TRUE)

  residuals_out[finite_rows] <- winsorized_residuals

  list(
    e_wins = residuals_out,
    lower = residual_center + lower_cut,
    upper = residual_center + upper_cut,
    n_wins_low = as.integer(n_low),
    n_wins_high = as.integer(n_high),
    applied = TRUE,
    status = "ok"
  )
}


# -----------------------------------------------------------------------------
# Freedman Lane component cache construction
# -----------------------------------------------------------------------------

# Build the per-taxon reusable cache (reduced fit, residuals,
# full QR, (Z^T Z)^-1, observed Wald).
prepare_component_fl_cache <- function(
    response_matrix,
    metadata,
    formula,
    tested_term,
    tested_vars = NULL,
    fit_args = list(),
    neutral_stat = 0,
    verbose = FALSE,
    component = c("prevalence", "abundance"),
    cache_label = "Abundance"
) {
  component <- match.arg(component)
  response_matrix <- as_numeric_matrix(response_matrix, "response_matrix")

  # fit_args knobs are component specific. Prevalence uses min_nonzero_n = 1
  # and no residual winsorization;
  min_nonzero_n         <- fit_args$min_nonzero_n
  resid_winsorize       <- isTRUE(fit_args$resid_winsorize)
  resid_winsor_quantile <- fit_args$resid_winsor_quantile
  resid_winsor_min_n    <- fit_args$resid_winsor_min_n

  if (is.null(min_nonzero_n)) min_nonzero_n <- 1L
  if (is.null(resid_winsor_quantile)) resid_winsor_quantile <- 0.05
  if (is.null(resid_winsor_min_n)) resid_winsor_min_n <- 20L

  formula_full    <- make_response_formula(formula, response_name = ".response")
  formula_info    <- match_tested_term(formula_full, tested_term = tested_term)
  tested_term_idx <- formula_info$tested_term_idx

  n_taxa <- ncol(response_matrix)
  taxa_names <- colnames(response_matrix)
  if (is.null(taxa_names)) taxa_names <- paste0("Taxon_", seq_len(n_taxa))

  # parallel summary vectors mirrored from per-taxon entries
  taxon_entries <- vector("list", n_taxa)
  observed_stat        <- rep(NA_real_, n_taxa)
  observed_estimate    <- rep(NA_real_, n_taxa)
  observed_se          <- rep(NA_real_, n_taxa)
  observed_df          <- rep(NA_real_, n_taxa)
  stat_for_permutation <- rep(neutral_stat, n_taxa)
  n_obs_vec            <- rep(NA_integer_, n_taxa)
  usable               <- rep(FALSE, n_taxa)
  status               <- rep(NA_character_, n_taxa)
  reason_unusable      <- rep(NA_character_, n_taxa)

  for (taxon_index in seq_len(n_taxa)) {
    taxon_response <- response_matrix[, taxon_index]

    # taxon-specific model frame: response column attached, rows missing any
    # formula variable dropped, factor levels relevelled
    fit_data <- metadata
    fit_data$.response <- as.numeric(taxon_response)

    needed_vars <- unique(c(".response", all.vars(formula_full)))
    complete_rows <- stats::complete.cases(fit_data[, needed_vars, drop = FALSE])

    taxon_fit_data <- fit_data[complete_rows, , drop = FALSE]
    taxon_fit_data <- droplevels(taxon_fit_data)

    row_idx <- which(complete_rows)               # full sample row indices this taxon uses
    observed_response <- taxon_fit_data$.response
    n_obs <- nrow(taxon_fit_data)

    # entry template with neutralized defaults; success path overwrites the FL fields
    entry <- list(
      taxon_index = taxon_index,
      taxon_name = taxa_names[taxon_index],
      row_idx = row_idx,
      n_obs = n_obs,
      qr_full = NULL, XtX_inv = NULL,
      yhat_red = NULL, e_red = NULL,
      p_nuis = NA_integer_, p_test = NA_integer_,
      observed_stat = NA_real_,
      observed_estimate = NA_real_,
      observed_se = NA_real_,
      observed_df = NA_real_,
      stat_for_permutation = neutral_stat,
      usable_for_testing = FALSE,
      status = NA_character_,
      reason_unusable = NA_character_,
      winsor_applied = FALSE, winsor_status = "not_used",
      winsor_lower = NA_real_, winsor_upper = NA_real_,
      n_wins_low = 0L, n_wins_high = 0L
    )

    # validity gates
    fail_reason <- NULL

    if (n_obs < min_nonzero_n) {
      fail_reason <- "too_few_nonzero_samples"
    } else if (length(unique(observed_response)) < 2L) {
      fail_reason <- "constant_response"
    } else if (!tested_term_has_variation(taxon_fit_data, tested_term)) {
      fail_reason <- "insufficient_tested_term_variation"
    }

    if (is.null(fail_reason)) {
      # build the model matrix; wrap in tryCatch since contrasts on collapsed factors can throw errors
      model_matrix <- tryCatch(
        stats::model.matrix(formula_info$terms_obj, data = taxon_fit_data),
        error = function(e) NULL
      )
      if (is.null(model_matrix)) {
        fail_reason <- "model_matrix_failed"
      }
    }

    if (is.null(fail_reason)) {
      design_split <- split_tested_and_nuisance_columns(model_matrix, tested_term_idx = tested_term_idx)
      if (length(design_split$tested_cols) < 1L) {
        fail_reason <- "tested_block_missing"
      }
    }

    if (is.null(fail_reason)) {
      nuisance_matrix <- as.matrix(model_matrix[, design_split$nuisance_cols, drop = FALSE])
      tested_matrix   <- as.matrix(model_matrix[, design_split$tested_cols,   drop = FALSE])

      # tested columns can be present but  flat if every
      # sample landed in the same level after droplevels. 1e-12 floor keeps variation.
      tested_sd <- apply(tested_matrix, 2L, function(x) stats::sd(x, na.rm = TRUE))
      if (all(!is.finite(tested_sd) | tested_sd < 1e-12)) {
        fail_reason <- "degenerate_test_block"
      }
    }

    if (!is.null(fail_reason)) {
      entry$status <- fail_reason
      entry$reason_unusable <- fail_reason
      taxon_entries[[taxon_index]] <- entry
      n_obs_vec[taxon_index] <- n_obs
      status[taxon_index] <- fail_reason
      reason_unusable[taxon_index] <- fail_reason
      next
    }

    # reduced model: FL permutes the residuals of THIS fit
    nuisance_qr <- qr(nuisance_matrix)
    yhat_red <- qr.fitted(nuisance_qr, observed_response)        # y_hat_red
    e_red_raw <- qr.resid(nuisance_qr, observed_response)        # r_red (raw)

    # Optional winsorization
    winsor_object <- if (isTRUE(resid_winsorize)) {
      winsorize_residuals(
        e = e_red_raw,
        quantile = resid_winsor_quantile,
        min_n = resid_winsor_min_n
      )
    } else {
      list(e_wins = e_red_raw, lower = NA_real_, upper = NA_real_,
           n_wins_low = 0L, n_wins_high = 0L,
           applied = FALSE, status = "not_used")
    }

    e_red <- winsor_object$e_wins
    y_obs_used <- yhat_red + e_red                               # observed response actually used for Wald

    # full design cache
    # Tested block goes LAST so linear_block_wald_from_precomp() can index it
    # at (p_nuis + 1):(p_nuis + p_test)
    full_model_matrix <- cbind(nuisance_matrix, tested_matrix)
    qr_full <- qr(full_model_matrix)

    # rank deficiency means collinear columns; Wald cannot be computed
    if (qr_full$rank < ncol(full_model_matrix)) {
      entry$status <- "design_rank_deficient"
      entry$reason_unusable <- "design_rank_deficient"
      taxon_entries[[taxon_index]] <- entry
      n_obs_vec[taxon_index] <- n_obs
      status[taxon_index] <- "design_rank_deficient"
      reason_unusable[taxon_index] <- "design_rank_deficient"
      next
    }

    # chol2inv(qr.R(...)) gives (Z^T Z)^-1 (cheaper + more stable than solve(crossprod(Z)))
    xtx_inverse <- tryCatch(chol2inv(qr.R(qr_full)), error = function(e) NULL)
    if (is.null(xtx_inverse)) {
      entry$status <- "xtx_inverse_failed"
      entry$reason_unusable <- "xtx_inverse_failed"
      taxon_entries[[taxon_index]] <- entry
      n_obs_vec[taxon_index] <- n_obs
      status[taxon_index] <- "xtx_inverse_failed"
      reason_unusable[taxon_index] <- "xtx_inverse_failed"
      next
    }

    # observed Wald, using the same cached pieces the null kernel reuses per perm
    observed_test <- linear_block_wald_from_precomp(
      y = y_obs_used,
      qr_full = qr_full,
      xtx_inverse = xtx_inverse,
      n_nuisance_coef = ncol(nuisance_matrix),
      n_test_coef = ncol(tested_matrix)
    )

    # fill FL cache fields on the entry
    entry$qr_full <- qr_full
    entry$XtX_inv <- xtx_inverse
    entry$yhat_red <- yhat_red
    entry$e_red <- e_red                         # post winsorization (what the kernel actually shuffles)
    entry$winsor_applied <- isTRUE(winsor_object$applied)
    entry$winsor_status <- winsor_object$status
    entry$winsor_lower <- winsor_object$lower
    entry$winsor_upper <- winsor_object$upper
    entry$n_wins_low <- winsor_object$n_wins_low
    entry$n_wins_high <- winsor_object$n_wins_high
    entry$p_nuis <- ncol(nuisance_matrix)
    entry$p_test <- ncol(tested_matrix)
    entry$usable_for_testing <- identical(observed_test$status, "ok")
    entry$status             <- observed_test$status
    entry$reason_unusable    <- if (entry$usable_for_testing) NA_character_ else observed_test$status

    if (entry$usable_for_testing) {
      entry$observed_stat <- observed_test$stat
      entry$observed_estimate <- observed_test$estimate
      entry$observed_se <- observed_test$se
      entry$observed_df <- observed_test$df
      entry$stat_for_permutation <- abs(observed_test$stat)
    }

    # mirror onto parallel summary vectors
    taxon_entries[[taxon_index]] <- entry
    observed_stat[taxon_index] <- entry$observed_stat
    observed_estimate[taxon_index] <- entry$observed_estimate
    observed_se[taxon_index] <- entry$observed_se
    observed_df[taxon_index] <- entry$observed_df
    stat_for_permutation[taxon_index] <- entry$stat_for_permutation
    n_obs_vec[taxon_index] <- n_obs
    usable[taxon_index] <- isTRUE(entry$usable_for_testing)
    status[taxon_index] <- entry$status
    reason_unusable[taxon_index] <- entry$reason_unusable
  }

  # per-taxon entries + parallel summary vectors + fit_args knobs.
  out <- list(
    taxon_entries = taxon_entries,
    taxon_name = taxa_names,
    observed_stat = observed_stat,
    observed_estimate = observed_estimate,
    observed_se = observed_se,
    observed_df = observed_df,
    stat_for_permutation = stat_for_permutation,
    n_obs = n_obs_vec,
    usable_for_testing = usable,
    status = status,
    reason_unusable = reason_unusable,
    n_taxa = n_taxa,
    neutral_stat = neutral_stat,
    component = component,
    tested_term = tested_term,
    tested_vars = tested_vars,
    statistic_type = "fl_block_wald",
    resid_winsorize = resid_winsorize,
    resid_winsor_quantile = resid_winsor_quantile,
    resid_winsor_min_n = resid_winsor_min_n
  )

  if (verbose) {
    message("Prepared ", cache_label, " FL cache.")
    message("  Taxa: ", n_taxa)
    message("  Usable for testing: ", sum(usable))
    message("  Too few nonzero: ", sum(status == "too_few_nonzero_samples", na.rm = TRUE))
    message("  Degenerate tested block: ", sum(status == "degenerate_test_block", na.rm = TRUE))
    if (isTRUE(resid_winsorize)) {
      message("  Taxa with winsorization applied: ",
              sum(vapply(taxon_entries, function(x) isTRUE(x$winsor_applied), logical(1))))
    }
    message(strrep("-", 50))
  }

  out
}


# -----------------------------------------------------------------------------
# Permutation object and row context helpers
# -----------------------------------------------------------------------------

# Recover the full sample cluster_id vector from a permutation design.
# Returns NULL when the design has no clusters (sample level with no blocking).
extract_full_cluster_id <- function(perm_design) {
  if (is.null(perm_design$cluster_id)) return(NULL)
  as.character(perm_design$cluster_id)
}

# Build the per taxon row context: full-sample row indices + cluster label
# per local row. Only used by the cluster mover; sample paths ignore it.
build_taxon_row_context <- function(row_idx, perm_design) {
  cluster_id <- extract_full_cluster_id(perm_design)
  list(
    row_idx = row_idx,
    cluster_id_row = if (is.null(cluster_id)) NULL else as.character(cluster_id[row_idx])
  )
}

# Apply a cluster permutation to per taxon residuals: for every target
# cluster, copy the residual block from its donor cluster into the target
# slots. Returns ok = FALSE with a reason string when a donor is missing
# or its local valid row count disagrees with the target's, so the kernel
# can neutralize that (taxon, perm) cell without aborting the whole run.
apply_cluster_map_to_residual_blocks <- function(e_red_full, row_context, perm_obj) {
  cluster_map    <- perm_obj$cluster_map
  row_idx        <- row_context$row_idx
  cluster_id_row <- row_context$cluster_id_row

  permuted_residuals <- rep(NA_real_, length(row_idx))
  target_rows_by_cluster <- split(seq_along(row_idx), cluster_id_row)

  for (target_cluster in names(target_rows_by_cluster)) {
    donor_cluster <- unname(cluster_map[target_cluster])

    if (length(donor_cluster) != 1L || is.na(donor_cluster) || !nzchar(donor_cluster)) {
      return(list(ok = FALSE, reason = paste0("missing_donor:", target_cluster), values = NULL))
    }

    target_pos <- target_rows_by_cluster[[target_cluster]]
    donor_pos  <- which(cluster_id_row == donor_cluster)

    # Different taxa have different complete-case rows; local cluster
    # sizes can disagree -> neutralize this (taxon, perm) cell.
    if (length(target_pos) != length(donor_pos)) {
      return(list(
        ok = FALSE,
        reason = paste0(
          "local_cluster_size_mismatch:", target_cluster, "->", donor_cluster,
          " (", length(target_pos), " vs ", length(donor_pos), ")"
        ),
        values = NULL
      ))
    }

    donor_rows_full <- row_idx[donor_pos]
    permuted_residuals[target_pos] <- e_red_full[donor_rows_full]
  }

  list(ok = TRUE, reason = "ok", values = permuted_residuals)
}


# -----------------------------------------------------------------------------
# Full length permutation cache layout
# -----------------------------------------------------------------------------

# Wrap the base FL cache with full-sample residual/fit vectors and row
# contexts so the null kernels can permute by unit even when taxa have
# different complete-case rows.
prepare_component_permutation_cache <- function(
    response_matrix,
    metadata,
    formula,
    tested_term,
    tested_vars = NULL,
    fit_args = list(),
    neutral_stat = 0,
    perm_design,
    verbose = FALSE,
    component = c("prevalence", "abundance"),
    cache_label = "Abundance"
) {
  if (!inherits(perm_design, "pursue_permutation_design")) {
    stop("perm_design must come from build_permutation_design().")
  }

  component <- match.arg(component)

  # build the base FL cache first (per taxon QR, residuals, observed Wald)
  base_cache <- prepare_component_fl_cache(
    response_matrix = response_matrix,
    metadata = metadata,
    formula = formula,
    tested_term = tested_term,
    tested_vars = tested_vars,
    fit_args = fit_args,
    neutral_stat = neutral_stat,
    verbose = verbose,
    component = component,
    cache_label = cache_label
  )

  n_full <- nrow(response_matrix)
  if (n_full != perm_design$n) {
    stop("nrow(response_matrix) must equal perm_design$n.")
  }

  # lift e_red into full-sample coordinates so the cluster mover can index by
  # full-table row; attach the cluster row context for the same reason
  for (taxon_index in seq_len(base_cache$n_taxa)) {
    entry <- base_cache$taxon_entries[[taxon_index]]

    # NA padded; neutralized taxa keep NA everywhere (kernel skips them via usable_for_testing)
    e_red_full <- rep(NA_real_, n_full)
    if (!is.null(entry$e_red)) e_red_full[entry$row_idx] <- entry$e_red

    entry$e_red_full  <- e_red_full
    entry$row_context <- build_taxon_row_context(row_idx = entry$row_idx, perm_design = perm_design)

    base_cache$taxon_entries[[taxon_index]] <- entry
  }

  base_cache$perm_design_unit_level <- perm_design$unit_level
  base_cache$perm_design_n <- perm_design$n
  base_cache$cache_layout <- "full_length"

  class(base_cache) <- c(
    paste0("pursue_", component, "_component_cache"),
    "pursue_component_permutation_cache",
    "list"
  )
  base_cache
}

# Restrict a global sample permutation to this taxon's row subset, returning
# a local permutation consistent with the global ordering of those rows.
induce_local_sample_permutation <- function(row_idx, perm_obj) {
  order(match(row_idx, perm_obj$index))
}


# -----------------------------------------------------------------------------
# Freedman Lane null kernels (sample and cluster level)
# -----------------------------------------------------------------------------

# Compute one sample level null column: stat per taxon for a single perm.
# Called once per perm_pos by compute_component_fl_null_run() in either a
# serial lapply or a parLapplyLB. The output list { b, stat, status } is
# stitched back into the (taxa x perm) matrix by the caller.
compute_component_fl_null_sample_one_perm <- function(
    perm_pos,
    component_cache,
    perm_list,
    perm_design,
    neutral_stat = 0
) {
  permutation <- perm_list[[perm_pos]]
  if (!identical(permutation$type, "sample")) {
    stop("Sample null kernel received non-sample permutation object.")
  }

  n_taxa <- component_cache$n_taxa
  stat_col   <- rep(neutral_stat, n_taxa)
  status_col <- rep("neutralized", n_taxa)

  for (taxon_index in seq_len(n_taxa)) {
    entry <- component_cache$taxon_entries[[taxon_index]]

    if (!isTRUE(entry$usable_for_testing)) {
      status_col[taxon_index] <- "unusable_component"
      next
    }

    # per-cell skip: rare numerical degeneracy in cached pieces
    if (any(!is.finite(entry$yhat_red)) || any(!is.finite(entry$e_red))) {
      status_col[taxon_index] <- "nonfinite_cached_residual_or_fit"
      next
    }

    # FL null pseudoresponse: y* = yhat_red + perm(r_red)
    local_index <- induce_local_sample_permutation(row_idx = entry$row_idx, perm_obj = permutation)
    permuted_residuals <- entry$e_red[local_index]
    y_star <- entry$yhat_red + permuted_residuals

    if (any(!is.finite(y_star))) {
      status_col[taxon_index] <- "nonfinite_y_star"
      next
    }

    permuted_test <- linear_block_wald_from_precomp(
      y = y_star,
      qr_full = entry$qr_full,
      xtx_inverse = entry$XtX_inv,
      n_nuisance_coef = entry$p_nuis,
      n_test_coef = entry$p_test
    )

    if (identical(permuted_test$status, "ok")) {
      stat_col[taxon_index]   <- abs(permuted_test$stat)    # |W| for two-sided null
      status_col[taxon_index] <- "ok"
    } else {
      status_col[taxon_index] <- permuted_test$status
    }
  }

  list(b = perm_pos, stat = stat_col, status = status_col)
}

# Shared null orchestration. Loops a per-perm kernel across all permutations
# (serial or parallel) and assembles the (taxa x perm) matrices.
# The kernel + the extra worker export names vary between sample and cluster
# paths; everything else is shared.
compute_component_fl_null_run <- function(
    component_cache,
    perm_list,
    perm_design,
    kernel_fn,
    kernel_name,
    extra_exports = character(0),
    neutral_stat = 0,
    parallelize_permutations = FALSE,
    permutation_n_cores = max(1L, parallel::detectCores() - 1L),
    verbose = FALSE,
    null_label = "Abundance"
) {
  if (!inherits(component_cache, "pursue_component_permutation_cache")) {
    stop("component_cache must come from prepare_component_permutation_cache().")
  }

  unit_level <- perm_design$unit_level
  n_taxa     <- component_cache$n_taxa
  n_perm     <- length(perm_list)

  # preallocate, dimnames so callers can index by taxon name/perm id
  stat_mat <- matrix(
    neutral_stat,
    nrow = n_taxa, ncol = n_perm,
    dimnames = list(component_cache$taxon_name, paste0("perm_", seq_len(n_perm)))
  )
  status_mat <- matrix("neutralized", nrow = n_taxa, ncol = n_perm, dimnames = dimnames(stat_mat))

  # do not parallelize when the permutation count is too small to amortize
  # worker startup overhead.
  use_parallel <- isTRUE(parallelize_permutations) &&
    n_perm > permutation_n_cores

  if (use_parallel) {
    n_workers <- min(as.integer(permutation_n_cores), n_perm)

    if (verbose) {
      message(null_label, " FL null (", unit_level, "): parallelizing ", n_perm,
              " permutations across ", n_workers, " workers.")
    }

    cl <- parallel::makeCluster(n_workers, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # export everything the kernel needs on the workers
    parallel::clusterExport(
      cl,
      varlist = c(
        "component_cache", "perm_list", "perm_design", "neutral_stat",
        kernel_name, extra_exports,
        "linear_block_wald_from_precomp", "wald_stat"
      ),
      envir = environment()
    )

    # workers pull next perm when free
    results <- parallel::parLapplyLB(cl, X = seq_len(n_perm), fun = function(perm_pos) {
      kernel_fn(
        perm_pos, component_cache = component_cache, perm_list = perm_list,
        perm_design = perm_design, neutral_stat = neutral_stat
      )
    })
  } else {
    # serial progress
    progress_step <- max(1L, floor(n_perm / 10L))
    results <- vector("list", n_perm)
    for (perm_pos in seq_len(n_perm)) {
      results[[perm_pos]] <- kernel_fn(
        perm_pos, component_cache = component_cache, perm_list = perm_list,
        perm_design = perm_design, neutral_stat = neutral_stat
      )

      if (verbose && (perm_pos %% progress_step == 0L || perm_pos == n_perm)) {
        message(null_label, " FL null (", unit_level, "): ", perm_pos, "/", n_perm)
      }
    }
  }

  # result$b carries the original perm_pos
  for (result in results) {
    stat_mat[, result$b]   <- result$stat
    status_mat[, result$b] <- result$status
  }

  out <- list(
    stat_mat = stat_mat,
    status_mat = status_mat,
    diagnostics = list(
      unit_level = unit_level,
      n_taxa = n_taxa, n_perm = n_perm,
      n_ok = sum(status_mat == "ok", na.rm = TRUE),
      n_neutralized = sum(status_mat != "ok", na.rm = TRUE),
      neutralization_reasons = sort(table(status_mat[status_mat != "ok"]), decreasing = TRUE)
    )
  )
  out
}


# Compute one cluster level null column: stat per taxon for a single perm.
# Same shape as the sample kernel but with one extra step: residuals are
# transplanted in cluster blocks (donor cluster to target slots) before
# the Wald refit. Called by compute_component_fl_null_run() per perm.
compute_component_fl_null_cluster_one_perm <- function(
    perm_pos,
    component_cache,
    perm_list,
    perm_design,
    neutral_stat = 0
) {
  permutation <- perm_list[[perm_pos]]

  n_taxa <- component_cache$n_taxa
  stat_col   <- rep(neutral_stat, n_taxa)
  status_col <- rep("neutralized", n_taxa)

  for (taxon_index in seq_len(n_taxa)) {
    entry <- component_cache$taxon_entries[[taxon_index]]

    if (!isTRUE(entry$usable_for_testing)) {
      status_col[taxon_index] <- "unusable_component"
      next
    }

    permuted_residual_object <- apply_cluster_map_to_residual_blocks(
      e_red_full = entry$e_red_full,
      row_context = entry$row_context,
      perm_obj = permutation
    )

    if (!isTRUE(permuted_residual_object$ok)) {
      # propagate the mover failure reason (e.g. local_cluster_size_mismatch)
      status_col[taxon_index] <- permuted_residual_object$reason
      next
    }

    # FL pseudoresponse (cluster blocks transplanted in place of single-row shuffles)
    permuted_residuals <- permuted_residual_object$values
    y_star <- entry$yhat_red + permuted_residuals

    if (any(!is.finite(y_star))) {
      status_col[taxon_index] <- "nonfinite_y_star"
      next
    }

    permuted_test <- linear_block_wald_from_precomp(
      y = y_star,
      qr_full = entry$qr_full,
      xtx_inverse = entry$XtX_inv,
      n_nuisance_coef = entry$p_nuis,
      n_test_coef = entry$p_test
    )

    if (identical(permuted_test$status, "ok")) {
      stat_col[taxon_index]   <- abs(permuted_test$stat)         # |W| for two sided null
      status_col[taxon_index] <- "ok"
    } else {
      status_col[taxon_index] <- permuted_test$status
    }
  }

  list(b = perm_pos, stat = stat_col, status = status_col)
}

# picks the per-perm kernel + extra worker exports by unit_level,
# then delegates the orchestration scaffolding to compute_component_fl_null_run().
compute_component_permutation_null <- function(
    component_cache,
    perm_list,
    perm_design,
    neutral_stat = 0,
    parallelize_permutations = FALSE,
    permutation_n_cores = max(1L, parallel::detectCores() - 1L),
    verbose = FALSE,
    null_label = "Abundance"
) {
  unit_level <- perm_design$unit_level

  if (identical(unit_level, "sample")) {
    kernel_fn     <- compute_component_fl_null_sample_one_perm
    kernel_name   <- "compute_component_fl_null_sample_one_perm"
    extra_exports <- "induce_local_sample_permutation"
  } else if (identical(unit_level, "cluster")) {
    kernel_fn     <- compute_component_fl_null_cluster_one_perm
    kernel_name   <- "compute_component_fl_null_cluster_one_perm"
    extra_exports <- "apply_cluster_map_to_residual_blocks"
  } else {
    stop("Unsupported perm_design$unit_level.")
  }

  compute_component_fl_null_run(
    component_cache = component_cache,
    perm_list = perm_list,
    perm_design = perm_design,
    kernel_fn = kernel_fn,
    kernel_name = kernel_name,
    extra_exports = extra_exports,
    neutral_stat = neutral_stat,
    parallelize_permutations = parallelize_permutations,
    permutation_n_cores = permutation_n_cores,
    verbose = verbose,
    null_label = null_label
  )
}
