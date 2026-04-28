# -----------------------------------------------------------------------------
# Permutation design and index generation
# -----------------------------------------------------------------------------

# Define the exchangeability structure used by downstream null generation
build_permutation_design <- function(
    metadata,
    tested_vars = NULL,
    tested_term = NULL,
    cluster_id = NULL,
    strata = NULL,
    cluster_exposure_mode = c("auto", "sample", "cluster"),
    verbose = FALSE
) {
  metadata <- as_plain_data_frame(metadata)

  # Validate the tested variables and normalize grouping vectors before deciding the unit level.
  tested_vars <- check_tested_vars(metadata, tested_vars, tested_term)
  cluster_exposure_mode <- match.arg(cluster_exposure_mode)

  cluster_vector <- extract_cluster_vector(metadata, cluster_id)
  strata_vector <- extract_strata_vector(metadata, strata)

  # Detect whether the tested exposure is cluster-level rather than sample-level.
  cluster_test_constant <- if (!is.null(cluster_vector)) {
    detect_cluster_level_test(metadata, tested_vars, cluster_vector)
  } else {
    FALSE
  }

  # Resolve automatic unit selection from the within-cluster exposure structure.
  if (cluster_exposure_mode == "auto") {
    unit_level <- if (cluster_test_constant) {
      "cluster"
    } else {
      "sample"
    }
  } else {
    unit_level <- cluster_exposure_mode
  }

  # Block invalid cluster-level permutations before constructing a misleading design.
  if (identical(unit_level, "cluster") && !cluster_test_constant) {
    stop(
      "cluster_exposure_mode = 'cluster' requires each tested variable to be constant within cluster."
    )
  }

  # For sample-level tests, samples are permuted within cluster and optional stratum blocks.
  if (unit_level == "sample") {
    sample_blocks <- if (!is.null(cluster_vector)) {
      factor(cluster_vector)
    } else {
      factor(rep("all", nrow(metadata)))
    }

    if (!is.null(strata_vector)) {
      sample_blocks <- interaction(
        strata_vector,
        sample_blocks,
        drop = TRUE,
        lex.order = TRUE
      )
    }

    out <- list(
      unit_level = "sample",
      tested_vars = tested_vars,
      sample_blocks = factor(sample_blocks),
      cluster_id = cluster_vector,
      strata = strata_vector,
      n = nrow(metadata)
    )

    # For cluster-level tests, whole clusters are the permuted exchangeability units.
  } else {
    if (is.null(cluster_vector)) {
      stop("cluster_exposure_mode = 'cluster' requires cluster_id.")
    }

    cluster_factor <- factor(cluster_vector)
    cluster_levels <- levels(cluster_factor)

    # Collapse sample-level strata to cluster-level labels when stratified permutation is requested.
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

      cluster_strata <- vapply(
        strata_by_cluster,
        function(x) as.character(x[1L]),
        character(1)
      )
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

  # Tag the object so null-generation functions can reject incompatible inputs early.
  class(out) <- c("pursue_permutation_design", "list")

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

# Generate sample- or cluster-level permutation objects within exchangeability blocks
generate_permutation_indices <- function(
    perm_design,
    n_perm = 199,
    seed = NULL,
    include_identity = FALSE,
    verbose = FALSE
) {
  if (!inherits(perm_design, "pursue_permutation_design")) {
    stop("perm_design must come from build_permutation_design().")
  }

  # Seed only the permutation generation step, leaving downstream fitting deterministic.
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Allocate all permutation objects and optionally reserve the first slot for identity.
  permutations <- vector("list", n_perm + as.integer(include_identity))
  first_perm_index <- 1L

  # Include the identity map when observed-like diagnostics are needed.
  if (include_identity) {
    permutations[[1L]] <- if (perm_design$unit_level == "sample") {
      list(type = "sample", index = seq_len(perm_design$n))
    } else {
      identity_map <- setNames(perm_design$cluster_levels, perm_design$cluster_levels)
      list(type = "cluster", cluster_map = identity_map)
    }

    first_perm_index <- 2L
  }

  # Generate each permutation according to the unit level encoded in the design.
  for (perm_pos in first_perm_index:length(permutations)) {
    if (perm_design$unit_level == "sample") {
      # Shuffle samples independently inside each exchangeability block.
      block_factor <- perm_design$sample_blocks
      sample_index <- seq_len(perm_design$n)

      permuted_index <- sample_index
      block_indices <- split(sample_index, block_factor)

      for (block_name in names(block_indices)) {
        block_rows <- block_indices[[block_name]]

        if (length(block_rows) > 1L) {
          permuted_index[block_rows] <- sample(
            block_rows,
            length(block_rows),
            replace = FALSE
          )
        } else {
          permuted_index[block_rows] <- block_rows
        }
      }

      permutations[[perm_pos]] <- list(type = "sample", index = permuted_index)
    } else {
      # Shuffle cluster labels globally or within cluster-level strata.
      cluster_levels <- perm_design$cluster_levels

      if (is.null(perm_design$cluster_strata)) {
        permuted_clusters <- sample(cluster_levels, length(cluster_levels), replace = FALSE)
      } else {
        permuted_clusters <- cluster_levels
        clusters_by_stratum <- split(
          cluster_levels,
          perm_design$cluster_strata[cluster_levels]
        )

        for (stratum_name in names(clusters_by_stratum)) {
          stratum_clusters <- clusters_by_stratum[[stratum_name]]

          if (length(stratum_clusters) > 1L) {
            permuted_clusters[match(stratum_clusters, cluster_levels)] <- sample(
              stratum_clusters,
              length(stratum_clusters),
              replace = FALSE
            )
          } else {
            permuted_clusters[match(stratum_clusters, cluster_levels)] <- stratum_clusters
          }
        }
      }

      names(permuted_clusters) <- cluster_levels
      permutations[[perm_pos]] <- list(type = "cluster", cluster_map = permuted_clusters)
    }
  }

  # Store as a typed list so callers do not treat raw indices as a design object.
  class(permutations) <- c("pursue_permutation_list", "list")

  if (verbose) {
    message("Generated ", length(permutations), " permutation objects.")
  }

  permutations
}

# -----------------------------------------------------------------------------
# Fixed-effect and linear-statistic helpers
# -----------------------------------------------------------------------------

# Resolve the tested fixed-effect term inside the model formula
get_fixed_formula_terms <- function(formula_full, tested_term) {
  fixed_formula <- formula_full
  terms_obj <- stats::terms(fixed_formula)
  term_labels <- attr(terms_obj, "term.labels")

  # The tested term must map to one formula term before the design is split into blocks.
  tested_term_idx <- which(term_labels == tested_term)
  if (length(tested_term_idx) != 1L) {
    stop(
      "FL abundance cache requires tested_term to match exactly one fixed-effect term label in the formula. ",
      "Available fixed terms: ", paste(term_labels, collapse = ", ")
    )
  }

  list(
    fixed_formula = fixed_formula,
    terms_obj = terms_obj,
    term_labels = term_labels,
    tested_term_idx = tested_term_idx
  )
}

# Compute a tested-block Wald statistic from precomputed full-design pieces
linear_block_wald_from_precomp <- function(
    y,
    qr_full,
    xtx_inverse,
    n_nuisance_coef,
    n_test_coef
) {
  # Refit the full design against a supplied response without rebuilding the model matrix.
  coef_estimates <- qr.coef(qr_full, y)
  residuals <- qr.resid(qr_full, y)

  n_obs <- length(y)
  n_coef <- n_nuisance_coef + n_test_coef
  df_resid <- n_obs - n_coef

  # Neutralize statistics when residual degrees of freedom cannot support variance estimation.
  if (!is.finite(df_resid) || df_resid <= 0L) {
    return(list(
      stat = 0,
      stat_type = "df_resid_nonpositive",
      df = NA_real_,
      status = "df_resid_nonpositive"
    ))
  }

  rss <- sum(residuals^2)
  sigma2 <- rss / df_resid

  if (!is.finite(sigma2) || sigma2 < 0) {
    return(list(
      stat = 0,
      stat_type = "sigma2_invalid",
      df = NA_real_,
      status = "sigma2_invalid"
    ))
  }

  # Extract the tested block, which is placed after the nuisance block in the cached design.
  tested_coef_index <- seq.int(
    n_nuisance_coef + 1L,
    n_nuisance_coef + n_test_coef
  )

  coef_tested <- coef_estimates[tested_coef_index]
  vcov_tested <- sigma2 * xtx_inverse[tested_coef_index, tested_coef_index, drop = FALSE]

  wald_result <- wald_stat(
    coef_estimates = coef_tested,
    vcov_matrix = vcov_tested
  )

  list(
    stat = unname(wald_result$stat),
    stat_type = wald_result$stat_type,
    df = if (n_test_coef == 1L) df_resid else wald_result$df,
    status = wald_result$status
  )
}

# Winsorize residual tails after centering around the median
winsorize_residuals <- function(
    e,
    quantile = 0.05,
    min_n = 20L
) {
  # Work on a numeric copy so the original residual vector is never mutated by reference.
  residuals <- as.numeric(e)

  finite_rows <- is.finite(residuals)
  residuals_out <- residuals

  # Skip winsorization when too few finite residuals make empirical cutoffs unstable.
  if (sum(finite_rows) < as.integer(min_n)) {
    return(list(
      e_wins = residuals_out,
      lower = NA_real_,
      upper = NA_real_,
      n_wins_low = 0L,
      n_wins_high = 0L,
      applied = FALSE,
      status = "too_few_residuals"
    ))
  }

  if (!is.finite(quantile) || quantile <= 0 || quantile >= 0.5) {
    stop("quantile must be in (0, 0.5).")
  }

  # Center before clipping so asymmetric residual location does not move the cutoffs.
  finite_residuals <- residuals[finite_rows]
  residual_center <- stats::median(finite_residuals)
  centered_residuals <- finite_residuals - residual_center

  lower_cut <- as.numeric(stats::quantile(
    centered_residuals,
    probs = quantile,
    names = FALSE,
    na.rm = TRUE
  ))
  upper_cut <- as.numeric(stats::quantile(
    centered_residuals,
    probs = 1 - quantile,
    names = FALSE,
    na.rm = TRUE
  ))

  # Clip only the empirical tails and then add the residual center back.
  clipped_centered_residuals <- pmin(
    pmax(centered_residuals, lower_cut),
    upper_cut
  )
  winsorized_residuals <- residual_center + clipped_centered_residuals

  n_low <- sum(centered_residuals < lower_cut, na.rm = TRUE)
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
# Cache commit helpers
# -----------------------------------------------------------------------------

# Mark a taxon cache entry as neutralized for null generation
mark_cache_entry_failed <- function(entry, status) {
  entry$observed_stat <- 0
  entry$usable_for_testing <- FALSE
  entry$status <- status
  entry
}

# Commit a taxon cache entry and synchronize the summary vectors
save_cache_entry <- function(
    taxon_entries,
    observed_stat,
    usable_for_testing,
    status,
    taxon_index,
    entry,
    neutral_stat
) {
  taxon_entries[[taxon_index]] <- entry

  observed_stat[taxon_index] <- if (isTRUE(entry$usable_for_testing)) {
    entry$observed_stat
  } else {
    neutral_stat
  }

  usable_for_testing[taxon_index] <- isTRUE(entry$usable_for_testing)
  status[taxon_index] <- entry$status

  list(
    taxon_entries = taxon_entries,
    observed_stat = observed_stat,
    usable_for_testing = usable_for_testing,
    status = status
  )
}

# Commit a failed cache entry using the component-neutral statistic
save_failed_cache_entry <- function(
    taxon_entries,
    observed_stat,
    usable_for_testing,
    status,
    taxon_index,
    entry,
    failure_status,
    neutral_stat
) {
  save_cache_entry(
    taxon_entries = taxon_entries,
    observed_stat = observed_stat,
    usable_for_testing = usable_for_testing,
    status = status,
    taxon_index = taxon_index,
    entry = mark_cache_entry_failed(entry, failure_status),
    neutral_stat = neutral_stat
  )
}

# -----------------------------------------------------------------------------
# Freedman-Lane taxon cache construction
# -----------------------------------------------------------------------------

# Prepare taxon-level fitted values, residuals, and design pieces for FL nulls
prepare_abundance_fl_cache <- function(
    abundance_matrix,
    metadata,
    formula,
    tested_term,
    tested_vars = NULL,
    observed_fit = NULL,
    fit_args = list(),
    neutral_stat = 0,
    verbose = FALSE,
    cache_label = "Abundance"
) {
  # Normalize inputs before any taxon-specific model frames are constructed.
  abundance_matrix <- as_numeric_matrix(abundance_matrix, "abundance_matrix")
  metadata <- align_metadata_to_matrix(metadata, abundance_matrix)

  # Read optional fitting controls from fit_args while keeping stable defaults.
  min_nonzero_n <- if (is.null(fit_args$min_nonzero_n)) {
    3L
  } else {
    as.integer(fit_args$min_nonzero_n)
  }

  resid_winsorize <- isTRUE(fit_args$resid_winsorize)

  resid_winsor_quantile <- if (is.null(fit_args$resid_winsor_quantile)) {
    0.05
  } else {
    as.numeric(fit_args$resid_winsor_quantile)
  }

  resid_winsor_min_n <- if (is.null(fit_args$resid_winsor_min_n)) {
    20L
  } else {
    as.integer(fit_args$resid_winsor_min_n)
  }

  # Attach the synthetic response column and resolve the tested fixed-effect block.
  formula_full <- make_response_formula(
    formula = formula,
    response_name = ".response"
  )

  formula_info <- get_fixed_formula_terms(formula_full, tested_term = tested_term)
  fixed_formula <- formula_info$fixed_formula
  tested_term_idx <- formula_info$tested_term_idx

  # Allocate one cache entry and one summary status per taxon.
  n_taxa <- ncol(abundance_matrix)
  taxa_names <- colnames(abundance_matrix)

  if (is.null(taxa_names)) {
    taxa_names <- paste0("Taxon_", seq_len(n_taxa))
  }

  taxon_entries <- vector("list", n_taxa)
  observed_stat <- rep(neutral_stat, n_taxa)
  usable <- rep(FALSE, n_taxa)
  status <- rep(NA_character_, n_taxa)

  # Build each taxon cache independently so failed taxa can be neutralized cleanly.
  for (taxon_index in seq_len(n_taxa)) {
    taxon_response <- abundance_matrix[, taxon_index]

    # Build the complete-case model frame for this taxon.
    fit_data <- as_plain_data_frame(metadata)
    fit_data$.response <- as.numeric(taxon_response)

    needed_vars <- unique(c(".response", all.vars(fixed_formula)))
    complete_rows <- stats::complete.cases(fit_data[, needed_vars, drop = FALSE])

    taxon_fit_data <- fit_data[complete_rows, , drop = FALSE]
    taxon_fit_data <- droplevels(taxon_fit_data)

    row_idx <- which(complete_rows)
    observed_response <- taxon_fit_data$.response
    n_obs <- nrow(taxon_fit_data)
    group_counts <- tested_term_group_counts(taxon_fit_data, tested_term)

    # Initialize a neutral cache entry and fill diagnostics only after each check passes.
    entry <- list(
      taxon_index = taxon_index,
      taxon_name = taxa_names[taxon_index],
      row_idx = row_idx,
      y = observed_response,
      n_obs = n_obs,
      group_counts = group_counts,
      Z_nuis = NULL,
      X_test = NULL,
      X_full = NULL,
      qr_full = NULL,
      XtX_inv = NULL,
      yhat_red = NULL,
      e_red = NULL,
      p_nuis = NA_integer_,
      p_test = NA_integer_,
      observed_stat = neutral_stat,
      usable_for_testing = FALSE,
      status = NA_character_,
      winsor_applied = FALSE,
      winsor_status = "not_used",
      winsor_lower = NA_real_,
      winsor_upper = NA_real_,
      n_wins_low = 0L,
      n_wins_high = 0L
    )

    # Neutralize taxa without enough observations after taxon-specific missingness.
    if (n_obs < min_nonzero_n) {
      commit <- save_failed_cache_entry(
        taxon_entries = taxon_entries,
        observed_stat = observed_stat,
        usable_for_testing = usable,
        status = status,
        taxon_index = taxon_index,
        entry = entry,
        failure_status = "too_few_nonzero_samples",
        neutral_stat = neutral_stat
      )
      taxon_entries <- commit$taxon_entries
      observed_stat <- commit$observed_stat
      usable <- commit$usable_for_testing
      status <- commit$status
      next
    }

    # Neutralize taxa whose response cannot produce a tested association statistic.
    if (length(unique(observed_response)) < 2L) {
      commit <- save_failed_cache_entry(
        taxon_entries = taxon_entries,
        observed_stat = observed_stat,
        usable_for_testing = usable,
        status = status,
        taxon_index = taxon_index,
        entry = entry,
        failure_status = "constant_response",
        neutral_stat = neutral_stat
      )
      taxon_entries <- commit$taxon_entries
      observed_stat <- commit$observed_stat
      usable <- commit$usable_for_testing
      status <- commit$status
      next
    }

    # Neutralize taxa whose complete-case subset collapses the tested variable.
    if (!tested_term_has_variation(taxon_fit_data, tested_term)) {
      commit <- save_failed_cache_entry(
        taxon_entries = taxon_entries,
        observed_stat = observed_stat,
        usable_for_testing = usable,
        status = status,
        taxon_index = taxon_index,
        entry = entry,
        failure_status = "insufficient_tested_term_variation",
        neutral_stat = neutral_stat
      )
      taxon_entries <- commit$taxon_entries
      observed_stat <- commit$observed_stat
      usable <- commit$usable_for_testing
      status <- commit$status
      next
    }

    # Build the fixed-effect design on the same rows used by the taxon response.
    model_matrix <- tryCatch(
      stats::model.matrix(formula_info$terms_obj, data = taxon_fit_data),
      error = function(e) NULL
    )

    if (is.null(model_matrix)) {
      commit <- save_failed_cache_entry(
        taxon_entries = taxon_entries,
        observed_stat = observed_stat,
        usable_for_testing = usable,
        status = status,
        taxon_index = taxon_index,
        entry = entry,
        failure_status = "model_matrix_failed",
        neutral_stat = neutral_stat
      )
      taxon_entries <- commit$taxon_entries
      observed_stat <- commit$observed_stat
      usable <- commit$usable_for_testing
      status <- commit$status
      next
    }

    # Split the model matrix into nuisance and tested-term columns.
    design_split <- split_tested_and_nuisance_columns(
      model_matrix,
      tested_term_idx = tested_term_idx
    )

    # Neutralize when formula expansion produced no columns for the tested block.
    if (length(design_split$tested_cols) < 1L) {
      commit <- save_failed_cache_entry(
        taxon_entries = taxon_entries,
        observed_stat = observed_stat,
        usable_for_testing = usable,
        status = status,
        taxon_index = taxon_index,
        entry = entry,
        failure_status = "tested_block_missing",
        neutral_stat = neutral_stat
      )
      taxon_entries <- commit$taxon_entries
      observed_stat <- commit$observed_stat
      usable <- commit$usable_for_testing
      status <- commit$status
      next
    }

    nuisance_matrix <- as.matrix(model_matrix[, design_split$nuisance_cols, drop = FALSE])
    tested_matrix <- as.matrix(model_matrix[, design_split$tested_cols, drop = FALSE])

    # Reject tested blocks that are technically present but numerically flat.
    tested_sd <- apply(tested_matrix, 2L, function(x) stats::sd(x, na.rm = TRUE))
    if (all(!is.finite(tested_sd) | tested_sd < 1e-12)) {
      commit <- save_failed_cache_entry(
        taxon_entries = taxon_entries,
        observed_stat = observed_stat,
        usable_for_testing = usable,
        status = status,
        taxon_index = taxon_index,
        entry = entry,
        failure_status = "degenerate_test_block",
        neutral_stat = neutral_stat
      )
      taxon_entries <- commit$taxon_entries
      observed_stat <- commit$observed_stat
      usable <- commit$usable_for_testing
      status <- commit$status
      next
    }

    # Fit the reduced model and cache its fitted values and residuals for FL resampling.
    nuisance_qr <- qr(nuisance_matrix)
    yhat_red <- qr.fitted(nuisance_qr, observed_response)
    e_red_raw <- qr.resid(nuisance_qr, observed_response)

    # Optionally dampen extreme reduced-model residuals before the null is generated.
    winsor_object <- if (isTRUE(resid_winsorize)) {
      winsorize_residuals(
        e = e_red_raw,
        quantile = resid_winsor_quantile,
        min_n = resid_winsor_min_n
      )
    } else {
      list(
        e_wins = e_red_raw,
        lower = NA_real_,
        upper = NA_real_,
        n_wins_low = 0L,
        n_wins_high = 0L,
        applied = FALSE,
        status = "not_used"
      )
    }

    e_red <- winsor_object$e_wins
    y_obs_used <- yhat_red + e_red

    # Precompute the full design so each permutation only changes the response vector.
    full_model_matrix <- cbind(nuisance_matrix, tested_matrix)
    qr_full <- qr(full_model_matrix)

    # A rank-deficient full design cannot support a stable tested-block statistic.
    if (qr_full$rank < ncol(full_model_matrix)) {
      commit <- save_failed_cache_entry(
        taxon_entries = taxon_entries,
        observed_stat = observed_stat,
        usable_for_testing = usable,
        status = status,
        taxon_index = taxon_index,
        entry = entry,
        failure_status = "design_rank_deficient",
        neutral_stat = neutral_stat
      )
      taxon_entries <- commit$taxon_entries
      observed_stat <- commit$observed_stat
      usable <- commit$usable_for_testing
      status <- commit$status
      next
    }

    # Cache the coefficient covariance kernel used by the repeated Wald computations.
    xtx_inverse <- tryCatch(
      chol2inv(qr.R(qr_full)),
      error = function(e) NULL
    )

    if (is.null(xtx_inverse)) {
      commit <- save_failed_cache_entry(
        taxon_entries = taxon_entries,
        observed_stat = observed_stat,
        usable_for_testing = usable,
        status = status,
        taxon_index = taxon_index,
        entry = entry,
        failure_status = "xtx_inverse_failed",
        neutral_stat = neutral_stat
      )
      taxon_entries <- commit$taxon_entries
      observed_stat <- commit$observed_stat
      usable <- commit$usable_for_testing
      status <- commit$status
      next
    }

    # Compute the observed statistic on the same residualized response used by the null.
    observed_test <- linear_block_wald_from_precomp(
      y = y_obs_used,
      qr_full = qr_full,
      xtx_inverse = xtx_inverse,
      n_nuisance_coef = ncol(nuisance_matrix),
      n_test_coef = ncol(tested_matrix)
    )

    # Store all pieces needed by one-permutation kernels and downstream diagnostics.
    entry$Z_nuis <- nuisance_matrix
    entry$X_test <- tested_matrix
    entry$X_full <- full_model_matrix
    entry$qr_full <- qr_full
    entry$XtX_inv <- xtx_inverse
    entry$yhat_red <- yhat_red
    entry$e_red_raw <- e_red_raw
    entry$e_red <- e_red
    entry$winsor_applied <- isTRUE(winsor_object$applied)
    entry$winsor_status <- winsor_object$status
    entry$winsor_lower <- winsor_object$lower
    entry$winsor_upper <- winsor_object$upper
    entry$n_wins_low <- winsor_object$n_wins_low
    entry$n_wins_high <- winsor_object$n_wins_high
    entry$p_nuis <- ncol(nuisance_matrix)
    entry$p_test <- ncol(tested_matrix)
    entry$observed_stat <- observed_test$stat
    entry$usable_for_testing <- identical(observed_test$status, "ok")
    entry$status <- observed_test$status

    # Commit successful and unsuccessful observed tests through the same summary path.
    commit <- save_cache_entry(
      taxon_entries = taxon_entries,
      observed_stat = observed_stat,
      usable_for_testing = usable,
      status = status,
      taxon_index = taxon_index,
      entry = entry,
      neutral_stat = neutral_stat
    )
    taxon_entries <- commit$taxon_entries
    observed_stat <- commit$observed_stat
    usable <- commit$usable_for_testing
    status <- commit$status
  }

  # Return the cache as the observed-statistic source and null-generation input.
  out <- list(
    taxon_entries = taxon_entries,
    taxon_name = taxa_names,
    observed_stat = observed_stat,
    usable_for_testing = usable,
    status = status,
    n_taxa = n_taxa,
    neutral_stat = neutral_stat,
    tested_term = tested_term,
    tested_vars = tested_vars,
    statistic_type = "fl_block_wald",
    resid_winsorize = resid_winsorize,
    resid_winsor_quantile = resid_winsor_quantile,
    resid_winsor_min_n = resid_winsor_min_n
  )

  class(out) <- c("pursue_abundance_fl_cache", "list")

  if (verbose) {
    message("Prepared ", cache_label, " FL cache.")
    message("  Taxa: ", n_taxa)
    message("  Usable for testing: ", sum(usable))
    message("  Too few nonzero: ", sum(status == "too_few_nonzero_samples", na.rm = TRUE))
    message("  Degenerate tested block: ", sum(status == "degenerate_test_block", na.rm = TRUE))
    if (isTRUE(resid_winsorize)) {
      message("  Taxa with winsorization applied: ", sum(vapply(taxon_entries, function(x) isTRUE(x$winsor_applied), logical(1))))
    }
  }

  out
}

# -----------------------------------------------------------------------------
# Permutation object and row-context helpers
# -----------------------------------------------------------------------------

# Extract the statistic matrix from raw or wrapped null results
extract_null_stat_matrix <- function(x) {
  if (is.matrix(x)) {
    return(x)
  }

  if (is.list(x) && !is.null(x$stat_mat)) {
    return(x$stat_mat)
  }

  stop("Cannot extract stat_mat from object.")
}

# Extract the full-sample cluster vector from a permutation design
extract_full_cluster_id <- function(perm_design) {
  if (!inherits(perm_design, "pursue_permutation_design")) {
    stop("perm_design must inherit from 'pursue_permutation_design'.")
  }

  cluster_id <- perm_design$cluster_id
  if (is.null(cluster_id)) {
    return(NULL)
  }

  cluster_id <- as.character(cluster_id)
  if (length(cluster_id) != perm_design$n) {
    stop("perm_design$cluster_id has wrong length.")
  }

  cluster_id
}

# Convert sample or cluster permutation objects into a named unit map
normalize_permutation_unit_map <- function(perm_obj, perm_design) {
  if (!inherits(perm_design, "pursue_permutation_design")) {
    stop("perm_design must inherit from 'pursue_permutation_design'.")
  }

  if (!is.list(perm_obj) || is.null(perm_obj$type)) {
    stop("perm_obj must be a permutation object from generate_permutation_indices().")
  }

  # Sample permutations map target sample IDs to donor sample IDs.
  if (identical(perm_obj$type, "sample")) {
    sample_index <- perm_obj$index

    if (is.null(sample_index)) {
      stop("Sample permutation object is missing index.")
    }

    if (length(sample_index) != perm_design$n) {
      stop("Sample permutation index has wrong length.")
    }

    unit_ids <- as.character(seq_len(perm_design$n))
    unit_map <- as.character(sample_index)
    names(unit_map) <- unit_ids

    return(list(
      unit_type = "sample",
      unit_ids = unit_ids,
      unit_map = unit_map
    ))
  }

  # Cluster permutations map target cluster IDs to donor cluster IDs.
  if (identical(perm_obj$type, "cluster")) {
    cluster_map <- perm_obj$cluster_map

    if (is.null(cluster_map)) {
      stop("Cluster permutation object is missing cluster_map.")
    }

    if (is.null(names(cluster_map))) {
      cluster_levels <- perm_design$cluster_levels

      if (is.null(cluster_levels)) {
        stop("Unnamed cluster_map requires perm_design$cluster_levels.")
      }

      if (length(cluster_map) != length(cluster_levels)) {
        stop("Unnamed cluster_map has wrong length.")
      }

      names(cluster_map) <- as.character(cluster_levels)
    }

    unit_ids <- as.character(names(cluster_map))
    unit_map <- as.character(cluster_map)

    return(list(
      unit_type = "cluster",
      unit_ids = unit_ids,
      unit_map = unit_map
    ))
  }

  stop("Unknown perm_obj$type: ", perm_obj$type)
}

# Map taxon-valid sample rows to the active permutation unit IDs
map_rows_to_unit_ids <- function(row_idx, perm_design, unit_level = perm_design$unit_level) {
  if (!inherits(perm_design, "pursue_permutation_design")) {
    stop("perm_design must inherit from 'pursue_permutation_design'.")
  }

  row_idx <- as.integer(row_idx)

  if (length(row_idx) == 0L) {
    return(character(0))
  }

  if (any(!is.finite(row_idx)) || any(row_idx < 1L) || any(row_idx > perm_design$n)) {
    stop("row_idx contains invalid sample rows.")
  }

  unit_level <- match.arg(unit_level, c("sample", "cluster"))

  # For sample-level designs, each valid row is its own exchangeability unit.
  if (identical(unit_level, "sample")) {
    return(as.character(row_idx))
  }

  # For cluster-level designs, each valid row inherits its full-sample cluster label.
  cluster_id <- extract_full_cluster_id(perm_design)
  if (is.null(cluster_id)) {
    stop("cluster_id is required for cluster-level row-to-unit mapping.")
  }

  as.character(cluster_id[row_idx])
}

# Store row-level sample, cluster, and stratum labels for one taxon
build_taxon_row_context <- function(row_idx, perm_design) {
  if (!inherits(perm_design, "pursue_permutation_design")) {
    stop("perm_design must inherit from 'pursue_permutation_design'.")
  }

  row_idx <- as.integer(row_idx)

  if (length(row_idx) == 0L) {
    return(list(
      row_idx = integer(0),
      sample_id_row = character(0),
      cluster_id_row = NULL,
      strata_row = NULL,
      unit_id_row = character(0)
    ))
  }

  # Keep labels in full-sample coordinates so local taxon subsets can be permuted safely.
  sample_id_row <- as.character(row_idx)

  cluster_id <- extract_full_cluster_id(perm_design)
  cluster_id_row <- if (is.null(cluster_id)) {
    NULL
  } else {
    as.character(cluster_id[row_idx])
  }

  strata_vector <- perm_design$strata
  strata_row <- if (is.null(strata_vector)) {
    NULL
  } else {
    as.character(strata_vector[row_idx])
  }

  unit_id_row <- map_rows_to_unit_ids(
    row_idx = row_idx,
    perm_design = perm_design
  )

  list(
    row_idx = row_idx,
    sample_id_row = sample_id_row,
    cluster_id_row = cluster_id_row,
    strata_row = strata_row,
    unit_id_row = unit_id_row
  )
}

# Summarize valid rows by permutation unit for one taxon
build_taxon_unit_context <- function(row_idx, perm_design) {
  # Group valid local row positions by their active permutation unit.
  row_context <- build_taxon_row_context(
    row_idx = row_idx,
    perm_design = perm_design
  )

  row_pos_by_unit <- split(seq_along(row_context$row_idx), row_context$unit_id_row)

  list(
    unit_level = perm_design$unit_level,
    unit_id_row = row_context$unit_id_row,
    row_pos_by_unit = row_pos_by_unit,
    unit_levels_local = names(row_pos_by_unit)
  )
}

# Apply a cluster permutation to taxon residuals while preserving local block sizes
apply_cluster_map_to_residual_blocks <- function(
    e_red_full,
    row_context,
    perm_obj,
    perm_design,
    strict = TRUE
) {
  if (!inherits(perm_design, "pursue_permutation_design")) {
    stop("perm_design must inherit from 'pursue_permutation_design'.")
  }

  if (!is.numeric(e_red_full) || length(e_red_full) != perm_design$n) {
    stop("e_red_full must be a numeric vector of length perm_design$n.")
  }

  if (!is.list(row_context) || is.null(row_context$row_idx)) {
    stop("row_context must come from build_taxon_row_context().")
  }

  if (!identical(perm_obj$type, "cluster")) {
    stop("apply_cluster_map_to_residual_blocks() requires perm_obj$type == 'cluster'.")
  }

  # Normalize the cluster map so target and donor cluster labels are explicit.
  unit_info <- normalize_permutation_unit_map(perm_obj, perm_design)
  if (!identical(unit_info$unit_type, "cluster")) {
    stop("Normalized unit type is not 'cluster'.")
  }

  row_idx <- row_context$row_idx
  cluster_id_row <- row_context$cluster_id_row

  if (is.null(cluster_id_row)) {
    stop("row_context$cluster_id_row is required for cluster residual-block permutation.")
  }

  # Fill target positions cluster by cluster from donor residual blocks.
  permuted_residuals <- rep(NA_real_, length(row_idx))
  target_rows_by_cluster <- split(seq_along(row_idx), cluster_id_row)

  for (target_cluster in names(target_rows_by_cluster)) {
    # The map is target cluster -> donor cluster for residual transfer.
    donor_cluster <- unname(unit_info$unit_map[target_cluster])

    if (length(donor_cluster) != 1L || is.na(donor_cluster) || !nzchar(donor_cluster)) {
      if (strict) {
        stop("Missing donor cluster for target cluster: ", target_cluster)
      }

      return(list(
        ok = FALSE,
        reason = paste0("missing_donor:", target_cluster),
        values = NULL
      ))
    }

    target_pos <- target_rows_by_cluster[[target_cluster]]
    donor_pos <- which(cluster_id_row == donor_cluster)

    # Local valid-row counts must match or the residual block cannot be transplanted.
    if (length(target_pos) != length(donor_pos)) {
      if (strict) {
        stop(
          "Local valid-row cluster size mismatch: target ", target_cluster,
          " has ", length(target_pos),
          " rows, donor ", donor_cluster,
          " has ", length(donor_pos), " rows."
        )
      }

      return(list(
        ok = FALSE,
        reason = paste0(
          "local_cluster_size_mismatch:",
          target_cluster, "->", donor_cluster,
          " (", length(target_pos), " vs ", length(donor_pos), ")"
        ),
        values = NULL
      ))
    }

    # Pull donor residuals from the full-length vector but write them into local target slots.
    donor_rows_full <- row_idx[donor_pos]
    permuted_residuals[target_pos] <- e_red_full[donor_rows_full]
  }

  if (strict) {
    return(permuted_residuals)
  }

  list(ok = TRUE, reason = "ok", values = permuted_residuals)
}

# -----------------------------------------------------------------------------
# Full-length permutation cache layout
# -----------------------------------------------------------------------------

# Add full-sample residual vectors and row contexts to the FL cache
prepare_abundance_permutation_cache <- function(
    abundance_matrix,
    metadata,
    formula,
    tested_term,
    tested_vars = NULL,
    observed_fit = NULL,
    fit_args = list(),
    neutral_stat = 0,
    perm_design,
    verbose = FALSE,
    cache_label = "Abundance"
) {
  if (missing(perm_design) || is.null(perm_design)) {
    stop("prepare_abundance_permutation_cache() requires perm_design.")
  }

  if (!inherits(perm_design, "pursue_permutation_design")) {
    stop("perm_design must inherit from 'pursue_permutation_design'.")
  }

  # First build the statistical FL cache without assuming a permutation unit layout.
  base_cache <- prepare_abundance_fl_cache(
    abundance_matrix = abundance_matrix,
    metadata = metadata,
    formula = formula,
    tested_term = tested_term,
    tested_vars = tested_vars,
    observed_fit = observed_fit,
    fit_args = fit_args,
    neutral_stat = neutral_stat,
    verbose = verbose,
    cache_label = cache_label
  )

  n_full <- nrow(abundance_matrix)
  if (n_full != perm_design$n) {
    stop("nrow(abundance_matrix) must equal perm_design$n.")
  }

  # Expand each taxon entry to full-sample coordinates for row-safe permutation.
  for (taxon_index in seq_len(base_cache$n_taxa)) {
    entry <- base_cache$taxon_entries[[taxon_index]]

    # Keep a full-length validity mask because taxa can have different complete-case rows.
    valid_mask_full <- rep(FALSE, n_full)
    valid_mask_full[entry$row_idx] <- TRUE

    yhat_red_full <- rep(NA_real_, n_full)
    e_red_full <- rep(NA_real_, n_full)
    y_full <- rep(NA_real_, n_full)

    if (!is.null(entry$yhat_red)) {
      yhat_red_full[entry$row_idx] <- entry$yhat_red
    }

    if (!is.null(entry$e_red)) {
      e_red_full[entry$row_idx] <- entry$e_red
    }

    if (!is.null(entry$y)) {
      y_full[entry$row_idx] <- entry$y
    }

    # Attach row and unit contexts derived from this taxon's valid sample subset.
    entry$valid_mask_full <- valid_mask_full
    entry$row_context <- build_taxon_row_context(
      row_idx = entry$row_idx,
      perm_design = perm_design
    )
    entry$unit_context <- build_taxon_unit_context(
      row_idx = entry$row_idx,
      perm_design = perm_design
    )

    entry$y_full <- y_full
    entry$yhat_red_full <- yhat_red_full
    entry$e_red_full <- e_red_full
    entry$perm_unit_level <- perm_design$unit_level

    base_cache$taxon_entries[[taxon_index]] <- entry
  }

  # Mark the cache as permutation-ready while retaining the original FL cache class.
  base_cache$perm_design_unit_level <- perm_design$unit_level
  base_cache$perm_design_n <- perm_design$n
  base_cache$cache_layout <- "full_length"

  class(base_cache) <- c("pursue_abundance_permutation_cache", class(base_cache))
  base_cache
}

# Induce a taxon-local row order from a global sample permutation
induce_local_sample_permutation <- function(row_idx, perm_obj) {
  if (!is.list(perm_obj) || is.null(perm_obj$type) || perm_obj$type != "sample") {
    stop("Local sample permutation currently supports only sample-level permutation objects.")
  }

  global_index <- perm_obj$index

  if (length(global_index) < max(row_idx)) {
    stop("Permutation index is shorter than max(row_idx).")
  }

  # Convert the global permutation into the ordering induced on this taxon's valid rows.
  position_in_global <- match(row_idx, global_index)

  if (anyNA(position_in_global)) {
    stop("Failed to induce a subset permutation from the global sample permutation.")
  }

  order(position_in_global)
}

# -----------------------------------------------------------------------------
# Sample-level Freedman-Lane null kernels
# -----------------------------------------------------------------------------

# Compute one sample-level FL null column across all taxa
compute_abundance_fl_null_sample_one_perm <- function(
    perm_pos,
    abund_cache,
    perm_list,
    perm_design,
    neutral_stat = 0
) {
  permutation <- perm_list[[perm_pos]]

  if (!identical(permutation$type, "sample")) {
    stop("Sample null kernel received non-sample permutation object.")
  }

  # Initialize all taxa as neutralized; usable taxa overwrite their slot after a valid test.
  n_taxa <- abund_cache$n_taxa
  stat_col <- rep(neutral_stat, n_taxa)
  status_col <- rep("neutralized", n_taxa)

  for (taxon_index in seq_len(n_taxa)) {
    entry <- abund_cache$taxon_entries[[taxon_index]]

    # Skip taxa that already failed observed-cache validation.
    if (!isTRUE(entry$usable_for_testing)) {
      status_col[taxon_index] <- "unusable_component"
      next
    }

    if (any(!is.finite(entry$yhat_red)) || any(!is.finite(entry$e_red))) {
      status_col[taxon_index] <- "nonfinite_cached_residual_or_fit"
      next
    }

    # Restrict the global sample permutation to this taxon's complete-case rows.
    local_index <- induce_local_sample_permutation(
      row_idx = entry$row_idx,
      perm_obj = permutation
    )

    # Freedman-Lane response: reduced fitted values plus permuted reduced residuals.
    permuted_residuals <- entry$e_red[local_index]
    y_star <- entry$yhat_red + permuted_residuals

    if (any(!is.finite(y_star))) {
      status_col[taxon_index] <- "nonfinite_y_star"
      next
    }

    # Refit only the response against the cached full design and score the tested block.
    permuted_test <- linear_block_wald_from_precomp(
      y = y_star,
      qr_full = entry$qr_full,
      xtx_inverse = entry$XtX_inv,
      n_nuisance_coef = entry$p_nuis,
      n_test_coef = entry$p_test
    )

    if (identical(permuted_test$status, "ok")) {
      stat_col[taxon_index] <- abs(permuted_test$stat)
      status_col[taxon_index] <- "ok"
    } else {
      status_col[taxon_index] <- permuted_test$status
    }
  }

  list(b = perm_pos, stat = stat_col, status = status_col)
}

# Compute the full sample-level FL null statistic matrix
compute_abundance_fl_null_sample <- function(
    abund_cache,
    perm_list,
    perm_design,
    neutral_stat = 0,
    parallelize_permutations = FALSE,
    permutation_n_cores = max(1L, parallel::detectCores() - 1L),
    verbose = FALSE,
    null_label = "Abundance"
) {
  if (!inherits(abund_cache, "pursue_abundance_permutation_cache")) {
    stop("abund_cache must come from prepare_abundance_permutation_cache().")
  }

  if (!inherits(perm_design, "pursue_permutation_design")) {
    stop("perm_design must come from build_permutation_design().")
  }

  if (!identical(perm_design$unit_level, "sample")) {
    stop("compute_abundance_fl_null_sample() requires perm_design$unit_level == 'sample'.")
  }

  n_taxa <- abund_cache$n_taxa
  n_perm <- length(perm_list)

  stat_mat <- matrix(
    neutral_stat,
    nrow = n_taxa,
    ncol = n_perm,
    dimnames = list(abund_cache$taxon_name, paste0("perm_", seq_len(n_perm)))
  )

  status_mat <- matrix(
    "neutralized",
    nrow = n_taxa,
    ncol = n_perm,
    dimnames = dimnames(stat_mat)
  )

  # Parallelize over permutation columns only; each worker computes all taxa for one column.
  use_parallel <- isTRUE(parallelize_permutations) &&
    is.numeric(permutation_n_cores) &&
    is.finite(permutation_n_cores) &&
    as.integer(permutation_n_cores) > 1L &&
    n_perm > 1L

  # Export only the cache, permutation list, and kernel helpers needed by PSOCK workers.
  if (use_parallel) {
    n_workers <- min(as.integer(permutation_n_cores), n_perm)

    if (verbose) {
      message(null_label, " FL null (sample): parallelizing ", n_perm,
              " permutations across ", n_workers, " workers.")
    }

    cl <- parallel::makeCluster(n_workers, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterExport(
      cl,
      varlist = c(
        "abund_cache", "perm_list", "perm_design", "neutral_stat",
        "compute_abundance_fl_null_sample_one_perm",
        "induce_local_sample_permutation",
        "linear_block_wald_from_precomp", "wald_stat"
      ),
      envir = environment()
    )

    results <- parallel::parLapplyLB(
      cl,
      X = seq_len(n_perm),
      fun = function(perm_pos) {
        compute_abundance_fl_null_sample_one_perm(
          perm_pos,
          abund_cache = abund_cache,
          perm_list = perm_list,
          perm_design = perm_design,
          neutral_stat = neutral_stat
        )
      }
    )
  } else {
    results <- lapply(seq_len(n_perm), function(perm_pos) {
      compute_abundance_fl_null_sample_one_perm(
        perm_pos,
        abund_cache = abund_cache,
        perm_list = perm_list,
        perm_design = perm_design,
        neutral_stat = neutral_stat
      )
    })
  }

  # Reassemble worker results into taxa-by-permutation matrices.
  for (result in results) {
    stat_mat[, result$b] <- result$stat
    status_mat[, result$b] <- result$status
  }

  if (verbose && !use_parallel && n_perm > 0L) {
    message(null_label, " FL null (sample): completed ", n_perm, " permutations.")
  }

  # Return both null statistics and per-cell status diagnostics.
  out <- list(
    stat_mat = stat_mat,
    status_mat = status_mat,
    diagnostics = list(
      unit_level = "sample",
      n_taxa = n_taxa,
      n_perm = n_perm,
      n_ok = sum(status_mat == "ok", na.rm = TRUE),
      n_neutralized = sum(status_mat != "ok", na.rm = TRUE)
    )
  )

  class(out) <- c("pursue_abundance_null", "pursue_component_null")
  out
}

# -----------------------------------------------------------------------------
# Cluster-level Freedman-Lane null kernels
# -----------------------------------------------------------------------------

# Compute one cluster-level FL null column across all taxa
compute_abundance_fl_null_cluster_one_perm <- function(
    perm_pos,
    abund_cache,
    perm_list,
    perm_design,
    neutral_stat = 0
) {
  permutation <- perm_list[[perm_pos]]

  # Require a named cluster map before any taxon residual blocks are moved.
  unit_info <- normalize_permutation_unit_map(permutation, perm_design)
  if (!identical(unit_info$unit_type, "cluster")) {
    stop("Cluster null kernel received non-cluster permutation object.")
  }

  # Initialize all taxa as neutralized; usable taxa overwrite their slot after a valid test.
  n_taxa <- abund_cache$n_taxa
  stat_col <- rep(neutral_stat, n_taxa)
  status_col <- rep("neutralized", n_taxa)

  for (taxon_index in seq_len(n_taxa)) {
    entry <- abund_cache$taxon_entries[[taxon_index]]

    # Skip taxa that already failed observed-cache validation.
    if (!isTRUE(entry$usable_for_testing)) {
      status_col[taxon_index] <- "unusable_component"
      next
    }

    # Transfer reduced residuals as cluster blocks inside this taxon's valid-row subset.
    permuted_residual_object <- apply_cluster_map_to_residual_blocks(
      e_red_full = entry$e_red_full,
      row_context = entry$row_context,
      perm_obj = permutation,
      perm_design = perm_design,
      strict = FALSE
    )

    if (!isTRUE(permuted_residual_object$ok)) {
      status_col[taxon_index] <- permuted_residual_object$reason
      next
    }

    # Freedman-Lane response: reduced fitted values plus cluster-permuted residual blocks.
    permuted_residuals <- permuted_residual_object$values
    y_star <- entry$yhat_red + permuted_residuals

    if (any(!is.finite(y_star))) {
      status_col[taxon_index] <- "nonfinite_y_star"
      next
    }

    # Refit only the response against the cached full design and score the tested block.
    permuted_test <- linear_block_wald_from_precomp(
      y = y_star,
      qr_full = entry$qr_full,
      xtx_inverse = entry$XtX_inv,
      n_nuisance_coef = entry$p_nuis,
      n_test_coef = entry$p_test
    )

    if (identical(permuted_test$status, "ok")) {
      stat_col[taxon_index] <- abs(permuted_test$stat)
      status_col[taxon_index] <- "ok"
    } else {
      status_col[taxon_index] <- permuted_test$status
    }
  }

  list(b = perm_pos, stat = stat_col, status = status_col)
}

# Compute the full cluster-level FL null statistic matrix
compute_abundance_fl_null_cluster <- function(
    abund_cache,
    perm_list,
    perm_design,
    neutral_stat = 0,
    parallelize_permutations = FALSE,
    permutation_n_cores = max(1L, parallel::detectCores() - 1L),
    verbose = FALSE,
    null_label = "Abundance"
) {
  if (!inherits(abund_cache, "pursue_abundance_permutation_cache")) {
    stop("abund_cache must come from prepare_abundance_permutation_cache().")
  }

  if (!inherits(perm_design, "pursue_permutation_design")) {
    stop("perm_design must come from build_permutation_design().")
  }

  if (!identical(perm_design$unit_level, "cluster")) {
    stop("compute_abundance_fl_null_cluster() requires perm_design$unit_level == 'cluster'.")
  }

  n_taxa <- abund_cache$n_taxa
  n_perm <- length(perm_list)

  stat_mat <- matrix(
    neutral_stat,
    nrow = n_taxa,
    ncol = n_perm,
    dimnames = list(abund_cache$taxon_name, paste0("perm_", seq_len(n_perm)))
  )

  status_mat <- matrix(
    "neutralized",
    nrow = n_taxa,
    ncol = n_perm,
    dimnames = dimnames(stat_mat)
  )

  # Parallelize over permutation columns only; each worker computes all taxa for one column.
  use_parallel <- isTRUE(parallelize_permutations) &&
    is.numeric(permutation_n_cores) &&
    is.finite(permutation_n_cores) &&
    as.integer(permutation_n_cores) > 1L &&
    n_perm > 1L

  # Export only the cache, permutation list, and kernel helpers needed by PSOCK workers.
  if (use_parallel) {
    n_workers <- min(as.integer(permutation_n_cores), n_perm)

    if (verbose) {
      message(null_label, " FL null (cluster): parallelizing ", n_perm,
              " permutations across ", n_workers, " workers.")
    }

    cl <- parallel::makeCluster(n_workers, type = "PSOCK")
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterExport(
      cl,
      varlist = c(
        "abund_cache", "perm_list", "perm_design", "neutral_stat",
        "compute_abundance_fl_null_cluster_one_perm",
        "normalize_permutation_unit_map",
        "apply_cluster_map_to_residual_blocks",
        "linear_block_wald_from_precomp", "wald_stat"
      ),
      envir = environment()
    )

    results <- parallel::parLapplyLB(
      cl,
      X = seq_len(n_perm),
      fun = function(perm_pos) {
        compute_abundance_fl_null_cluster_one_perm(
          perm_pos,
          abund_cache = abund_cache,
          perm_list = perm_list,
          perm_design = perm_design,
          neutral_stat = neutral_stat
        )
      }
    )
  } else {
    results <- lapply(seq_len(n_perm), function(perm_pos) {
      compute_abundance_fl_null_cluster_one_perm(
        perm_pos,
        abund_cache = abund_cache,
        perm_list = perm_list,
        perm_design = perm_design,
        neutral_stat = neutral_stat
      )
    })
  }

  # Reassemble worker results into taxa-by-permutation matrices.
  for (result in results) {
    stat_mat[, result$b] <- result$stat
    status_mat[, result$b] <- result$status
  }

  if (verbose && !use_parallel && n_perm > 0L) {
    message(null_label, " FL null (cluster): completed ", n_perm, " permutations.")
  }

  # Return both null statistics and per-cell status diagnostics.
  out <- list(
    stat_mat = stat_mat,
    status_mat = status_mat,
    diagnostics = list(
      unit_level = "cluster",
      n_taxa = n_taxa,
      n_perm = n_perm,
      n_ok = sum(status_mat == "ok", na.rm = TRUE),
      n_neutralized = sum(status_mat != "ok", na.rm = TRUE),
      neutralization_reasons = sort(table(status_mat[status_mat != "ok"]), decreasing = TRUE)
    )
  )

  class(out) <- c("pursue_abundance_null", "pursue_component_null")
  out
}

# -----------------------------------------------------------------------------
# Null dispatcher
# -----------------------------------------------------------------------------

# Dispatch abundance null generation to the active permutation unit level
compute_abundance_permutation_null <- function(
    abund_cache,
    perm_list,
    perm_design,
    neutral_stat = 0,
    parallelize_permutations = FALSE,
    permutation_n_cores = max(1L, parallel::detectCores() - 1L),
    verbose = FALSE,
    null_label = "Abundance"
) {
  if (identical(perm_design$unit_level, "sample")) {
    return(compute_abundance_fl_null_sample(
      abund_cache = abund_cache,
      perm_list = perm_list,
      perm_design = perm_design,
      neutral_stat = neutral_stat,
      parallelize_permutations = parallelize_permutations,
      permutation_n_cores = permutation_n_cores,
      verbose = verbose,
      null_label = null_label
    ))
  }

  if (identical(perm_design$unit_level, "cluster")) {
    return(compute_abundance_fl_null_cluster(
      abund_cache = abund_cache,
      perm_list = perm_list,
      perm_design = perm_design,
      neutral_stat = neutral_stat,
      parallelize_permutations = parallelize_permutations,
      permutation_n_cores = permutation_n_cores,
      verbose = verbose,
      null_label = null_label
    ))
  }

  stop("Unsupported perm_design$unit_level.")
}
