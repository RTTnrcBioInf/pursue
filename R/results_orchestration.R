# -----------------------------------------------------------------------------
# Component permutation orchestrator
# -----------------------------------------------------------------------------

# Fit observed components, generate shared nulls, and assemble permutation diagnostics
run_component_permutations <- function(
    prevalence_matrix,
    abundance_matrix,
    metadata,
    formula,
    tested_term,
    prevalence_formula = NULL,
    tested_vars = NULL,
    observed_prevalence_fit = NULL,
    observed_abundance_fit = NULL,
    n_perm = 199,
    cluster_id = NULL,
    strata = NULL,
    cluster_exposure_mode = c("auto", "sample", "cluster"),
    prevalence_fit_args = list(),
    abundance_fit_args = list(),
    neutral_stat = 0,
    parallelize_permutations = FALSE,
    permutation_n_cores = max(1L, parallel::detectCores() - 1L),
    seed = NULL,
    verbose = FALSE
) {
  n_perm <- as.integer(n_perm)

  if (!is.finite(n_perm) || n_perm < 1L) {
    stop("n_perm must be a positive integer.")
  }

  metadata <- as_plain_data_frame(metadata)
  cluster_exposure_mode <- match.arg(cluster_exposure_mode)

  # Resolve the tested variables and prevalence formula before any component fitting.
  tested_vars <- check_tested_vars(metadata, tested_vars, tested_term)
  prevalence_formula_use <- if (is.null(prevalence_formula)) {
    formula
  } else {
    prevalence_formula
  }

  # Fit the observed prevalence arm unless the caller supplied a reusable fit.
  if (is.null(observed_prevalence_fit)) {
    prevalence_fit_call <- list(
      prevalence_matrix = prevalence_matrix,
      metadata = metadata,
      formula = prevalence_formula_use,
      tested_term = tested_term
    )

    observed_prevalence_fit <- do.call(
      fit_prevalence_expected_rarefied_arm,
      c(prevalence_fit_call, prevalence_fit_args)
    )
  }

  # Fit the observed abundance arm unless the caller supplied a reusable fit.
  if (is.null(observed_abundance_fit)) {
    abundance_fit_call <- list(
      abundance_matrix = abundance_matrix,
      metadata = metadata,
      formula = formula,
      tested_term = tested_term
    )

    observed_abundance_fit <- do.call(
      fit_abundance_arm,
      c(abundance_fit_call, abundance_fit_args)
    )
  }

  # Join observed component statistics before attaching permutation-specific metadata.
  observed_component_table <- combine_component_tables(
    prevalence_fit = observed_prevalence_fit,
    abundance_fit = observed_abundance_fit,
    neutral_stat = neutral_stat
  )

  # Build one exchangeability design shared by both model components.
  perm_design <- build_permutation_design(
    metadata = metadata,
    tested_vars = tested_vars,
    tested_term = tested_term,
    cluster_id = cluster_id,
    strata = strata,
    cluster_exposure_mode = cluster_exposure_mode,
    verbose = FALSE
  )

  # Draw the concrete permutation maps used by both prevalence and abundance nulls.
  perm_list <- generate_permutation_indices(
    perm_design = perm_design,
    n_perm = n_perm,
    seed = seed,
    include_identity = FALSE,
    verbose = FALSE
  )

  # Prepare the prevalence FL cache and compute its null statistic matrix.
  prev_cache <- prepare_expected_rarefied_prevalence_cache(
    prevalence_matrix = prevalence_matrix,
    metadata = metadata,
    formula = prevalence_formula_use,
    tested_term = tested_term,
    tested_vars = tested_vars,
    observed_fit = observed_prevalence_fit,
    fit_args = prevalence_fit_args,
    neutral_stat = neutral_stat,
    perm_design = perm_design,
    verbose = verbose
  )

  prev_null <- compute_expected_rarefied_prevalence_null(
    prev_cache = prev_cache,
    perm_list = perm_list,
    perm_design = perm_design,
    neutral_stat = neutral_stat,
    parallelize_permutations = parallelize_permutations,
    permutation_n_cores = permutation_n_cores,
    verbose = verbose
  )

  prev_perm_stats <- extract_null_stat_matrix(prev_null)

  if (verbose) message(strrep("-", 50))

  # Prepare the abundance FL cache and compute its null statistic matrix.
  abund_cache <- prepare_abundance_permutation_cache(
    abundance_matrix = abundance_matrix,
    metadata = metadata,
    formula = formula,
    tested_term = tested_term,
    tested_vars = tested_vars,
    observed_fit = observed_abundance_fit,
    fit_args = abundance_fit_args,
    neutral_stat = neutral_stat,
    perm_design = perm_design,
    verbose = verbose,
    cache_label = "Abundance"
  )

  abund_null <- compute_abundance_permutation_null(
    abund_cache = abund_cache,
    perm_list = perm_list,
    perm_design = perm_design,
    neutral_stat = neutral_stat,
    parallelize_permutations = parallelize_permutations,
    permutation_n_cores = permutation_n_cores,
    verbose = verbose,
    null_label = "Abundance"
  )

  abund_perm_stats <- extract_null_stat_matrix(abund_null)

  if (verbose) message(strrep("-", 50))

  # Attach cache-derived usability and observed-statistic diagnostics to the component table.
  observed_component_table$prev_usable_for_testing <- prev_cache$usable_for_testing
  observed_component_table$prev_status <- prev_cache$status
  observed_component_table$prev_stat_for_permutation <- prev_cache$observed_stat
  observed_component_table$prev_reason_unusable <- ifelse(
    prev_cache$usable_for_testing,
    NA_character_,
    prev_cache$status
  )

  observed_component_table$abund_usable_for_testing <- abund_cache$usable_for_testing
  observed_component_table$abund_status <- abund_cache$status
  observed_component_table$abund_stat_for_permutation <- abund_cache$observed_stat
  observed_component_table$abund_reason_unusable <- ifelse(
    abund_cache$usable_for_testing,
    NA_character_,
    abund_cache$status
  )

  # Empirical p-values use two-sided absolute component statistics.
  observed_prev_stat <- abs(prev_cache$observed_stat)
  observed_abund_stat <- abs(abund_cache$observed_stat)

  # Assert that both null engines returned matrices aligned to the observed taxon table.
  n_taxa <- nrow(observed_component_table)
  expected_dim <- c(as.integer(n_taxa), as.integer(n_perm))

  if (!is.matrix(prev_perm_stats) || !all(dim(prev_perm_stats) == expected_dim)) {
    stop(
      "Prevalence null engine must return a matrix of dimension n_taxa x n_perm. ",
      "Got: ", paste(dim(prev_perm_stats), collapse = " x "),
      "; expected: ", paste(expected_dim, collapse = " x ")
    )
  }

  if (!is.matrix(abund_perm_stats) || !all(dim(abund_perm_stats) == expected_dim)) {
    stop(
      "Abundance null engine must return a matrix of dimension n_taxa x n_perm. ",
      "Got: ", paste(dim(abund_perm_stats), collapse = " x "),
      "; expected: ", paste(expected_dim, collapse = " x ")
    )
  }

  # Keep observed-table order as the source of truth even when returned row names disagree.
  if (!is.null(rownames(prev_perm_stats))) {
    if (!identical(rownames(prev_perm_stats), observed_component_table$taxon_name)) {
      warning("Row names of prev_perm_stats do not match observed taxon order; using observed order as source of truth.")
    }
  }

  if (!is.null(rownames(abund_perm_stats))) {
    if (!identical(rownames(abund_perm_stats), observed_component_table$taxon_name)) {
      warning("Row names of abund_perm_stats do not match observed taxon order; using observed order as source of truth.")
    }
  }

  # Normalize null matrix names for downstream extraction and debugging.
  colnames(prev_perm_stats) <- paste0("perm_", seq_len(n_perm))
  colnames(abund_perm_stats) <- paste0("perm_", seq_len(n_perm))
  rownames(prev_perm_stats) <- observed_component_table$taxon_name
  rownames(abund_perm_stats) <- observed_component_table$taxon_name

  out <- list(
    observed_component_table = observed_component_table,
    observed_prev_stat = observed_prev_stat,
    observed_abund_stat = observed_abund_stat,
    prev_perm_stats = prev_perm_stats,
    abund_perm_stats = abund_perm_stats,
    perm_design = perm_design,
    perm_list = perm_list,
    tested_term = tested_term,
    tested_vars = tested_vars,
    n_perm = n_perm,
    prev_cache = prev_cache,
    abund_cache = abund_cache,
    prev_null = prev_null,
    abund_null = abund_null
  )

  class(out) <- c("pursue_component_permutations", "list")

  if (verbose) {
    message("Component resampling complete.")
    message("  Taxa: ", n_taxa)
    message("  Permutations: ", n_perm)
    message("  Prev usable: ", sum(observed_component_table$prev_usable_for_testing, na.rm = TRUE))
    message("  Abund usable: ", sum(observed_component_table$abund_usable_for_testing, na.rm = TRUE))
    message("  Permutation unit level: ", perm_design$unit_level)
  }

  out
}


# -----------------------------------------------------------------------------
# Component p-value and union-combination helpers
# -----------------------------------------------------------------------------

# Turn observed component statistics and null matrices into empirical p-values
compute_empirical_component_pvalues <- function(
    perm_obj,
    plus_one = TRUE,
    force_unusable_to_one = TRUE
) {
  if (!inherits(perm_obj, "pursue_component_permutations")) {
    stop("perm_obj must come from run_component_permutations().")
  }

  observed_table <- perm_obj$observed_component_table
  n_perm <- perm_obj$n_perm
  plus_one_addition <- if (plus_one) 1 else 0
  denominator <- n_perm + plus_one_addition

  prev_p_values <- rep(NA_real_, nrow(observed_table))
  abund_p_values <- rep(NA_real_, nrow(observed_table))

  # Compare each observed absolute statistic against its taxon-specific null row.
  for (taxon_pos in seq_len(nrow(observed_table))) {
    observed_prev <- perm_obj$observed_prev_stat[taxon_pos]
    observed_abund <- perm_obj$observed_abund_stat[taxon_pos]

    prev_null <- perm_obj$prev_perm_stats[taxon_pos, ]
    abund_null <- perm_obj$abund_perm_stats[taxon_pos, ]

    prev_p_values[taxon_pos] <- (
      plus_one_addition + sum(prev_null >= observed_prev, na.rm = TRUE)
    ) / denominator

    abund_p_values[taxon_pos] <- (
      plus_one_addition + sum(abund_null >= observed_abund, na.rm = TRUE)
    ) / denominator
  }

  # Neutralize components that failed support checks before union combination.
  if (force_unusable_to_one) {
    if ("prev_usable_for_testing" %in% names(observed_table)) {
      prev_p_values[!observed_table$prev_usable_for_testing] <- 1
    }

    if ("abund_usable_for_testing" %in% names(observed_table)) {
      abund_p_values[!observed_table$abund_usable_for_testing] <- 1
    }
  }

  out <- observed_table
  out$prev_empirical_p <- prev_p_values
  out$abund_empirical_p <- abund_p_values

  out
}


# Combine valid p-values with the Cauchy combination test
acat_single <- function(p, weights = NULL, clip = 1e-15) {
  p_values <- as.numeric(p)

  valid_values <- is.finite(p_values)
  if (!any(valid_values)) {
    return(NA_real_)
  }

  p_values <- p_values[valid_values]

  if (is.null(weights)) {
    weights_use <- rep(1 / length(p_values), length(p_values))
  } else {
    weights_use <- as.numeric(weights)[valid_values]
    weights_use <- weights_use / sum(weights_use)
  }

  # Avoid infinite tangent values at exact 0 or 1.
  p_values <- pmin(pmax(p_values, clip), 1 - clip)

  cauchy_stat <- sum(weights_use * tan((0.5 - p_values) * pi))
  acat_p <- 0.5 - atan(cauchy_stat) / pi
  acat_p <- min(max(acat_p, 0), 1)

  acat_p
}


# Combine the two component p-values into one taxon-level union p-value
combine_observed_components <- function(
    component_p_table,
    acat_weights = c(0.5, 0.5)
) {
  needed <- c(
    "taxon_index", "taxon_name",
    "prev_empirical_p", "abund_empirical_p",
    "prev_usable_for_testing", "abund_usable_for_testing"
  )

  missing_cols <- setdiff(needed, names(component_p_table))
  if (length(missing_cols) > 0L) {
    stop(
      "component_p_table is missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (length(acat_weights) != 2L) {
    stop("acat_weights must have length 2.")
  }

  n_taxa <- nrow(component_p_table)
  union_p <- rep(NA_real_, n_taxa)
  union_stat_label <- rep("acat", n_taxa)

  for (taxon_pos in seq_len(n_taxa)) {
    component_p_values <- c(
      component_p_table$prev_empirical_p[taxon_pos],
      component_p_table$abund_empirical_p[taxon_pos]
    )

    active_components <- c(
      component_p_table$prev_usable_for_testing[taxon_pos],
      component_p_table$abund_usable_for_testing[taxon_pos]
    )

    component_weights <- acat_weights

    # Drop unusable components and renormalize the remaining ACAT weights.
    if (sum(active_components) > 0L) {
      component_weights <- ifelse(active_components, component_weights, 0)

      if (sum(component_weights) <= 0) {
        component_weights <- as.numeric(active_components)
      }

      component_weights <- component_weights / sum(component_weights)
    } else {
      union_p[taxon_pos] <- 1
      next
    }

    union_p[taxon_pos] <- acat_single(
      component_p_values,
      weights = component_weights
    )
  }

  out <- component_p_table
  out$union_method <- "acat"
  out$union_stat_label <- union_stat_label
  out$union_p_value <- union_p

  out
}


# -----------------------------------------------------------------------------
# Independent filtering and final union layer
# -----------------------------------------------------------------------------

# Extract the statistic used for independent filtering
extract_independent_filter_stat <- function(union_table) {
  filter_stat <- if ("abund_n_obs" %in% names(union_table)) {
    as.numeric(union_table$abund_n_obs)
  } else {
    rep(NA_real_, nrow(union_table))
  }

  filter_stat[!is.finite(filter_stat)] <- NA_real_
  filter_stat
}


# Pick a filter cutoff and run BH on the surviving taxa
adaptive_independent_filter_bh <- function(
    p_value,
    filter_stat,
    alpha = 0.05,
    theta_grid = seq(0, 0.8, by = 0.02),
    method = "BH",
    filtered_q_value = NA_real_
) {
  if (length(p_value) != length(filter_stat)) {
    stop("p_value and filter_stat must have the same length.")
  }

  if (
    !is.numeric(alpha) ||
    length(alpha) != 1L ||
    !is.finite(alpha) ||
    alpha <= 0 ||
    alpha >= 1
  ) {
    stop("alpha must be a single finite number in (0, 1).")
  }

  # Restrict cutoff learning to taxa with usable p-values and filter statistics.
  valid_rows <- is.finite(p_value) & is.finite(filter_stat)
  p_values <- as.numeric(p_value[valid_rows])
  filter_values <- as.numeric(filter_stat[valid_rows])

  if (length(p_values) == 0L) {
    return(list(
      q_value = rep(filtered_q_value, length(p_value)),
      keep = rep(FALSE, length(p_value)),
      threshold = NA_real_,
      scan_table = data.frame(),
      filter_stat = filter_stat
    ))
  }

  # Convert the theta grid into concrete filter-statistic cutoffs.
  theta_grid <- unique(as.numeric(theta_grid[is.finite(theta_grid)]))
  theta_grid <- theta_grid[theta_grid >= 0 & theta_grid <= 1]

  if (length(theta_grid) == 0L) {
    theta_grid <- 0
  }

  candidate_cutoffs <- unique(as.numeric(stats::quantile(
    filter_values,
    probs = theta_grid,
    na.rm = TRUE,
    names = FALSE,
    type = 7
  )))

  candidate_cutoffs <- candidate_cutoffs[is.finite(candidate_cutoffs)]

  if (length(candidate_cutoffs) == 0L) {
    candidate_cutoffs <- min(filter_values, na.rm = TRUE)
  }

  candidate_cutoffs <- sort(unique(candidate_cutoffs))

  # Scan each cutoff by applying BH only to taxa that survive the filter.
  scan_list <- lapply(candidate_cutoffs, function(cutoff) {
    keep_rows <- filter_values >= cutoff
    q_values <- rep(filtered_q_value, length(p_values))

    if (sum(keep_rows) > 0L) {
      q_values[keep_rows] <- stats::p.adjust(
        p_values[keep_rows],
        method = method
      )
    }

    data.frame(
      threshold = cutoff,
      n_keep = sum(keep_rows),
      n_rej = sum(q_values <= alpha, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  scan_table <- do.call(rbind, scan_list)

  # Prefer the smallest cutoff among thresholds with the maximum rejection count.
  best_rows <- which(scan_table$n_rej == max(scan_table$n_rej, na.rm = TRUE))
  best_rows <- best_rows[which.min(scan_table$threshold[best_rows])]
  best_cutoff <- scan_table$threshold[best_rows]

  # Recompute the final q-values at the selected cutoff.
  keep_rows <- filter_values >= best_cutoff
  best_q_values <- rep(filtered_q_value, length(p_values))

  if (sum(keep_rows) > 0L) {
    best_q_values[keep_rows] <- stats::p.adjust(
      p_values[keep_rows],
      method = method
    )
  }

  # Scatter filtered q-values back into the original taxon order.
  q_out <- rep(filtered_q_value, length(p_value))
  keep_out <- rep(FALSE, length(p_value))

  q_out[valid_rows] <- best_q_values
  keep_out[valid_rows] <- keep_rows

  list(
    q_value = q_out,
    keep = keep_out,
    threshold = best_cutoff,
    scan_table = scan_table,
    filter_stat = filter_stat
  )
}


# Build the final union-layer statistics and adjusted p-values
compute_union_layer <- function(
    perm_obj = NULL,
    component_p_table = NULL,
    acat_weights = c(0.5, 0.5),
    independent_filtering = FALSE,
    independent_filter_alpha = 0.05,
    independent_filter_theta_grid = seq(0, 0.8, by = 0.02),
    independent_filter_filtered_q_value = NA_real_
) {
  # Allow callers to start from either a permutation object or a precomputed component table.
  if (is.null(component_p_table) && !is.null(perm_obj)) {
    if (!inherits(perm_obj, "pursue_component_permutations")) {
      stop("perm_obj must come from run_component_permutations().")
    }

    component_p_table <- compute_empirical_component_pvalues(perm_obj)
  }

  if (is.null(component_p_table)) {
    stop("For analytic union methods, supply component_p_table or perm_obj.")
  }

  # Combine prevalence and abundance evidence before any multiple-testing correction.
  union_table <- combine_observed_components(
    component_p_table = component_p_table,
    acat_weights = acat_weights
  )

  union_table$union_combined_stat_obs <- union_table$union_p_value

  union_q_raw <- stats::p.adjust(union_table$union_p_value, method = "BH")

  union_table$union_q_value_raw_bh <- union_q_raw
  union_table$if_enabled <- isTRUE(independent_filtering)
  union_table$if_filter_stat_name <- "abund_n_obs"
  union_table$if_filter_stat <- NA_real_
  union_table$if_keep <- NA
  union_table$if_threshold <- NA_real_

  union_q_final <- union_q_raw

  # Optionally learn an independent-filter threshold and use its BH-adjusted q-values.
  if (isTRUE(independent_filtering)) {
    filter_stat <- extract_independent_filter_stat(union_table = union_table)

    filter_result <- adaptive_independent_filter_bh(
      p_value = union_table$union_p_value,
      filter_stat = filter_stat,
      alpha = independent_filter_alpha,
      theta_grid = independent_filter_theta_grid,
      method = "BH",
      filtered_q_value = independent_filter_filtered_q_value
    )

    union_q_final <- filter_result$q_value
    union_table$if_filter_stat <- filter_result$filter_stat
    union_table$if_keep <- filter_result$keep
    union_table$if_threshold <- filter_result$threshold
  }

  # Store both the raw and final rejection calls for comparison during benchmarking.
  union_table$union_q_value <- union_q_final
  union_table$rejected_bh_0.05_raw <- ifelse(is.na(union_q_raw), FALSE, union_q_raw <= 0.05)
  union_table$rejected_bh_0.05 <- ifelse(is.na(union_q_final), FALSE, union_q_final <= 0.05)

  out <- list(
    result_table = union_table,
    method = "acat",
    independent_filtering = isTRUE(independent_filtering),
    independent_filter_stat = "abund_n_obs",
    independent_filter_alpha = independent_filter_alpha,
    independent_filter_theta_grid = independent_filter_theta_grid
  )

  out
}
