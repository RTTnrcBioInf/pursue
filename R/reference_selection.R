# -----------------------------------------------------------------------------
# Reference selector orchestrator
# -----------------------------------------------------------------------------

# Select a stable reference prefix using full-data scoring, bootstrap stability, and mass checks
select_reference_bootstrap_prefix <- function(
    X,
    metadata,
    formula,
    tested_term,
    p_max = 500L,
    pseudocount = 1,
    bootstrap_B = 100L,
    bootstrap_stability_threshold = 0.80,
    minimum_TA = NULL,
    minimum_rel_TA = NULL,
    hostage_prevention = TRUE,
    bulk_rel_quantile = 0.10,
    low_count_floor = 5L,
    tail_zero_tol = 1L,
    tail_low_prop_tol = 0.05,
    hostage_max_size_inflation = 0.50,
    hostage_max_score_inflation = 0.25,
    hostage_max_stability_drop = 0.05,
    return_ratio_matrix = FALSE,
    seed = NULL,
    verbose = FALSE
) {
  # Score the full data once; this ranking defines the canonical prefix path.
  score_obj <- score_reference_candidates(
    X = X,
    metadata = metadata,
    formula = formula,
    tested_term = tested_term,
    p_max = p_max,
    pseudocount = pseudocount,
    return_ratio_matrix = return_ratio_matrix
  )

  # Re-rank candidates on bootstrap-resampled samples to assess rank stability.
  bootstrap_rankings <- bootstrap_reference_rankings(
    X = X,
    metadata = metadata,
    formula = formula,
    tested_term = tested_term,
    p_max = p_max,
    pseudocount = pseudocount,
    n_bootstrap = bootstrap_B,
    seed = seed
  )

  # Convert bootstrap rankings into a stability score for every full-data prefix size.
  stability_table <- compute_prefix_stability(
    full_sorted_candidate_taxa = score_obj$sorted_candidate_taxa,
    bootstrap_rankings = bootstrap_rankings
  )

  # Summarize absolute and relative reference mass along the same prefix path.
  adequacy_table <- compute_prefix_adequacy(
    X = X,
    full_sorted_candidate_taxa = score_obj$sorted_candidate_taxa,
    rel_quantile = bulk_rel_quantile,
    low_count_floor = low_count_floor,
    use_relative = !is.null(minimum_rel_TA)
  )

  # Resolve the smallest prefix that satisfies stability and denominator-mass criteria.
  selection_obj <- find_bulk_feasible_prefix(
    prefix_stability = stability_table,
    prefix_adequacy = adequacy_table,
    stability_threshold = bootstrap_stability_threshold,
    minimum_TA = minimum_TA,
    minimum_rel_TA = minimum_rel_TA
  )

  selection_table <- selection_obj$selection_table

  # Optionally grow the prefix when the feasible cutoff still leaves fragile low-mass tails.
  if (hostage_prevention) {
    hostage_obj <- extend_fragile_reference_prefix(
      selection_table = selection_table,
      sorted_scores = score_obj$sorted_scores,
      n_samples = nrow(X),
      k0 = selection_obj$selected_k,
      low_count_floor = low_count_floor,
      tail_zero_tol = tail_zero_tol,
      tail_low_prop_tol = tail_low_prop_tol,
      hostage_max_size_inflation = hostage_max_size_inflation,
      hostage_max_score_inflation = hostage_max_score_inflation,
      hostage_max_stability_drop = hostage_max_stability_drop
    )

    selected_k <- hostage_obj$selected_k
  } else {
    hostage_obj <- list(
      selected_k = selection_obj$selected_k,
      hostage_base_k = selection_obj$selected_k,
      hostage_reason = "not_used"
    )

    selected_k <- selection_obj$selected_k
  }

  # Freeze final reference indices and names after any hostage-prevention expansion.
  selected_references <- score_obj$sorted_candidate_taxa[seq_len(selected_k)]
  selected_reference_names <- colnames(X)[selected_references]

  # Trace prefix-level prevalence and abundance diagnostics for downstream reporting.
  paths <- compute_reference_paths(
    X = X,
    sorted_candidate_taxa = score_obj$sorted_candidate_taxa
  )

  # Summarize the final selected denominator mass across samples.
  mass_diag <- summarize_reference_mass(
    X = X,
    ref_idx = selected_references,
    include_zero_count = FALSE
  )

  # Return the selected set together with scoring, stability, and mass diagnostics.
  out <- list(
    X = X,
    method = "bootstrap_prefix",
    selected_references = selected_references,
    selected_reference_names = selected_reference_names,
    selected_cutpoint = selected_k,
    selected_MinAbundance = unname(min(rowSums(X[, selected_references, drop = FALSE]))),

    scores = score_obj$scores_full,
    sorted_candidate_taxa = score_obj$sorted_candidate_taxa,
    sorted_candidate_names = score_obj$sorted_candidate_names,
    sorted_scores = score_obj$sorted_scores,

    ratio_matrix = score_obj$ratio_matrix,
    mean_prevalence_over_the_sorted = unname(paths$mean_prevalence_over_the_sorted),
    min_abundance_over_the_sorted = unname(paths$min_abundance_over_the_sorted),
    which_is_min_abundance_over_the_sorted = unname(paths$which_is_min_abundance_over_the_sorted),
    reference_mass = mass_diag$reference_mass,
    reference_mass_min = mass_diag$reference_mass_min,
    reference_mass_median = mass_diag$reference_mass_median,
    reference_mass_max = mass_diag$reference_mass_max,

    bootstrap_B = bootstrap_B,
    bootstrap_stability_threshold = bootstrap_stability_threshold,
    minimum_TA = minimum_TA,
    minimum_rel_TA = minimum_rel_TA,
    bootstrap_hostage_prevention = hostage_prevention,
    bootstrap_bulk_rel_quantile = bulk_rel_quantile,
    bootstrap_low_count_floor = low_count_floor,
    bootstrap_tail_zero_tol = tail_zero_tol,
    bootstrap_tail_low_prop_tol = tail_low_prop_tol,
    bootstrap_hostage_max_size_inflation = hostage_max_size_inflation,
    bootstrap_hostage_max_score_inflation = hostage_max_score_inflation,
    bootstrap_hostage_max_stability_drop = hostage_max_stability_drop,
    bootstrap_bulk_base_k = selection_obj$selected_k,
    bootstrap_hostage_base_k = hostage_obj$hostage_base_k,
    bootstrap_hostage_reason = hostage_obj$hostage_reason,
    bootstrap_selection_table = selection_table,

    n_samples = nrow(X),
    n_taxa = ncol(X),
    n_candidates = length(score_obj$candidate_idx),
    pseudocount_input = pseudocount
  )

  if (verbose) {
    message("Bootstrap-prefix reference selection complete.")
    message("  Selected reference size: ", length(selected_references))

    if (hostage_prevention) {
      message("  Hostage prevention base k: ", selection_obj$selected_k)
      message("  Hostage prevention reason: ", hostage_obj$hostage_reason)
    }
  }

  invisible(out)
}


# -----------------------------------------------------------------------------
# Reference-selection component functions
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Post-selection cleanup and prefix diagnostics
# -----------------------------------------------------------------------------

# Remove the most group-shifted references after primary reference selection
cleanup_reference_group_shift <- function(
    X,
    ref_idx,
    metadata,
    tested_term,
    frac = 0.10,
    min_ref = 10L,
    pseudocount = 1
) {
  # Normalize the selected references before applying any cleanup rule.
  ref_idx <- as.integer(unique(ref_idx))

  # Keep the selected set unchanged when cleanup would leave too few references.
  if (length(ref_idx) <= min_ref) {
    return(list(
      selected_references = ref_idx,
      removed_references = integer(0),
      group_cleanup_scores = setNames(rep(NA_real_, length(ref_idx)), colnames(X)[ref_idx])
    ))
  }

  # Score each reference against the remaining reference mass across tested groups.
  scores <- score_reference_group_shift(
    X = X,
    ref_idx = ref_idx,
    metadata = metadata,
    tested_term = tested_term,
    pseudocount = pseudocount
  )

  # Cap the number of dropped references so the retained set respects min_ref.
  n_drop <- min(
    floor(length(ref_idx) * frac),
    length(ref_idx) - min_ref
  )

  # Return the original set if the requested fraction rounds to no removable taxa.
  if (n_drop <= 0L) {
    return(list(
      selected_references = ref_idx,
      removed_references = integer(0),
      group_cleanup_scores = scores
    ))
  }

  # Drop the references with the largest apparent tested-group shift.
  score_order <- order(scores, decreasing = TRUE, na.last = TRUE)
  drop_idx <- ref_idx[score_order[seq_len(n_drop)]]
  keep_idx <- setdiff(ref_idx, drop_idx)

  list(
    selected_references = keep_idx,
    removed_references = drop_idx,
    group_cleanup_scores = scores
  )
}

# Trace prevalence and abundance diagnostics across the sorted reference prefix
compute_reference_paths <- function(X, sorted_candidate_taxa) {
  candidate_counts <- X[, sorted_candidate_taxa, drop = FALSE]

  # Track the weakest sample-level denominator as references are added.
  cumulative_abundance <- t(apply(candidate_counts, 1L, cumsum))
  min_abundance_over_the_sorted <- apply(cumulative_abundance, 2L, min)
  which_is_min_abundance_over_the_sorted <- apply(cumulative_abundance, 2L, which.min)

  # Track how many samples have at least one selected reference observed.
  presence <- candidate_counts > 0
  cumulative_presence <- t(apply(presence, 1L, cummax))
  mean_prevalence_over_the_sorted <- colMeans(cumulative_presence)

  list(
    mean_prevalence_over_the_sorted = mean_prevalence_over_the_sorted,
    min_abundance_over_the_sorted = min_abundance_over_the_sorted,
    which_is_min_abundance_over_the_sorted = which_is_min_abundance_over_the_sorted
  )
}


# -----------------------------------------------------------------------------
# Candidate scoring and bootstrap ranking
# -----------------------------------------------------------------------------

# Cap the candidate pool before quadratic pairwise reference scoring
select_reference_candidates <- function(X, p_max = 500L) {
  candidate_idx <- seq_len(ncol(X))

  # Keep the highest-total-count taxa when all taxa cannot be scored pairwise.
  if (length(candidate_idx) > p_max) {
    abundance_order <- order(
      colSums(X[, candidate_idx, drop = FALSE]),
      decreasing = TRUE
    )

    candidate_idx <- candidate_idx[abundance_order[seq_len(p_max)]]
  }

  candidate_idx
}

# Rank candidate references by nuisance-adjusted pairwise log-ratio stability
score_reference_candidates <- function(
    X,
    metadata,
    formula,
    tested_term,
    p_max = 500L,
    pseudocount = 1,
    return_ratio_matrix = TRUE
) {
  # Standardize count and metadata inputs before formula-derived operations.
  X <- as_numeric_matrix(X, "X")
  metadata <- as_plain_data_frame(metadata)

  if (nrow(X) != nrow(metadata)) {
    stop("nrow(X) must equal nrow(metadata).")
  }

  # Add stable taxon and sample labels when the caller did not provide them.
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("Taxon_", seq_len(ncol(X)))
  }

  if (is.null(rownames(X))) {
    rownames(X) <- paste0("Sample_", seq_len(nrow(X)))
  }

  # Build the nuisance-only design used to remove covariate-driven ratio variation.
  nuisance_matrix <- make_nuisance_matrix(
    metadata = metadata,
    formula = formula,
    tested_term = tested_term
  )

  # Restrict scoring to the highest-abundance candidates when the feature space is large.
  candidate_idx <- select_reference_candidates(X, p_max)
  candidate_counts <- X[, candidate_idx, drop = FALSE]

  # Work on taxa x samples log counts so each row is a candidate reference.
  candidate_table <- t(candidate_counts) + pseudocount
  log_candidate_counts <- log(candidate_table)

  n_candidates <- nrow(log_candidate_counts)
  n_samples <- ncol(log_candidate_counts)

  # Precompute the nuisance projection reused for every pairwise log-ratio model.
  xtx_inverse <- tryCatch(
    solve(crossprod(nuisance_matrix)),
    error = function(e) qr.solve(crossprod(nuisance_matrix))
  )

  projection_solver <- xtx_inverse %*% t(nuisance_matrix)
  residual_df <- n_samples - ncol(nuisance_matrix)

  # Allocate the symmetric pairwise instability matrix.
  pairwise_sd <- matrix(
    NA_real_,
    nrow = n_candidates,
    ncol = n_candidates,
    dimnames = list(colnames(candidate_counts), colnames(candidate_counts))
  )

  # Estimate nuisance-adjusted residual SD for every candidate pair.
  for (candidate_pos in seq_len(n_candidates - 1L)) {
    # Compare one candidate against all candidates ranked after it to fill the upper triangle.
    comparison_block <- log_candidate_counts[(candidate_pos + 1L):n_candidates, , drop = FALSE]

    # Build candidate-vs-comparator log ratios in sample x pair orientation.
    log_ratio_block <- matrix(
      log_candidate_counts[candidate_pos, ],
      nrow = nrow(comparison_block),
      ncol = ncol(comparison_block),
      byrow = TRUE
    ) - comparison_block

    log_ratio_block <- t(log_ratio_block)

    # Residualize each log-ratio vector against nuisance covariates.
    coef_estimates <- projection_solver %*% log_ratio_block
    residual_matrix <- log_ratio_block - nuisance_matrix %*% coef_estimates
    residual_sd <- sqrt(colSums(residual_matrix^2) / residual_df)

    # Mirror the block because pairwise instability is symmetric.
    pairwise_sd[candidate_pos, (candidate_pos + 1L):n_candidates] <- residual_sd
    pairwise_sd[(candidate_pos + 1L):n_candidates, candidate_pos] <- residual_sd
  }

  diag(pairwise_sd) <- NA_real_

  # Rank candidates by their median instability against all other candidates.
  median_pairwise_score <- apply(pairwise_sd, 1L, median, na.rm = TRUE)
  score_order <- order(median_pairwise_score, na.last = TRUE)

  sorted_candidate_taxa <- candidate_idx[score_order]
  sorted_candidate_names <- colnames(X)[sorted_candidate_taxa]
  sorted_scores <- median_pairwise_score[score_order]

  # Expand candidate-only scores back to the full taxon space for reporting.
  scores_full <- rep(Inf, ncol(X))
  names(scores_full) <- colnames(X)
  scores_full[candidate_idx] <- median_pairwise_score

  list(
    X = X,
    candidate_idx = candidate_idx,
    ratio_matrix = if (return_ratio_matrix) pairwise_sd else NULL,
    median_pairwise_score = median_pairwise_score,
    sorted_candidate_taxa = sorted_candidate_taxa,
    sorted_candidate_names = sorted_candidate_names,
    sorted_scores = sorted_scores,
    scores_full = scores_full,
    pseudocount_input = pseudocount
  )
}

# Re-rank reference candidates after bootstrap resampling of samples
bootstrap_reference_rankings <- function(
    X,
    metadata,
    formula,
    tested_term,
    p_max = 500L,
    pseudocount = 1,
    n_bootstrap = 100L,
    seed = NULL
) {
  # Use an optional seed so bootstrap rankings can be reproduced exactly.
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n_samples <- nrow(X)
  bootstrap_rankings <- vector("list", n_bootstrap)

  # Each bootstrap replicate gets an independent candidate ranking.
  for (boot_pos in seq_len(n_bootstrap)) {
    sample_rows <- sample.int(n_samples, size = n_samples, replace = TRUE)

    # Resample counts and metadata together to preserve sample-level covariates.
    boot_counts <- X[sample_rows, , drop = FALSE]
    boot_metadata <- metadata[sample_rows, , drop = FALSE]

    # Force unique row names after sampling with replacement.
    rownames(boot_counts) <- paste0("boot_", boot_pos, "_", seq_len(nrow(boot_counts)))
    rownames(boot_metadata) <- rownames(boot_counts)

    # Score the bootstrap replicate without storing its full pairwise matrix.
    boot_score <- score_reference_candidates(
      X = boot_counts,
      metadata = boot_metadata,
      formula = formula,
      tested_term = tested_term,
      p_max = p_max,
      pseudocount = pseudocount,
      return_ratio_matrix = FALSE
    )

    bootstrap_rankings[[boot_pos]] <- list(
      sorted_candidate_taxa = boot_score$sorted_candidate_taxa,
      sorted_scores = boot_score$sorted_scores
    )
  }

  bootstrap_rankings
}


# -----------------------------------------------------------------------------
# Prefix feasibility and hostage prevention
# -----------------------------------------------------------------------------

# Estimate how often each full-data ranked prefix survives bootstrap resampling
compute_prefix_stability <- function(
    full_sorted_candidate_taxa,
    bootstrap_rankings,
    k_grid = NULL
) {
  n_candidates <- length(full_sorted_candidate_taxa)

  # Evaluate every prefix size unless the caller provides a reduced grid.
  if (is.null(k_grid)) {
    k_grid <- seq_len(n_candidates)
  }

  n_bootstrap <- length(bootstrap_rankings)
  stability <- numeric(length(k_grid))

  # For each prefix size, average per-reference bootstrap inclusion rates.
  for (grid_pos in seq_along(k_grid)) {
    prefix_size <- k_grid[grid_pos]
    full_prefix <- full_sorted_candidate_taxa[seq_len(prefix_size)]
    inclusion_rate <- numeric(prefix_size)

    # Estimate how often each full-data prefix member stays inside the same-size bootstrap prefix.
    for (prefix_pos in seq_along(full_prefix)) {
      taxon_id <- full_prefix[prefix_pos]

      inclusion_rate[prefix_pos] <- mean(vapply(
        bootstrap_rankings,
        function(boot_entry) {
          boot_prefix_size <- min(prefix_size, length(boot_entry$sorted_candidate_taxa))
          taxon_id %in% boot_entry$sorted_candidate_taxa[seq_len(boot_prefix_size)]
        },
        logical(1L)
      ))
    }

    # Collapse per-reference inclusion into one stability score for this prefix size.
    stability[grid_pos] <- mean(inclusion_rate)
  }

  data.frame(
    k = k_grid,
    stability = stability,
    n_boot = n_bootstrap,
    stringsAsFactors = FALSE
  )
}

# Measure how much reference mass each ranked prefix carries
compute_prefix_adequacy <- function(
    X,
    full_sorted_candidate_taxa,
    k_grid = NULL,
    rel_quantile = 0.10,
    low_count_floor = 5L,
    use_relative = TRUE
) {
  n_candidates <- length(full_sorted_candidate_taxa)

  if (is.null(k_grid)) {
    k_grid <- seq_len(n_candidates)
  }

  sample_depth <- rowSums(X)
  adequacy_rows <- vector("list", length(k_grid))

  # Summarize absolute and optional relative mass for every candidate prefix size.
  for (grid_pos in seq_along(k_grid)) {
    prefix_size <- k_grid[grid_pos]
    reference_indices <- full_sorted_candidate_taxa[seq_len(prefix_size)]

    # Compute the denominator supplied by this prefix in every sample.
    reference_mass <- rowSums(X[, reference_indices, drop = FALSE])
    relative_reference_mass <- if (use_relative) {
      reference_mass / pmax(sample_depth, 1)
    } else {
      rep(NA_real_, length(reference_mass))
    }

    # Store both lower-tail mass and explicit zero/low-count failure counts.
    adequacy_rows[[grid_pos]] <- data.frame(
      k = prefix_size,
      ref_mass_min = min(reference_mass),
      ref_mass_median = unname(stats::median(reference_mass)),
      ref_mass_q05 = as.numeric(stats::quantile(reference_mass, 0.05, names = FALSE)),
      rel_mass_q = if (use_relative) {
        as.numeric(stats::quantile(relative_reference_mass, rel_quantile, names = FALSE))
      } else {
        NA_real_
      },
      n_zero_ref_mass = sum(reference_mass == 0),
      n_low_ref_mass = sum(reference_mass < low_count_floor),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, adequacy_rows)
}

# Pick the smallest prefix that clears the stability and adequacy thresholds
find_bulk_feasible_prefix <- function(
    prefix_stability,
    prefix_adequacy,
    stability_threshold = 0.80,
    minimum_TA = NULL,
    minimum_rel_TA = NULL
) {
  # Align stability and mass summaries by prefix size.
  selection_table <- merge(
    prefix_stability,
    prefix_adequacy,
    by = "k",
    sort = TRUE
  )

  # Apply stability first, then optional absolute and relative mass constraints.
  feasible <- rep(TRUE, nrow(selection_table))
  feasible <- feasible &
    is.finite(selection_table$stability) &
    (selection_table$stability >= stability_threshold)

  if (!is.null(minimum_TA)) {
    feasible <- feasible & (selection_table$ref_mass_min >= minimum_TA)
  }

  if (!is.null(minimum_rel_TA)) {
    feasible <- feasible &
      is.finite(selection_table$rel_mass_q) &
      (selection_table$rel_mass_q >= minimum_rel_TA)
  }

  feasible_rows <- which(feasible)

  # Fail loudly when every candidate prefix is either unstable or too low-mass.
  if (length(feasible_rows) == 0L) {
    stop("No prefix satisfied the bootstrap stability and bulk adequacy criteria.")
  }

  # Choose the earliest feasible prefix to keep the denominator compact.
  selected_row <- feasible_rows[1L]

  list(
    selected_k = selection_table$k[selected_row],
    selection_table = selection_table,
    selected_row = selection_table[selected_row, , drop = FALSE]
  )
}

# Extend the selected prefix when the initial cutoff leaves fragile low-mass tails
extend_fragile_reference_prefix <- function(
    selection_table,
    sorted_scores,
    n_samples,
    k0,
    low_count_floor = 5L,
    tail_zero_tol = 1L,
    tail_low_prop_tol = 0.05,
    hostage_max_size_inflation = 0.50,
    hostage_max_score_inflation = 0.25,
    hostage_max_stability_drop = 0.05
) {
  # Use the bulk-feasible prefix as the anchor for bounded expansion.
  base_row <- selection_table[selection_table$k == k0, , drop = FALSE]

  if (nrow(base_row) != 1L) {
    stop("k0 not found in selection_table.")
  }

  # Anchor score and stability define the collateral-damage budget.
  base_stability <- base_row$stability
  base_score <- sorted_scores[k0]

  if (!is.finite(base_score) || base_score <= 0) {
    base_score <- 1e-8
  }

  # Convert the allowed low-mass tail fraction into a sample count.
  max_low_tail_samples <- max(1L, ceiling(tail_low_prop_tol * n_samples))
  selected_k <- k0
  selection_reason <- "bulk_feasible"

  # Keep the bulk-feasible cutoff if it already has acceptable low-mass tails.
  tail_ok_at_start <- (base_row$n_zero_ref_mass <= tail_zero_tol) &&
    (base_row$n_low_ref_mass <= max_low_tail_samples)

  if (tail_ok_at_start) {
    return(list(
      selected_k = selected_k,
      hostage_base_k = k0,
      hostage_reason = "tail_ok_at_bulk_cutoff"
    ))
  }

  later_prefix_sizes <- selection_table$k[selection_table$k > k0]

  # Walk forward until tails are rescued or the added references become too costly.
  for (prefix_size in later_prefix_sizes) {
    candidate_row <- selection_table[selection_table$k == prefix_size, , drop = FALSE]

    # Check whether this larger prefix fixes zero or low reference-mass samples.
    tail_ok <- (candidate_row$n_zero_ref_mass <= tail_zero_tol) &&
      (candidate_row$n_low_ref_mass <= max_low_tail_samples)

    # Bound expansion by size growth, added candidate instability, and stability loss.
    size_inflation <- (prefix_size - k0) / max(k0, 1)
    score_inflation <- (sorted_scores[prefix_size] - base_score) / base_score
    stability_drop <- base_stability - candidate_row$stability

    collateral_ok <- (size_inflation <= hostage_max_size_inflation) &&
      (score_inflation <= hostage_max_score_inflation) &&
      (stability_drop <= hostage_max_stability_drop)

    if (!collateral_ok) {
      selection_reason <- "hostage_stop_collateral_damage"
      break
    }

    # Accept the larger prefix provisionally while it stays within collateral limits.
    selected_k <- prefix_size

    if (tail_ok) {
      selection_reason <- "tail_rescued"
      break
    }
  }

  list(
    selected_k = selected_k,
    hostage_base_k = k0,
    hostage_reason = selection_reason
  )
}


# -----------------------------------------------------------------------------
# Group-aware reference cleanup
# -----------------------------------------------------------------------------

# Convert the tested term into a binary grouping factor for cleanup scoring
make_group_cleanup_factor <- function(metadata, tested_term) {
  # Read the tested column as the grouping variable used only for cleanup scoring.
  metadata <- as_plain_data_frame(metadata)
  group_factor <- metadata[[tested_term]]

  if (is.null(group_factor)) {
    stop("tested_term not found in metadata for group-aware reference cleanup.")
  }

  # Cleanup currently uses a two-group contrast rather than a general multi-df term.
  group_factor <- droplevels(as.factor(group_factor))

  if (nlevels(group_factor) != 2L) {
    stop("group-aware reference cleanup expects a binary tested_term.")
  }

  group_factor
}

# Estimate how much each reference shifts across tested-term groups
score_reference_group_shift <- function(
    X,
    ref_idx,
    metadata,
    tested_term,
    pseudocount = 1
) {
  # Convert the tested term to the binary grouping used for mean-shift scoring.
  group_factor <- make_group_cleanup_factor(metadata, tested_term)
  group_levels <- levels(group_factor)

  group1_rows <- which(group_factor == group_levels[1L])
  group2_rows <- which(group_factor == group_levels[2L])

  scores <- numeric(length(ref_idx))
  names(scores) <- colnames(X)[ref_idx]

  # Score each reference as its group difference against the remaining reference mass.
  for (ref_pos in seq_along(ref_idx)) {
    reference_index <- ref_idx[ref_pos]
    other_references <- setdiff(ref_idx, reference_index)

    if (length(other_references) < 1L) {
      scores[ref_pos] <- 0
      next
    }

    # Compare each reference against the pooled mass of the other selected references.
    other_reference_mass <- rowSums(X[, other_references, drop = FALSE])
    log_reference_ratio <- log(
      (X[, reference_index] + pseudocount) /
        (other_reference_mass + pseudocount)
    )

    # Use the absolute two-group mean difference as the cleanup removal score.
    scores[ref_pos] <- abs(
      mean(log_reference_ratio[group1_rows], na.rm = TRUE) -
        mean(log_reference_ratio[group2_rows], na.rm = TRUE)
    )
  }

  scores
}
