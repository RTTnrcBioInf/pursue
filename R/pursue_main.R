#' Run the PURSUE workflow
#'
#' `run_pursue()` fits the full PURSUE workflow for microbial differential
#' abundance analysis. The procedure filters low-information features,
#' selects a stable reference set, constructs reference-normalized prevalence and
#' abundance responses, fits the two model components separately, computes
#' empirical p-values using a Freedman--Lane residual permutation scheme,
#' optionally applies residual winsorization within the abundance permutation
#' engine, and combines component-level evidence into union-level results with
#' optional adaptive independent filtering.
#'
#' @param otu A feature count matrix. Rows correspond to samples and columns
#'   correspond to microbial features.
#' @param meta A sample metadata data frame. Rows must correspond to the samples
#'   in `otu`.
#' @param formula A model formula specifying the tested term and any nuisance
#'   covariates.
#' @param tested_term A single term from `formula` to be tested.
#' @param tested_vars Optional character vector specifying the variable or
#'   variables to be permuted in the Freedman--Lane scheme. If `NULL`, defaults
#'   to `tested_term`.
#' @param taxa_keep Optional full length logical or integer index specifying features to
#'   retain before model fitting. Default is `NULL`.
#' @param min_prevalence_filter Optional prevalence threshold for filtering
#'   low-prevalence features. Default is `NULL`.
#' @param min_total_count_filter Optional total-count threshold for filtering
#'   low-count features. Default is `NULL`.
#' @param reference_pseudocount Pseudocount used in log-ratio reference scoring.
#' @param ref_n_max Maximum number of candidate reference features considered
#'   during reference scoring.
#' @param bootstrap_B Number of bootstrap resamples used for bootstrap-prefix
#'   reference selection.
#' @param bootstrap_stability_threshold Minimum prefix stability required for
#'   bootstrap-prefix reference selection.
#' @param minimum_TA Minimum total reference abundance threshold used in
#'   bootstrap-prefix reference selection. Default is `10`.
#' @param minimum_rel_TA Minimum relative reference abundance threshold used in
#'   bootstrap-prefix reference selection. Default is `0.01`.
#' @param hostage_prevention Logical; if `TRUE`, extends the selected reference
#'   prefix when needed to reduce fragile low-mass tails.
#' @param bulk_rel_quantile Quantile used when evaluating relative reference
#'   adequacy.
#' @param low_count_floor Count threshold used in bootstrap adequacy and hostage
#'   prevention diagnostics.
#' @param tail_zero_tol Maximum tolerated number of zero reference-mass samples
#'   in hostage prevention.
#' @param tail_low_prop_tol Maximum tolerated proportion of low reference-mass
#'   samples in hostage prevention.
#' @param hostage_max_size_inflation Maximum allowed reference-set size inflation
#'   during hostage prevention.
#' @param hostage_max_score_inflation Maximum allowed score inflation during
#'   hostage prevention.
#' @param hostage_max_stability_drop Maximum allowed stability loss during
#'   hostage prevention.
#' @param reference_group_cleanup Logical; if `TRUE`, performs post-selection
#'   group-shift cleanup of the reference set.
#' @param reference_group_cleanup_frac Fraction of selected references removed
#'   during group-shift cleanup.
#' @param reference_group_cleanup_min_ref Minimum number of references retained
#'   after group-shift cleanup.
#' @param reference_group_cleanup_pseudocount Pseudocount used during
#'   group-shift cleanup.
#' @param depth_adjust Logical; if `TRUE`, appends `log(depth)` to the model
#'   formula. Default is `TRUE`.
#' @param prevalence_expected_rarefied_eps Numerical clipping constant used in
#'   expected-rarefied prevalence construction.
#' @param parallelize_permutations Logical; if `TRUE`, permutations are
#'   parallelized.
#' @param permutation_n_cores Number of cores used when
#'   `parallelize_permutations = TRUE`.
#' @param abundance_min_nonzero_n Minimum number of observed nonzero abundance
#'   values required to fit an abundance model for a feature.
#' @param abund_resid_winsorize Logical; if `TRUE`, winsorizes abundance-model
#'   residuals in the permutation engine.
#' @param abund_resid_winsor_quantile Quantile used for residual winsorization.
#' @param abund_resid_winsor_min_n Minimum sample size required before residual
#'   winsorization is applied.
#' @param n_perm Number of permutations used for empirical p-value estimation.
#' @param cluster_id Optional clustering variable for clustered permutation.
#'   Useful for repeated-measures or paired designs. Default is `NULL`.
#' @param strata Optional stratification variable restricting permutations to
#'   occur within strata. Default is `NULL`.
#' @param cluster_exposure_mode Controls whether permutation is performed at the
#'   sample or cluster level. One of `"auto"`, `"sample"`, or `"cluster"`.
#' @param seed Optional random seed.
#' @param union_acat_weights Numeric vector of length two giving prevalence and
#'   abundance weights for ACAT union p-value combination.
#' @param independent_filtering Logical; if `TRUE`, applies adaptive independent
#'   filtering before multiple-testing correction.
#' @param independent_filter_alpha Target FDR level used by adaptive independent
#'   filtering.
#' @param independent_filter_theta_grid Candidate filtering thresholds used by
#'   the adaptive independent filtering procedure.
#' @param independent_filter_filtered_q_value Q-value assigned to filtered-out
#'   features.
#' @param keep_diagnostics Logical; if `TRUE`, returns intermediate diagnostics
#'   in addition to the main results table.
#' @param verbose Logical; if `TRUE`, prints progress messages.
#'
#' @details
#' PURSUE models microbial differential signals using two reference-normalized
#' components: a prevalence component and an abundance component. Reference
#' features are selected before response construction, after which the two
#' components are fitted separately and evaluated using empirical p-values from
#' a Freedman--Lane residual permutation scheme. Component-level evidence is
#' then combined into union-level results, with optional adaptive independent
#' filtering before multiple-testing correction.
#'
#' @return
#' If `keep_diagnostics = FALSE`, a list with `results_table`. If
#' `keep_diagnostics = TRUE`, a list containing the final `results_table` plus
#' the filter object, reference-selection output, constructed responses,
#' observed component fits, permutation outputs, and union-layer output.
#'
#' @examples
#' set.seed(1)
#'
#' otu <- matrix(
#'   c(
#'     12, 0, 3, 5,
#'      8, 1, 2, 4,
#'     15, 0, 1, 6,
#'      3, 7, 0, 8,
#'      2, 9, 1, 10,
#'      1, 11, 0, 9
#'   ),
#'   nrow = 6,
#'   byrow = TRUE
#' )
#'
#' colnames(otu) <- paste0("Taxon_", seq_len(ncol(otu)))
#' rownames(otu) <- paste0("Sample_", seq_len(nrow(otu)))
#'
#' meta <- data.frame(
#'   group = factor(c("A", "A", "A", "B", "B", "B")),
#'   row.names = rownames(otu)
#' )
#'
#' out <- run_pursue(
#'   otu = otu,
#'   meta = meta,
#'   formula = ~ group,
#'   tested_term = "group",
#'   n_perm = 99,
#'   verbose = TRUE
#' )
#'
#' head(out$results_table)
#'
#' @export

# -----------------------------------------------------------------------------
# PURSUE workflow entry point
# -----------------------------------------------------------------------------

# Run filtering, reference selection, response construction, component testing, and union calls
run_pursue <- function(
    otu,
    meta,
    formula,
    tested_term,
    tested_vars = NULL,

    taxa_keep = NULL,
    min_prevalence_filter = NULL,
    min_total_count_filter = NULL,

    reference_pseudocount = 1,
    ref_n_max = 100000000L,
    bootstrap_B = 100L,
    bootstrap_stability_threshold = 0.80,
    minimum_TA = 10,
    minimum_rel_TA = 0.01,
    hostage_prevention = TRUE,
    bulk_rel_quantile = 0.10,
    low_count_floor = 5L,
    tail_zero_tol = 1L,
    tail_low_prop_tol = 0.05,
    hostage_max_size_inflation = 0.50,
    hostage_max_score_inflation = 0.25,
    hostage_max_stability_drop = 0.05,
    reference_group_cleanup = TRUE,
    reference_group_cleanup_frac = 0.10,
    reference_group_cleanup_min_ref = 10L,
    reference_group_cleanup_pseudocount = 1,

    depth_adjust = TRUE,

    prevalence_expected_rarefied_eps = 1e-12,
    parallelize_permutations = FALSE,
    permutation_n_cores = max(1L, parallel::detectCores() - 1L),

    abundance_min_nonzero_n = 3L,
    abund_resid_winsorize = TRUE,
    abund_resid_winsor_quantile = 0.05,
    abund_resid_winsor_min_n = 20L,

    n_perm = 1999,
    cluster_id = NULL,
    strata = NULL,
    cluster_exposure_mode = c("auto", "sample", "cluster"),
    seed = NULL,

    union_acat_weights = c(0.5, 0.5),

    independent_filtering = FALSE,
    independent_filter_alpha = 0.05,
    independent_filter_theta_grid = seq(0, 0.8, by = 0.02),
    independent_filter_filtered_q_value = NA_real_,

    keep_diagnostics = FALSE,
    verbose = TRUE
) {
  # Normalize input containers before any row or column bookkeeping.
  otu <- as_numeric_matrix(otu, "otu")
  meta <- as_plain_data_frame(meta)

  if (nrow(otu) != nrow(meta)) {
    stop("nrow(otu) must equal nrow(meta).")
  }

  # Fill missing identifiers so downstream diagnostics always have stable labels.
  if (is.null(rownames(otu))) {
    rownames(otu) <- paste0("Sample_", seq_len(nrow(otu)))
  }

  if (is.null(colnames(otu))) {
    colnames(otu) <- paste0("Taxon_", seq_len(ncol(otu)))
  }

  if (is.null(rownames(meta))) {
    rownames(meta) <- rownames(otu)
  }

  # Align metadata to the count matrix when the same samples are present out of order.
  if (!identical(rownames(otu), rownames(meta))) {
    if (setequal(rownames(otu), rownames(meta))) {
      meta <- meta[rownames(otu), , drop = FALSE]
    } else {
      stop("Row names of otu and meta do not match.")
    }
  }

  # Optionally append log library size as a nuisance covariate.
  if (depth_adjust) {
    meta <- add_log_depth(
      otu = otu,
      metadata = meta,
      depth_var_name = "log_depth"
    )

    formula <- add_term_if_missing(formula, "log_depth")
  }

  cluster_exposure_mode <- match.arg(cluster_exposure_mode)

  # Apply explicit, prevalence, and total-count filters before reference selection.
  filter_obj <- filter_taxa(
    otu = otu,
    taxa_keep = taxa_keep,
    min_prevalence = min_prevalence_filter,
    min_total_count = min_total_count_filter,
    verbose = verbose
  )

  filtered_otu <- filter_obj$otu

  if (verbose) message(strrep("-", 50))

  # Select references on the filtered feature set using nuisance-adjusted bootstrap stability.
  reference_fit <- select_reference_bootstrap_prefix(
    X = filtered_otu,
    metadata = meta,
    formula = formula,
    tested_term = tested_term,
    p_max = ref_n_max,
    pseudocount = reference_pseudocount,
    bootstrap_B = bootstrap_B,
    bootstrap_stability_threshold = bootstrap_stability_threshold,
    minimum_TA = minimum_TA,
    minimum_rel_TA = minimum_rel_TA,
    hostage_prevention = hostage_prevention,
    bulk_rel_quantile = bulk_rel_quantile,
    low_count_floor = low_count_floor,
    tail_zero_tol = tail_zero_tol,
    tail_low_prop_tol = tail_low_prop_tol,
    hostage_max_size_inflation = hostage_max_size_inflation,
    hostage_max_score_inflation = hostage_max_score_inflation,
    hostage_max_stability_drop = hostage_max_stability_drop,
    return_ratio_matrix = keep_diagnostics,
    seed = seed,
    verbose = verbose
  )

  final_reference_idx <- reference_fit$selected_references
  final_reference_names <- reference_fit$selected_reference_names

  # Optionally remove selected references that show the largest tested-term group shift.
  if (isTRUE(reference_group_cleanup)) {
    cleanup_obj <- cleanup_reference_group_shift(
      X = filtered_otu,
      ref_idx = final_reference_idx,
      metadata = meta,
      tested_term = tested_term,
      frac = reference_group_cleanup_frac,
      min_ref = reference_group_cleanup_min_ref,
      pseudocount = reference_group_cleanup_pseudocount
    )

    final_reference_idx <- cleanup_obj$selected_references
    final_reference_names <- colnames(filtered_otu)[final_reference_idx]

    # Replace the selected reference set while retaining cleanup diagnostics.
    reference_fit$selected_references <- final_reference_idx
    reference_fit$selected_reference_names <- final_reference_names
    reference_fit$group_cleanup_removed_references <- cleanup_obj$removed_references
    reference_fit$group_cleanup_removed_names <- colnames(filtered_otu)[cleanup_obj$removed_references]
    reference_fit$group_cleanup_scores <- cleanup_obj$group_cleanup_scores
    reference_fit$group_cleanup_frac <- reference_group_cleanup_frac
    reference_fit$group_cleanup_min_ref <- reference_group_cleanup_min_ref
    reference_fit$group_cleanup_pseudocount <- reference_group_cleanup_pseudocount

    # Refresh reference-mass diagnostics after cleanup changes the denominator.
    mass_summary <- summarize_reference_mass(
      X = filtered_otu,
      ref_idx = final_reference_idx,
      include_zero_count = FALSE
    )

    reference_fit$reference_mass <- mass_summary$reference_mass
    reference_fit$reference_mass_min <- mass_summary$reference_mass_min
    reference_fit$reference_mass_median <- mass_summary$reference_mass_median
    reference_fit$reference_mass_max <- mass_summary$reference_mass_max

    if (verbose) {
      message("Group-aware reference cleanup complete.")
      message("  Removed references: ", length(cleanup_obj$removed_references))
      message("  Final reference size: ", length(final_reference_idx))
      message(
        "  Reference mass (min/median/max): ",
        round(reference_fit$reference_mass_min, 3), " / ",
        round(reference_fit$reference_mass_median, 3), " / ",
        round(reference_fit$reference_mass_max, 3)
      )
    }
  }

  if (verbose) message(strrep("-", 50))

  # Build the paired prevalence and abundance responses against the final references.
  responses <- build_reference_normalized_responses(
    X = filtered_otu,
    ref_idx = final_reference_idx,
    return_prevalence = TRUE,
    return_abundance = TRUE,
    verbose = verbose
  )

  if (verbose) message(strrep("-", 50))

  # Convert prevalence responses into the expected-rarefied matrix used by the model arm.
  prevalence_inputs <- build_prevalence_inputs_by_engine(
    otu = filtered_otu,
    taxon_idx = responses$taxon_idx,
    metadata = meta,
    formula = formula,
    prevalence_expected_rarefied_eps = prevalence_expected_rarefied_eps
  )

  prevalence_matrix_fit <- prevalence_inputs$prevalence_matrix
  prevalence_formula_fit <- prevalence_inputs$prevalence_formula

  # Fit observed prevalence statistics before the permutation engine reuses them.
  prevalence_fit <- fit_prevalence_expected_rarefied_arm(
    prevalence_matrix = prevalence_matrix_fit,
    metadata = meta,
    formula = prevalence_formula_fit,
    tested_term = tested_term,
    return_models = keep_diagnostics,
    verbose = verbose
  )

  if (verbose) message(strrep("-", 50))

  # Fit observed abundance statistics on reference-normalized abundance responses.
  abundance_fit <- fit_abundance_arm(
    abundance_matrix = responses$abundance_matrix,
    metadata = meta,
    formula = formula,
    tested_term = tested_term,
    min_nonzero_n = abundance_min_nonzero_n,
    return_models = keep_diagnostics,
    verbose = verbose
  )

  if (verbose) message(strrep("-", 50))

  # Generate shared Freedman--Lane nulls for both components.
  permutation_obj <- run_component_permutations(
    prevalence_matrix = prevalence_matrix_fit,
    abundance_matrix = responses$abundance_matrix,
    metadata = meta,
    formula = formula,
    prevalence_formula = prevalence_formula_fit,
    tested_term = tested_term,
    tested_vars = tested_vars,
    observed_prevalence_fit = prevalence_fit,
    observed_abundance_fit = abundance_fit,
    n_perm = n_perm,
    cluster_id = cluster_id,
    strata = strata,
    cluster_exposure_mode = cluster_exposure_mode,
    prevalence_fit_args = list(
      return_models = FALSE,
      verbose = FALSE
    ),
    abundance_fit_args = list(
      min_nonzero_n = abundance_min_nonzero_n,
      resid_winsorize = abund_resid_winsorize,
      resid_winsor_quantile = abund_resid_winsor_quantile,
      resid_winsor_min_n = abund_resid_winsor_min_n,
      return_models = FALSE,
      verbose = FALSE
    ),
    neutral_stat = 0,
    parallelize_permutations = parallelize_permutations,
    permutation_n_cores = permutation_n_cores,
    seed = seed,
    verbose = verbose
  )

  # Convert observed and null component statistics into empirical p-values.
  component_p_table <- compute_empirical_component_pvalues(permutation_obj)

  # Combine component evidence and optionally apply adaptive independent filtering.
  union_layer <- compute_union_layer(
    perm_obj = permutation_obj,
    component_p_table = component_p_table,
    acat_weights = union_acat_weights,
    independent_filtering = independent_filtering,
    independent_filter_alpha = independent_filter_alpha,
    independent_filter_theta_grid = independent_filter_theta_grid,
    independent_filter_filtered_q_value = independent_filter_filtered_q_value
  )

  results_table <- union_layer$result_table

  # Remap modeled taxon indices from the filtered matrix back to the original OTU matrix.
  results_table$taxon_index <- filter_obj$kept_idx[results_table$taxon_index]
  results_table$reason_NA <- NA_character_

  all_taxon_names <- colnames(otu)
  kept_taxon_names <- colnames(filtered_otu)
  filtered_taxon_names <- all_taxon_names[filter_obj$dropped_idx]
  filtered_taxon_index <- filter_obj$dropped_idx

  reference_taxon_names <- final_reference_names
  reference_taxon_index <- filter_obj$kept_idx[match(reference_taxon_names, kept_taxon_names)]

  # Add NA-statistic rows for taxa excluded from modeling as references or pre-filters.
  make_placeholder_rows <- function(taxon_names, taxon_index, reason, template) {
    n_add <- length(taxon_names)

    if (n_add == 0L) {
      return(template[FALSE, , drop = FALSE])
    }

    # Preserve all result columns by cloning the template shape before filling identifiers.
    placeholder_rows <- template[rep(NA_integer_, n_add), , drop = FALSE]
    placeholder_rows$taxon_name <- taxon_names
    placeholder_rows$taxon_index <- taxon_index
    placeholder_rows$reason_NA <- reason

    placeholder_rows
  }

  reference_rows <- make_placeholder_rows(
    taxon_names = reference_taxon_names,
    taxon_index = reference_taxon_index,
    reason = "reference",
    template = results_table
  )

  filtered_rows <- make_placeholder_rows(
    taxon_names = filtered_taxon_names,
    taxon_index = filtered_taxon_index,
    reason = "filtered",
    template = results_table
  )

  # Return modeled, reference, and filtered taxa in the original feature order.
  results_table <- rbind(results_table, reference_rows, filtered_rows)
  results_table <- results_table[order(results_table$taxon_index), , drop = FALSE]
  rownames(results_table) <- NULL

  # Move the main p-value, q-value, and exclusion fields to the front.
  preferred_cols <- c(
    "taxon_index",
    "taxon_name",
    "prev_empirical_p",
    "abund_empirical_p",
    "union_p_value",
    "union_q_value",
    "rejected_bh_0.05",
    "reason_NA"
  )

  existing_preferred <- preferred_cols[preferred_cols %in% names(results_table)]
  remaining_cols <- setdiff(names(results_table), existing_preferred)
  results_table <- results_table[, c(existing_preferred, remaining_cols), drop = FALSE]

  # Keep the diagnostic union object but avoid returning the pre-placeholder result table twice.
  union_layer$result_table <- NULL

  # Return only the main table unless intermediate objects were requested.
  if (!isTRUE(keep_diagnostics)) {
    out <- list(
      results_table = results_table
    )

    if (verbose) {
      message(strrep("-", 50))
      message("Finished!")
    }

    return(out)
  }

  # Return the full workflow state for debugging, inspection, or method development.
  out <- list(
    call = match.call(),
    filter = filter_obj,
    reference = list(
      selection = reference_fit
    ),
    responses = responses,
    observed = list(
      prevalence = prevalence_fit,
      abundance = abundance_fit
    ),
    permutations = permutation_obj,
    union = union_layer,
    results_table = results_table
  )

  if (verbose) {
    message(strrep("-", 50))
    message("Finished!")
  }

  out
}
