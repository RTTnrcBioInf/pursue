#' Run the PURSUE workflow
#'
#' `run_pursue()` fits the full PURSUE workflow for microbial differential
#' abundance analysis. The procedure optionally filters low information features,
#' selects a stable reference set, constructs reference-normalized abundance and
#' rarefied prevalence responses, fits the two model components separately, computes
#' empirical p-values using a Freedman-Lane residual permutation scheme,
#' optionally applies residual winsorization within the abundance permutation
#' engine, and combines component evidence into union level results with
#' optional independent filtering.
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
#' @param require_no_zero_mass Logical; if `TRUE`, requires that no sample have
#'   zero reference mass at the selected prefix. This is a hard correctness gate:
#'   zero reference mass produces undefined log-ratios. Default is `TRUE`.
#' @param min_bulk_count Minimum absolute reference count required at the lower
#'   tail of the sample distribution. Default is `10`.
#' @param bulk_count_quantile Lower-tail quantile at which `min_bulk_count` must
#'   be cleared. Default is `0.05` (i.e., the 5th percentile).
#' @param min_bulk_rel Minimum reference mass as a fraction of sample depth,
#'   required at the lower tail of the sample distribution. Default is `0.01`.
#' @param bulk_rel_quantile Lower-tail quantile at which `min_bulk_rel` must be
#'   cleared. Default is `0.10` (i.e., the 10th percentile).
#' @param hostage_prevention Logical; if `TRUE`, extends the selected reference
#'   prefix when needed to reduce fragile low-mass tails.
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
#' @param depth_adjust Logical; if `TRUE`, appends `log(depth)` to the model
#'   formula. Default is `TRUE`. **NOTE** that this is used by reference
#'   selection and the abundance component, not by the prevalence component.
#'   The prevalence component removes any column named `log_depth` from the
#'   formula.
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
#'   sample or cluster level. One of `"sample"` or `"cluster"`.
#' @param seed Optional random seed.
#' @param independent_filtering Logical; if `TRUE`, applies independent
#'   filtering before multiple-testing correction.
#' @param q_alpha Numeric scalar in `(0, 1)`. Q-value/FDR cutoff used to
#'   label BH rejections in the raw and final/independent-filtered result
#'   layers. Default is `0.05`.
#' @param independent_filter_quantile Quantile of the filter statistic below
#'   which low-information taxa are excluded before BH correction. Must be in
#'   [0, 1). Default is 0.20.
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
#' a Freedman-Lane residual permutation scheme. Component evidence is
#' then combined into union results, with optional independent
#' filtering before multiple testing correction.
#'
#' @return
#' If `keep_diagnostics = FALSE`, a list with `results_table`. If
#' `keep_diagnostics = TRUE`, a list containing the final `results_table` plus
#' the filter object, reference selection output, constructed responses,
#' cache-derived observed component table, permutation outputs, and union
#' layer output.
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
#'   bootstrap_B = 10,
#'   n_perm = 99,
#'   verbose = TRUE
#' )
#'
#' head(out$results_table)
#'
#' @export

# Top level runner and user entry point
# Calls module functions in sequence
run_pursue <- function(
    otu,
    meta,
    formula,
    tested_term,
    tested_vars = NULL,

    min_prevalence_filter = NULL,
    min_total_count_filter = NULL,

    reference_pseudocount = 1,
    ref_n_max = 1000L,
    bootstrap_B = 100L,
    bootstrap_stability_threshold = 0.80,
    require_no_zero_mass = TRUE,
    min_bulk_count = 10,
    bulk_count_quantile = 0.05,
    min_bulk_rel = 0.01,
    bulk_rel_quantile = 0.10,
    hostage_prevention = TRUE,
    hostage_max_size_inflation = 0.50,
    hostage_max_score_inflation = 0.25,
    hostage_max_stability_drop = 0.05,
    reference_group_cleanup = FALSE,
    reference_group_cleanup_frac = 0.15,
    reference_group_cleanup_min_ref = 10L,

    depth_adjust = TRUE,

    prevalence_expected_rarefied_eps = 1e-12,
    parallelize_permutations = FALSE,
    permutation_n_cores = max(1L, parallel::detectCores() - 1L),

    abundance_min_nonzero_n = 3L,
    abund_resid_winsorize = FALSE,
    abund_resid_winsor_quantile = 0.05,
    abund_resid_winsor_min_n = 20L,

    n_perm = 1999,
    cluster_id = NULL,
    strata = NULL,
    cluster_exposure_mode = c("sample", "cluster"),
    seed = 1,

    q_alpha = 0.05,

    independent_filtering = TRUE,
    independent_filter_quantile = 0.20,
    independent_filter_filtered_q_value = NA_real_,

    keep_diagnostics = FALSE,
    verbose = TRUE
) {

  # input coercion
  otu <- as_numeric_matrix(otu, "otu")
  meta <- as_plain_data_frame(meta)

  # otu table checks
  if (any(!is.finite(otu))) stop("otu must contain only finite counts.")
  if (any(otu < 0)) stop("otu must contain non-negative counts.")
  if (nrow(otu) < 2L) stop("otu must contain at least two samples.")
  if (ncol(otu) < 2L) stop("otu must contain at least two taxa.")
  if (any(rowSums(otu) <= 0)) stop("Each sample in otu must have positive total counts.")

  # parameter checks
  check_scalar <- function(x, name, lo = -Inf, hi = Inf,
                           lo_open = FALSE, hi_open = FALSE,
                           integer = FALSE, allow_null = FALSE) {
    if (allow_null && is.null(x)) return(invisible())
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        (lo_open && x <= lo) || (!lo_open && x < lo) ||
        (hi_open && x >= hi) || (!hi_open && x > hi) ||
        (integer && x != as.integer(x))) {
      stop(name, " must be ",
           if (allow_null) "NULL or ",
           "a single ", if (integer) "integer" else "numeric",
           " value in ",
           if (lo_open) "(" else "[", lo, ", ", hi, if (hi_open) ")" else "]",
           ".", call. = FALSE)
    }
  }

  check_scalar(reference_pseudocount,         "reference_pseudocount", lo = 0, lo_open = TRUE)
  check_scalar(bootstrap_B,                   "bootstrap_B", lo = 1, integer = TRUE)
  check_scalar(bootstrap_stability_threshold, "bootstrap_stability_threshold", lo = 0, hi = 1)
  check_scalar(bulk_count_quantile,           "bulk_count_quantile",           lo = 0, hi = 1)
  check_scalar(bulk_rel_quantile,             "bulk_rel_quantile",             lo = 0, hi = 1)
  check_scalar(min_prevalence_filter,         "min_prevalence_filter",         lo = 0, hi = 1, allow_null = TRUE)
  check_scalar(min_total_count_filter,        "min_total_count_filter",        lo = 0,         allow_null = TRUE)
  check_scalar(reference_group_cleanup_frac,  "reference_group_cleanup_frac",  lo = 0, hi = 1, hi_open = TRUE)
  check_scalar(abund_resid_winsor_quantile,   "abund_resid_winsor_quantile",   lo = 0, hi = 0.5, lo_open = TRUE, hi_open = TRUE)
  check_scalar(independent_filter_quantile,   "independent_filter_quantile",   lo = 0, hi = 1, hi_open = TRUE)
  check_scalar(q_alpha,                       "q_alpha",                       lo = 0, hi = 1, lo_open = TRUE, hi_open = TRUE)

  # input ID checks
  if (nrow(otu) != nrow(meta)) {
    stop("nrow(otu) must equal nrow(meta).")
  }

  # fill missing IDs
  if (is.null(rownames(otu))) {
    rownames(otu) <- paste0("Sample_", seq_len(nrow(otu)))
  }

  if (is.null(colnames(otu))) {
    colnames(otu) <- paste0("Taxon_", seq_len(ncol(otu)))
  }

  if (is.null(rownames(meta))) {
    rownames(meta) <- rownames(otu)
  }

  # reorder meta rows to match otu if they're a permutation
  if (!identical(rownames(otu), rownames(meta))) {
    if (setequal(rownames(otu), rownames(meta))) {
      meta <- meta[rownames(otu), , drop = FALSE]
    } else {
      stop("Row names of otu and meta do not match.")
    }
  }

  # Caller: depth adjustment as a covariate (when requested)
  if (depth_adjust) {
    meta <- add_log_depth(
      otu = otu,
      metadata = meta,
      depth_var_name = "log_depth"
    )

    formula <- add_term_if_missing(formula, "log_depth")
  }

  cluster_exposure_mode <- match.arg(cluster_exposure_mode)

  # Caller: filter taxa
  filter_obj <- filter_taxa(
    otu = otu,
    min_prevalence = min_prevalence_filter,
    min_total_count = min_total_count_filter,
    verbose = verbose
  )

  filtered_otu <- filter_obj$otu

  if (verbose) {
    message(strrep("-", 50))
    message("Bootstrap reference selection in progress...")
  }

  # Caller: reference selection
  reference_fit <- select_reference_bootstrap_prefix(
    X = filtered_otu,
    metadata = meta,
    formula = formula,
    tested_term = tested_term,
    p_max = ref_n_max,
    pseudocount = reference_pseudocount,
    bootstrap_B = bootstrap_B,
    bootstrap_stability_threshold = bootstrap_stability_threshold,
    require_no_zero_mass = require_no_zero_mass,
    min_bulk_count = min_bulk_count,
    bulk_count_quantile = bulk_count_quantile,
    min_bulk_rel = min_bulk_rel,
    bulk_rel_quantile = bulk_rel_quantile,
    hostage_prevention = hostage_prevention,
    hostage_max_size_inflation = hostage_max_size_inflation,
    hostage_max_score_inflation = hostage_max_score_inflation,
    hostage_max_stability_drop = hostage_max_stability_drop,
    return_ratio_matrix = keep_diagnostics,
    seed = seed,
    verbose = verbose
  )

  final_reference_idx <- reference_fit$selected_references # pass to group cleanup
  final_reference_names <- reference_fit$selected_reference_names

  # Caller: group shift cleanup
  if (isTRUE(reference_group_cleanup)) {
    cleanup_obj <- cleanup_reference_group_shift(
      X = filtered_otu,
      ref_idx = final_reference_idx,
      metadata = meta,
      tested_term = tested_term,
      frac = reference_group_cleanup_frac,
      min_ref = reference_group_cleanup_min_ref,
      pseudocount = reference_pseudocount
    )

    final_reference_idx <- cleanup_obj$selected_references
    final_reference_names <- colnames(filtered_otu)[final_reference_idx]

    # stash cleanup diagnostics on the reference fit
    reference_fit$selected_references <- final_reference_idx
    reference_fit$selected_reference_names <- final_reference_names
    reference_fit$group_cleanup_removed_references <- cleanup_obj$removed_references
    reference_fit$group_cleanup_removed_names <- colnames(filtered_otu)[cleanup_obj$removed_references]
    reference_fit$group_cleanup_scores <- cleanup_obj$group_cleanup_scores
    reference_fit$group_cleanup_frac <- reference_group_cleanup_frac
    reference_fit$group_cleanup_min_ref <- reference_group_cleanup_min_ref
    reference_fit$group_cleanup_pseudocount <- reference_pseudocount

    # refresh mass diagnostics (denominator changed)
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

  # Caller: response construction (abundance matrix)
  responses <- build_reference_normalized_responses(
    X = filtered_otu,
    ref_idx = final_reference_idx,
    verbose = verbose
  )

  if (verbose) message(strrep("-", 50))

  # filtered-table indices are used internally after taxon filtering
  # original-table indices are restored for final reporting
  modeled_filtered_idx <- responses$taxon_idx
  modeled_original_idx <- filter_obj$kept_idx[modeled_filtered_idx]

  # The prevalence component wants the original count subset; depth still
  # comes from the full table.
  otu_prevalence_taxon_counts <- otu[, modeled_original_idx, drop = FALSE]

  # keep prevalence cols aligned with abundance
  colnames(otu_prevalence_taxon_counts) <- responses$taxon_names
  rownames(otu_prevalence_taxon_counts) <- rownames(filtered_otu)

  # build rarefied prevalence responses for modeled taxa
  # rarefaction depth and full sample depth are based on the original count table
  prevalence_inputs <- build_prevalence_inputs_by_engine(
    otu_taxon_counts = otu_prevalence_taxon_counts,
    otu_depth_universe = otu,
    metadata = meta,
    formula = formula,
    prevalence_expected_rarefied_eps = prevalence_expected_rarefied_eps
  )

  # matrix and formula (with any depth terms removed)
  prevalence_matrix_fit <- prevalence_inputs$prevalence_matrix
  prevalence_formula_fit <- prevalence_inputs$prevalence_formula

  if (!identical(colnames(prevalence_matrix_fit), colnames(responses$abundance_matrix))) {
    stop("Prevalence and abundance response columns are not aligned.")
  }

  if (!identical(rownames(prevalence_matrix_fit), rownames(responses$abundance_matrix))) {
    stop("Prevalence and abundance response rows are not aligned.")
  }

  # stash rarefaction bookkeeping on responses
  responses$modeled_filtered_idx <- modeled_filtered_idx
  responses$modeled_original_idx <- modeled_original_idx
  responses$prevalence_rarefaction_depth <- prevalence_inputs$rarefaction_depth
  responses$prevalence_sample_depth <- prevalence_inputs$sample_depth
  responses$prevalence_mode <- prevalence_inputs$prevalence_mode

  # Caller: shared QR/FL component inference
  permutation_obj <- run_component_permutations(
    prevalence_matrix = prevalence_matrix_fit,
    abundance_matrix = responses$abundance_matrix,
    metadata = meta,
    formula = formula,
    prevalence_formula = prevalence_formula_fit,
    tested_term = tested_term,
    tested_vars = tested_vars,
    n_perm = n_perm,
    cluster_id = cluster_id,
    strata = strata,
    cluster_exposure_mode = cluster_exposure_mode,
    abundance_fit_args = list(
      min_nonzero_n = abundance_min_nonzero_n,
      resid_winsorize = abund_resid_winsorize,
      resid_winsor_quantile = abund_resid_winsor_quantile,
      resid_winsor_min_n = abund_resid_winsor_min_n
    ),
    neutral_stat = 0,
    parallelize_permutations = parallelize_permutations,
    permutation_n_cores = permutation_n_cores,
    seed = seed,
    verbose = verbose
  )

  # empirical p per component
  component_p_table <- compute_empirical_component_pvalues(
    permutation_obj,
    verbose = verbose
  )

  # IF stat
  #   F_j = log(1 + sum_i X_ij)
  independent_filter_stat <- log1p(
    colSums(filtered_otu[, responses$taxon_idx, drop = FALSE])
  )
  names(independent_filter_stat) <- responses$taxon_names

  # Caller: ACAT union + optional IF
  union_layer <- compute_union_layer(
    component_p_table = component_p_table,
    acat_weights = c(0.5,0.5), # hard coded equal weights for now
    q_alpha = q_alpha,
    independent_filtering = independent_filtering,
    independent_filter_quantile = independent_filter_quantile,
    independent_filter_stat = independent_filter_stat,
    independent_filter_filtered_q_value = independent_filter_filtered_q_value
  )

  results_table <- union_layer$result_table

  # remap modeled row taxon IDS back to the original table IDS
  results_table$taxon_index_model <- results_table$taxon_index
  # taxon_index is the original table index, taxon_index_model is only modeled taxa
  results_table$taxon_index <- modeled_original_idx[results_table$taxon_index_model]
  results_table$reason_NA <- NA_character_

  all_taxon_names <- colnames(otu)
  kept_taxon_names <- colnames(filtered_otu)
  filtered_taxon_names <- all_taxon_names[filter_obj$dropped_idx]
  filtered_taxon_index <- filter_obj$dropped_idx

  reference_taxon_names <- final_reference_names
  reference_taxon_index <- filter_obj$kept_idx[match(reference_taxon_names, kept_taxon_names)]

  # placeholder rows for excluded taxa (references + filtered).

  reference_rows <- if (length(reference_taxon_names) == 0L) {
    results_table[FALSE, , drop = FALSE]
  } else {
    # Make an NA template with the same columns as results_table,
    # then fill in taxon_name, taxon_index, and reason_NA
    rows <- results_table[rep(NA_integer_, length(reference_taxon_names)), , drop = FALSE]
    rows$taxon_name <- reference_taxon_names
    rows$taxon_index <- reference_taxon_index
    rows$reason_NA <- "reference"
    rows
  }

  filtered_rows <- if (length(filtered_taxon_names) == 0L) {
    results_table[FALSE, , drop = FALSE]
  } else {
    rows <- results_table[rep(NA_integer_, length(filtered_taxon_names)), , drop = FALSE]
    rows$taxon_name <- filtered_taxon_names
    rows$taxon_index <- filtered_taxon_index
    rows$reason_NA <- "filtered"
    rows
  }

  # reassemble in original feature order
  results_table <- rbind(results_table, reference_rows, filtered_rows)
  results_table <- results_table[order(results_table$taxon_index), , drop = FALSE]
  rownames(results_table) <- NULL

  # main fields to the front
  front_cols <- c(
    "taxon_index",
    "taxon_name",
    "taxon_index_model",
    "prev_empirical_p",
    "abund_empirical_p",
    "union_p_value",
    "union_q_value",
    "rejected_bh",
    "rejected_bh_raw",
    "reason_NA"
  )

  existing_preferred <- front_cols[front_cols %in% names(results_table)]
  remaining_cols <- setdiff(names(results_table), existing_preferred)
  results_table <- results_table[, c(existing_preferred, remaining_cols), drop = FALSE]

  # duplicate results table removal
  union_layer$result_table <- NULL

  # results table only return
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

  # + diagnostic return
  out <- list(
    call = match.call(),
    filter = filter_obj,
    reference = list(
      selection = reference_fit
    ),
    responses = responses,
    observed = list(
      components = permutation_obj$observed_component_table
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
