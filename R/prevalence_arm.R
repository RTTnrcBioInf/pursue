# -----------------------------------------------------------------------------
# Expected-rarefied prevalence arm
# -----------------------------------------------------------------------------

# Fit the expected-rarefied prevalence arm through the abundance-model interface
fit_prevalence_expected_rarefied_arm <- function(
    prevalence_matrix,
    metadata,
    formula,
    tested_term,
    taxa_subset = NULL,
    return_models = FALSE,
    verbose = FALSE
) {
  # Reuse the abundance fitting path, treating expected-rarefied prevalence as the response.
  abundance_fit <- fit_abundance_arm(
    abundance_matrix = prevalence_matrix,
    metadata = metadata,
    formula = formula,
    tested_term = tested_term,
    taxa_subset = taxa_subset,
    min_nonzero_n = 1L,
    return_models = return_models,
    verbose = verbose,
    progress_label = "Prevalence arm"
  )

  result_table <- abundance_fit$result_table
  prevalence_values <- as_numeric_matrix(prevalence_matrix, "prevalence_matrix")
  taxon_indices <- result_table$taxon_index

  n_present <- integer(nrow(result_table))
  n_absent <- integer(nrow(result_table))

  # Track prevalence-specific support diagnostics for each fitted taxon.
  for (row_pos in seq_len(nrow(result_table))) {
    taxon_prevalence <- prevalence_values[, taxon_indices[row_pos]]
    finite_rows <- is.finite(taxon_prevalence)

    n_present[row_pos] <- sum(taxon_prevalence[finite_rows] > 0)
    n_absent[row_pos] <- sum(taxon_prevalence[finite_rows] <= 0)
  }

  # Preserve the abundance-arm result shape while marking the prevalence-scale model.
  result_table$n_present <- n_present
  result_table$n_absent <- n_absent
  result_table$model_type <- "expected_rarefied_lm"

  out <- list(
    result_table = result_table,
    fits = abundance_fit$fits,
    fit_summaries = abundance_fit$fit_summaries,
    tested_term = tested_term,
    formula = formula,
    n_taxa_fit = abundance_fit$n_taxa_fit
  )

  class(out) <- c("pursue_prevalence_arm", "list")
  out
}


# -----------------------------------------------------------------------------
# Expected-rarefied prevalence permutation cache
# -----------------------------------------------------------------------------

# Prepare prevalence-scale permutation inputs through the abundance cache builder
prepare_expected_rarefied_prevalence_cache <- function(
    prevalence_matrix,
    metadata,
    formula,
    tested_term,
    tested_vars = NULL,
    observed_fit = NULL,
    fit_args = list(),
    neutral_stat = 0,
    perm_design,
    verbose = FALSE
) {
  if (missing(perm_design) || is.null(perm_design)) {
    stop("prepare_expected_rarefied_prevalence_cache() requires perm_design.")
  }

  if (!inherits(perm_design, "pursue_permutation_design")) {
    stop("perm_design must inherit from 'pursue_permutation_design'.")
  }

  # Use abundance-cache machinery but force prevalence-appropriate fit settings.
  cache <- prepare_abundance_permutation_cache(
    abundance_matrix = prevalence_matrix,
    metadata = metadata,
    formula = formula,
    tested_term = tested_term,
    tested_vars = tested_vars,
    observed_fit = observed_fit,
    fit_args = c(
      fit_args,
      list(
        min_nonzero_n = 1L,
        resid_winsorize = FALSE
      )
    ),
    neutral_stat = neutral_stat,
    perm_design = perm_design,
    verbose = verbose,
    cache_label = "Prevalence"
  )

  # Relabel the cache so downstream code can distinguish the prevalence null path.
  cache$statistic_type <- "expected_rarefied_linear_perm"
  class(cache) <- c(
    "pursue_prevalence_expected_rarefied_cache",
    "pursue_abundance_permutation_cache",
    "list"
  )

  cache
}


# -----------------------------------------------------------------------------
# Expected-rarefied prevalence permutation null
# -----------------------------------------------------------------------------

# Compute the expected-rarefied prevalence null through the abundance null engine
compute_expected_rarefied_prevalence_null <- function(
    prev_cache,
    perm_list,
    perm_design,
    neutral_stat = 0,
    parallelize_permutations = FALSE,
    permutation_n_cores = max(1L, parallel::detectCores() - 1L),
    verbose = FALSE
) {
  compute_abundance_permutation_null(
    abund_cache = prev_cache,
    perm_list = perm_list,
    perm_design = perm_design,
    neutral_stat = neutral_stat,
    parallelize_permutations = parallelize_permutations,
    permutation_n_cores = permutation_n_cores,
    verbose = verbose,
    null_label = "Prevalence"
  )
}
