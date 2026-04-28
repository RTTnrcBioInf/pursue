# -----------------------------------------------------------------------------
# Response input checks and reference summaries
# -----------------------------------------------------------------------------

# Validate the count matrix and normalize reference indices
check_response_inputs <- function(
    X,
    ref_idx
) {
  X <- as_numeric_matrix(X, "X")

  if (nrow(X) < 1L || ncol(X) < 1L) {
    stop("X must have at least 1 row and 1 column.")
  }

  if (any(!is.finite(X))) {
    stop("X contains non-finite values.")
  }

  if (any(X < 0)) {
    stop("X must be non-negative.")
  }

  if (max(abs(X - round(X))) > 1e-8) {
    warning("X is not integer-valued within tolerance; proceeding anyway.")
  }

  # Accept logical reference masks as a convenience, then work with indices.
  if (is.logical(ref_idx)) {
    if (length(ref_idx) != ncol(X)) {
      stop("Logical ref_idx must have length ncol(X).")
    }

    ref_idx <- which(ref_idx)
  }

  ref_idx <- as.integer(unique(ref_idx))

  if (length(ref_idx) < 1L) {
    stop("ref_idx must contain at least one column index.")
  }

  if (any(ref_idx < 1L | ref_idx > ncol(X))) {
    stop("ref_idx contains invalid column indices.")
  }

  # Ensure downstream matrices and summaries can always carry stable names.
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("Taxon_", seq_len(ncol(X)))
  }

  if (is.null(rownames(X))) {
    rownames(X) <- paste0("Sample_", seq_len(nrow(X)))
  }

  taxon_idx <- setdiff(seq_len(ncol(X)), ref_idx)

  if (length(taxon_idx) < 1L) {
    stop("No taxa left after excluding reference taxa.")
  }

  list(
    X = X,
    ref_idx = ref_idx,
    taxon_idx = taxon_idx
  )
}

# Summarize per-sample mass contributed by the selected reference set
summarize_reference_mass <- function(
    X,
    ref_idx,
    include_zero_count = TRUE
) {
  reference_mass <- rowSums(X[, ref_idx, drop = FALSE])

  out <- list(
    reference_mass = reference_mass,
    reference_mass_min = unname(min(reference_mass)),
    reference_mass_median = unname(stats::median(reference_mass)),
    reference_mass_max = unname(max(reference_mass))
  )

  if (include_zero_count) {
    out$n_zero_reference_mass <- sum(reference_mass == 0)
  }

  out
}


# -----------------------------------------------------------------------------
# Prevalence response construction
# -----------------------------------------------------------------------------

# Convert non-reference counts into binary observed-prevalence responses
build_prevalence_matrix <- function(X, taxon_idx) {
  prevalence_matrix <- (X[, taxon_idx, drop = FALSE] > 0) * 1L
  storage.mode(prevalence_matrix) <- "integer"

  positive_count_per_taxon <- colSums(prevalence_matrix)
  prevalence_per_taxon <- positive_count_per_taxon / nrow(prevalence_matrix)

  list(
    prevalence_matrix = prevalence_matrix,
    positive_count_per_taxon = positive_count_per_taxon,
    prevalence_per_taxon = prevalence_per_taxon
  )
}

# Use the minimum sample depth as the common rarefaction target
choose_min_rarefaction_depth <- function(otu_full_filtered) {
  otu_full_filtered <- as_numeric_matrix(otu_full_filtered, "otu_full_filtered")
  full_depth <- rowSums(otu_full_filtered)
  rarefaction_depth <- as.integer(min(full_depth))

  if (!is.finite(rarefaction_depth) || rarefaction_depth < 1L) {
    stop("Minimum prevalence rarefaction depth must be >= 1.")
  }

  rarefaction_depth
}

# Compute expected taxon presence after rarefying each sample to a common depth
compute_expected_rarefied_prevalence <- function(
    otu_full_filtered,
    rarefaction_depth,
    eps = 1e-12
) {
  otu_full_filtered <- as_numeric_matrix(otu_full_filtered, "otu_full_filtered")

  rarefaction_depth <- as.integer(rarefaction_depth)
  full_depth <- rowSums(otu_full_filtered)

  if (rarefaction_depth > min(full_depth)) {
    stop("Expected-rarefied prevalence requires rarefaction_depth <= the minimum sample depth.")
  }

  prevalence_matrix <- matrix(
    0,
    nrow = nrow(otu_full_filtered),
    ncol = ncol(otu_full_filtered),
    dimnames = dimnames(otu_full_filtered)
  )

  for (sample_pos in seq_len(nrow(otu_full_filtered))) {
    sample_depth <- full_depth[sample_pos]
    sample_counts <- otu_full_filtered[sample_pos, ]
    present_taxa <- which(sample_counts > 0)

    if (length(present_taxa) == 0L) {
      next
    }

    # Hypergeometric complement: probability at least one read survives rarefaction.
    log_denom <- lchoose(sample_depth, rarefaction_depth)
    log_num <- lchoose(sample_depth - sample_counts[present_taxa], rarefaction_depth)

    expected_presence <- 1 - exp(log_num - log_denom)
    expected_presence <- pmin(pmax(expected_presence, eps), 1 - eps)

    prevalence_matrix[sample_pos, present_taxa] <- expected_presence
  }

  list(
    prevalence_matrix = prevalence_matrix,
    positive_count_per_taxon = colSums(otu_full_filtered > 0),
    prevalence_per_taxon = colMeans(otu_full_filtered > 0),
    rarefaction_depth = rarefaction_depth,
    full_depth = full_depth
  )
}

# Remove explicit depth adjustment from the prevalence model formula
remove_depth_from_prevalence_formula <- function(formula) {
  terms_obj <- stats::terms(formula)
  term_labels <- attr(terms_obj, "term.labels")
  term_labels <- setdiff(term_labels, "log_depth")

  right_side <- if (length(term_labels) == 0L) {
    "1"
  } else {
    paste(term_labels, collapse = " + ")
  }

  stats::as.formula(paste("~", right_side), env = environment(formula))
}

# Build prevalence inputs expected by the prevalence fitting engine
build_prevalence_inputs_by_engine <- function(
    otu,
    taxon_idx,
    metadata,
    formula,
    prevalence_expected_rarefied_eps = 1e-12
) {
  # Prevalence is depth-normalized through expected rarefaction, not log-depth adjustment.
  prevalence_formula_nondepth <- remove_depth_from_prevalence_formula(
    formula = formula
  )

  # Rarefy on the full filtered table before subsetting returned taxa.
  rarefaction_depth <- choose_min_rarefaction_depth(
    otu_full_filtered = otu
  )

  prevalence_obj <- compute_expected_rarefied_prevalence(
    otu_full_filtered = otu,
    rarefaction_depth = rarefaction_depth,
    eps = prevalence_expected_rarefied_eps
  )

  list(
    prevalence_matrix = prevalence_obj$prevalence_matrix[, taxon_idx, drop = FALSE],
    prevalence_formula = prevalence_formula_nondepth,
    rarefaction_depth = rarefaction_depth,
    prevalence_mode = "expected_rarefied_min_all_samples"
  )
}


# -----------------------------------------------------------------------------
# Abundance response construction
# -----------------------------------------------------------------------------

# Build log count-to-reference abundance responses for returned taxa
build_abundance_matrix <- function(
    X,
    taxon_idx,
    ref_idx,
    ref_mass
) {
  taxon_counts <- X[, taxon_idx, drop = FALSE]
  n_samples <- nrow(taxon_counts)
  n_taxa <- ncol(taxon_counts)
  returned_is_reference_taxon <- taxon_idx %in% ref_idx

  abundance_matrix <- matrix(
    NA_real_,
    nrow = n_samples,
    ncol = n_taxa,
    dimnames = list(rownames(taxon_counts), colnames(taxon_counts))
  )

  # Start each denominator from the per-sample reference mass.
  denominator_matrix <- matrix(
    ref_mass,
    nrow = n_samples,
    ncol = n_taxa,
    dimnames = list(rownames(taxon_counts), colnames(taxon_counts))
  )

  # Reference taxa are normalized leave-one-out to avoid self-comparison.
  if (any(returned_is_reference_taxon)) {
    denominator_matrix[, returned_is_reference_taxon] <-
      denominator_matrix[, returned_is_reference_taxon, drop = FALSE] -
      taxon_counts[, returned_is_reference_taxon, drop = FALSE]
  }

  # Only positive taxon counts with positive reference denominators yield log-ratios.
  has_taxon_count <- taxon_counts > 0
  valid_denominator <- denominator_matrix > 0
  usable_cells <- has_taxon_count & valid_denominator

  if (any(usable_cells)) {
    abundance_matrix[usable_cells] <-
      log(taxon_counts[usable_cells]) -
      log(denominator_matrix[usable_cells])
  }

  # Track missingness and invalid-denominator diagnostics for downstream reporting.
  abundance_n_observed_per_taxon <- colSums(!is.na(abundance_matrix))
  abundance_na_count <- sum(is.na(abundance_matrix))
  abundance_inf_count <- sum(is.infinite(abundance_matrix))
  abundance_nan_count <- sum(is.nan(abundance_matrix))

  zero_taxon_count_mask <- taxon_counts == 0
  nonpositive_den_mask <- !valid_denominator & has_taxon_count
  leave_one_out_den_mask <- matrix(
    FALSE,
    nrow = n_samples,
    ncol = n_taxa,
    dimnames = list(rownames(taxon_counts), colnames(taxon_counts))
  )

  if (any(returned_is_reference_taxon)) {
    leave_one_out_den_mask[, returned_is_reference_taxon] <- TRUE
  }

  list(
    abundance_matrix = abundance_matrix,
    abundance_n_observed_per_taxon = abundance_n_observed_per_taxon,
    abundance_na_count = abundance_na_count,
    abundance_inf_count = abundance_inf_count,
    abundance_nan_count = abundance_nan_count,
    n_zero_taxon_count_cells = sum(zero_taxon_count_mask),
    n_nonpositive_den_cells = sum(nonpositive_den_mask),
    returned_is_reference_taxon = returned_is_reference_taxon,
    n_reference_taxa_leave_one_out = sum(returned_is_reference_taxon),
    n_nonpositive_leave_one_out_den_cells = sum(nonpositive_den_mask & leave_one_out_den_mask)
  )
}


# -----------------------------------------------------------------------------
# Combined response object
# -----------------------------------------------------------------------------

# Construct reference-normalized prevalence and abundance responses from a count table
build_reference_normalized_responses <- function(
    X,
    ref_idx,
    return_prevalence = TRUE,
    return_abundance = TRUE,
    verbose = FALSE
) {
  # Validate the input matrix and resolve the taxa returned for testing.
  input_info <- check_response_inputs(
    X = X,
    ref_idx = ref_idx
  )

  X <- input_info$X
  ref_idx <- input_info$ref_idx
  taxon_idx <- input_info$taxon_idx

  if (!is.logical(return_prevalence) || length(return_prevalence) != 1L) {
    stop("return_prevalence must be TRUE or FALSE.")
  }

  if (!is.logical(return_abundance) || length(return_abundance) != 1L) {
    stop("return_abundance must be TRUE or FALSE.")
  }

  # Compute reference mass once and reuse it across abundance construction and diagnostics.
  reference_summary <- summarize_reference_mass(
    X = X,
    ref_idx = ref_idx,
    include_zero_count = TRUE
  )

  # Construct component responses only when requested by the caller.
  prevalence_obj <- if (return_prevalence) {
    build_prevalence_matrix(X, taxon_idx)
  } else {
    NULL
  }

  abundance_obj <- if (return_abundance) {
    build_abundance_matrix(
      X = X,
      taxon_idx = taxon_idx,
      ref_idx = ref_idx,
      ref_mass = reference_summary$reference_mass
    )
  } else {
    NULL
  }

  # Keep reference and returned-taxon metadata aligned with the generated matrices.
  is_reference_taxon <- seq_len(ncol(X)) %in% ref_idx

  reference_info <- list(
    ref_idx = ref_idx,
    ref_names = colnames(X)[ref_idx],
    taxon_idx = taxon_idx,
    taxon_names = colnames(X)[taxon_idx],
    is_reference_taxon = is_reference_taxon,
    reference_mass = reference_summary$reference_mass,
    reference_mass_min = reference_summary$reference_mass_min,
    reference_mass_median = reference_summary$reference_mass_median,
    reference_mass_max = reference_summary$reference_mass_max,
    n_zero_reference_mass = reference_summary$n_zero_reference_mass
  )

  prevalence_info <- list(
    prevalence_matrix = if (return_prevalence) prevalence_obj$prevalence_matrix else NULL,
    positive_count_per_taxon = if (return_prevalence) prevalence_obj$positive_count_per_taxon else NULL,
    prevalence_per_taxon = if (return_prevalence) prevalence_obj$prevalence_per_taxon else NULL
  )

  abundance_info <- list(
    abundance_matrix = if (return_abundance) abundance_obj$abundance_matrix else NULL,
    abundance_n_observed_per_taxon = if (return_abundance) abundance_obj$abundance_n_observed_per_taxon else NULL,
    abundance_na_count = if (return_abundance) abundance_obj$abundance_na_count else NA_integer_,
    abundance_inf_count = if (return_abundance) abundance_obj$abundance_inf_count else NA_integer_,
    abundance_nan_count = if (return_abundance) abundance_obj$abundance_nan_count else NA_integer_,
    n_zero_taxon_count_cells = if (return_abundance) abundance_obj$n_zero_taxon_count_cells else NA_integer_,
    n_nonpositive_den_cells = if (return_abundance) abundance_obj$n_nonpositive_den_cells else NA_integer_
  )

  run_info <- list(
    n_samples = nrow(X),
    n_taxa_total = ncol(X),
    n_taxa_returned = length(taxon_idx)
  )

  out <- c(
    reference_info,
    prevalence_info,
    abundance_info,
    run_info
  )

  if (verbose) {
    message("Reference-normalized response construction complete.")
    message("  Returned taxa: ", out$n_taxa_returned, " / ", out$n_taxa_total)
    message("  Reference size: ", length(out$ref_idx))
    message("  Zero reference-mass samples: ", out$n_zero_reference_mass)
    if (return_abundance) {
      message("  Abundance NA cells: ", out$abundance_na_count)
      message("  Abundance Inf cells: ", out$abundance_inf_count)
    }
  }

  out
}
