# Function definitions for the two response matrices built off the reference set

# -----------------------------------------------------------------------------
# Response input checks and reference summaries
# -----------------------------------------------------------------------------

# Validate the count matrix and normalize ref_idx into an integer vector.
# Called from build_reference_normalized_responses() once at the top so the
# rest of that wrapper can assume X is a clean numeric matrix and ref_idx
# is a deduplicated integer index vector. Returns the partition (ref, taxon).
check_response_inputs <- function(X, ref_idx) {
  X <- as_numeric_matrix(X, "X")

  ref_idx <- as.integer(ref_idx)

  if (length(ref_idx) < 1L) stop("ref_idx must contain at least one column index.")

  # the modeled non reference taxa: everything not in ref_idx.
  # downstream code lives on this partition (taxon_idx for response building,
  # ref_idx for denominator construction).
  taxon_idx <- setdiff(seq_len(ncol(X)), ref_idx)
  if (length(taxon_idx) < 1L) stop("No taxa left after excluding reference taxa.")

  list(X = X, ref_idx = ref_idx, taxon_idx = taxon_idx)
}

# Per sample reference mass plus summary stats.
# Called from build_reference_normalized_responses() (for the response build)
# and from select_reference_bootstrap_prefix() / pursue_main (after cleanup)
# for diagnostics on the final reference set.
# include_zero_count adds the count of samples with M_i = 0, which is the
# correctness gate for the log ratio downstream.
summarize_reference_mass <- function(X, ref_idx, include_zero_count = TRUE) {
  reference_mass <- rowSums(X[, ref_idx, drop = FALSE])   # M_i = sum_{r in R} X_ir

  out <- list(
    reference_mass = reference_mass,
    reference_mass_min = unname(min(reference_mass)),
    reference_mass_median = unname(stats::median(reference_mass)),
    reference_mass_max = unname(max(reference_mass))
  )

  if (include_zero_count) out$n_zero_reference_mass <- sum(reference_mass == 0)
  out
}


# -----------------------------------------------------------------------------
# Prevalence response construction
# -----------------------------------------------------------------------------

# Expected presence per (sample, taxon) with closed form rarefaction
# to a common depth
# Computed on the log scale via lchoose for numerical stability (choose()
# itself overflows for moderate N); the final exp() recovers a probability.
# Result clipped into (eps, 1 - eps) so downstream log transforms
# stay finite for cells where the probability is essentially 0 or 1.
compute_expected_rarefied_prevalence <- function(
    otu_taxon_counts,
    sample_depth,
    rarefaction_depth,
    eps = 1e-12
) {
  otu_taxon_counts <- as_numeric_matrix(otu_taxon_counts, "otu_taxon_counts")
  sample_depth <- as.numeric(sample_depth)

  rarefaction_depth <- as.integer(rarefaction_depth)

  if (rarefaction_depth > min(sample_depth)) {
    stop("Expected-rarefied prevalence requires rarefaction_depth <= the minimum sample depth.")
  }

  # output matrix; cells stay 0 for absent taxa (X_ij = 0 -> P_ij^ER = 0)
  prevalence_matrix <- matrix(
    0,
    nrow = nrow(otu_taxon_counts), ncol = ncol(otu_taxon_counts),
    dimnames = dimnames(otu_taxon_counts)
  )

  # loop is per sample because each row has its own N_i, but we vectorize
  # across taxa within a sample via lchoose on the X_ij vector
  for (sample_pos in seq_len(nrow(otu_taxon_counts))) {
    this_depth <- sample_depth[sample_pos]                  # N_i
    sample_counts <- otu_taxon_counts[sample_pos, ]         # X_i. (row vector)
    present_taxa <- which(sample_counts > 0)                # only taxa with X_ij > 0 need computing

    if (length(present_taxa) == 0L) next                    # nothing observed; the row of zeros stays

    log_denom <- lchoose(this_depth, rarefaction_depth)                                 # log choose(N_i, d)
    log_num   <- lchoose(this_depth - sample_counts[present_taxa], rarefaction_depth)   # log choose(N_i - X_ij, d)

    expected_presence <- 1 - exp(log_num - log_denom)                                   # P_ij^ER = 1 - choose(N_i - X_ij, d) / choose(N_i, d)
    expected_presence <- pmin(pmax(expected_presence, eps), 1 - eps)                    # clip

    prevalence_matrix[sample_pos, present_taxa] <- expected_presence
  }

  list(
    prevalence_matrix = prevalence_matrix,
    positive_count_per_taxon = colSums(otu_taxon_counts > 0),
    prevalence_per_taxon = colMeans(otu_taxon_counts > 0),
    rarefaction_depth = rarefaction_depth,
    sample_depth = sample_depth
  )
}

# Strip log_depth from a formula
remove_depth_from_prevalence_formula <- function(formula) {
  terms_obj <- stats::terms(formula)
  term_labels <- attr(terms_obj, "term.labels")
  term_labels <- setdiff(term_labels, "log_depth")

  right_side <- if (length(term_labels) == 0L) "1" else paste(term_labels, collapse = " + ")
  stats::as.formula(paste("~", right_side), env = environment(formula))
}

# Pack the prevalence response into the matrix/formula pair expected by the
# component inference engine, plus rarefaction bookkeeping.
build_prevalence_inputs_by_engine <- function(
    otu_taxon_counts,
    otu_depth_universe,
    metadata,
    formula,
    prevalence_expected_rarefied_eps = 1e-12
) {
  prevalence_formula_nondepth <- remove_depth_from_prevalence_formula(formula)

  # d = min_i N_i
  sample_depth <- rowSums(otu_depth_universe)              # N_i = sum_{j=1}^m X_ij
  rarefaction_depth <- as.integer(min(sample_depth))       # d

  if (!is.finite(rarefaction_depth) || rarefaction_depth < 1L) {
    stop("Minimum prevalence rarefaction depth must be >= 1.")
  }

  prevalence_obj <- compute_expected_rarefied_prevalence(
    otu_taxon_counts = otu_taxon_counts,
    sample_depth = sample_depth,
    rarefaction_depth = rarefaction_depth,
    eps = prevalence_expected_rarefied_eps
  )

  list(
    prevalence_matrix = prevalence_obj$prevalence_matrix,
    prevalence_formula = prevalence_formula_nondepth,
    rarefaction_depth = rarefaction_depth,
    prevalence_mode = "expected_rarefied_min_all_samples",
    sample_depth = sample_depth
  )
}


# -----------------------------------------------------------------------------
# Abundance response construction
# -----------------------------------------------------------------------------

# Abundance response
# Cell is NA when X_ij = 0 or M_i = 0.
# Called from build_reference_normalized_responses() with ref_mass = M_i
# precomputed by summarize_reference_mass() so we do not recompute it.
build_abundance_matrix <- function(X, taxon_idx, ref_idx, ref_mass) {
  taxon_counts <- X[, taxon_idx, drop = FALSE]               # X_ij over non-reference j only
  n_samples <- nrow(taxon_counts)
  n_taxa <- ncol(taxon_counts)

  # initialize as NA; only usable cells get overwritten
  abundance_matrix <- matrix(
    NA_real_,
    nrow = n_samples, ncol = n_taxa,
    dimnames = list(rownames(taxon_counts), colnames(taxon_counts))
  )

  # broadcast M_i across all non-reference taxa columns so the divide is per cell
  denominator_matrix <- matrix(
    ref_mass,
    nrow = n_samples, ncol = n_taxa,
    dimnames = list(rownames(taxon_counts), colnames(taxon_counts))
  )

  has_taxon_count   <- taxon_counts > 0        # X_ij > 0
  valid_denominator <- denominator_matrix > 0  # M_i > 0
  usable_cells <- has_taxon_count & valid_denominator   # both conditions required

  if (any(usable_cells)) {
    # A_ij = log(X_ij) - log(M_i)
    abundance_matrix[usable_cells] <-
      log(taxon_counts[usable_cells]) - log(denominator_matrix[usable_cells])
  }

  list(
    abundance_matrix = abundance_matrix,
    abundance_n_observed_per_taxon = colSums(!is.na(abundance_matrix)),
    abundance_na_count = sum(is.na(abundance_matrix)),
    abundance_inf_count = sum(is.infinite(abundance_matrix)),
    abundance_nan_count = sum(is.nan(abundance_matrix)),
    n_zero_taxon_count_cells = sum(taxon_counts == 0),
    n_nonpositive_den_cells = sum(!valid_denominator & has_taxon_count)
  )
}


# -----------------------------------------------------------------------------
# Combined response object
# -----------------------------------------------------------------------------

# Combined wrapper. Builds the log ratio abundance response matrix for the
# modeled (non reference) taxa and packs it up with reference diagnostics.
# The expected rarefied prevalence variant is built separately by pursue_main
build_reference_normalized_responses <- function(
    X,
    ref_idx,
    verbose = FALSE
) {
  input_info <- check_response_inputs(X = X, ref_idx = ref_idx)
  X <- input_info$X
  ref_idx <- input_info$ref_idx
  taxon_idx <- input_info$taxon_idx

  reference_summary <- summarize_reference_mass(X = X, ref_idx = ref_idx, include_zero_count = TRUE)

  abundance_obj <- build_abundance_matrix(
    X = X, taxon_idx = taxon_idx, ref_idx = ref_idx,
    ref_mass = reference_summary$reference_mass
  )

  out <- list(
    ref_idx = ref_idx,
    ref_names = colnames(X)[ref_idx],
    taxon_idx = taxon_idx,
    taxon_names = colnames(X)[taxon_idx],
    is_reference_taxon = seq_len(ncol(X)) %in% ref_idx,
    reference_mass = reference_summary$reference_mass,
    reference_mass_min = reference_summary$reference_mass_min,
    reference_mass_median = reference_summary$reference_mass_median,
    reference_mass_max = reference_summary$reference_mass_max,
    n_zero_reference_mass = reference_summary$n_zero_reference_mass,
    abundance_matrix               = abundance_obj$abundance_matrix,
    abundance_n_observed_per_taxon = abundance_obj$abundance_n_observed_per_taxon,
    abundance_na_count             = abundance_obj$abundance_na_count,
    abundance_inf_count            = abundance_obj$abundance_inf_count,
    abundance_nan_count            = abundance_obj$abundance_nan_count,
    n_zero_taxon_count_cells       = abundance_obj$n_zero_taxon_count_cells,
    n_nonpositive_den_cells        = abundance_obj$n_nonpositive_den_cells,
    n_samples = nrow(X),
    n_taxa_total = ncol(X),
    n_taxa_returned = length(taxon_idx)
  )

  if (verbose) {
    message("Reference-normalized response construction complete.")
    message("  Returned taxa: ", out$n_taxa_returned, " / ", out$n_taxa_total)
    message("  Reference size: ", length(out$ref_idx))
    message("  Zero reference-mass samples: ", out$n_zero_reference_mass)
    message("  Abundance NA cells: ", out$abundance_na_count)
    message("  Abundance Inf cells: ", out$abundance_inf_count)
  }

  out
}
