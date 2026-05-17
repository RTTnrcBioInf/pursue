# -----------------------------------------------------------------------------
# Coercion and input helpers
# -----------------------------------------------------------------------------

# Coerce a matrix or data.frame to a numeric matrix with double storage
as_numeric_matrix <- function(x, name = "x") {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.matrix(x)) stop(name, " must be a matrix or coercible data.frame.")
  storage.mode(x) <- "double"
  x
}

# Coerce to a plain data.frame
as_plain_data_frame <- function(x) {
  if (!is.data.frame(x)) x <- as.data.frame(x, stringsAsFactors = FALSE)
  x
}

# Accept formulas as formula objects or as character strings.
# Several formula manipulators in this file then assume a real formula object,
# so this does the coercion.
as_model_formula <- function(formula, name = "formula") {
  if (is.character(formula)) formula <- stats::as.formula(formula)
  if (!inherits(formula, "formula")) stop(name, " must be a formula or coercible character string.")
  formula
}

# Filter taxa: prevalence plus total count
filter_taxa <- function(
    otu,
    min_prevalence = NULL,
    min_total_count = NULL,
    verbose = FALSE
) {
  otu <- as_numeric_matrix(otu, "otu")
  keep_taxa <- rep(TRUE, ncol(otu))

  # Filter 1: minimum proportion of samples the taxon is detected in.
  # colMeans(otu > 0) is the per taxon detection rate across samples.
  if (!is.null(min_prevalence)) {
    keep_taxa <- keep_taxa & (colMeans(otu > 0) >= min_prevalence)
  }

  # Filter 2: minimum total count summed across samples
  if (!is.null(min_total_count)) {
    keep_taxa <- keep_taxa & (colSums(otu) >= min_total_count)
  }

  kept_idx <- which(keep_taxa)
  filtered_otu <- otu[, kept_idx, drop = FALSE]

  # both kept_idx and dropped_idx are kept on the return value so callers
  # can map filtered table positions back to original columns and assemble
  # placeholder rows for the dropped taxa in the final result table
  out <- list(
    otu = filtered_otu,
    kept_idx = kept_idx,
    dropped_idx = which(!keep_taxa),
    n_taxa_input = ncol(otu),
    n_taxa_kept = ncol(filtered_otu),
    n_taxa_dropped = sum(!keep_taxa)
  )

  if (verbose) {
    message("Filtering complete.")
    message("  Input taxa: ", out$n_taxa_input)
    message("  Kept taxa: ", out$n_taxa_kept)
    message("  Dropped taxa: ", out$n_taxa_dropped)
  }

  out
}

# Stick log(library size) onto metadata as a depth covariate.
# Called from run_pursue() when depth_adjust = TRUE; the new column then
# enters the nuisance design via add_term_if_missing("log_depth").
add_log_depth <- function(otu, metadata, depth_var_name = "log_depth") {
  otu <- as_numeric_matrix(otu, "otu")
  metadata <- as_plain_data_frame(metadata)

  sample_depth <- rowSums(otu)                       # N_i
  # guard: log-depth is defined only for positive finite library sizes
  if (any(!is.finite(sample_depth)) || any(sample_depth <= 0)) {
    stop("Sample depths must be finite and positive.")
  }

  log_depth <- log(sample_depth)

  metadata[[depth_var_name]] <- log_depth
  metadata
}


# -----------------------------------------------------------------------------
# Formula and design matrix helpers
# -----------------------------------------------------------------------------

# Splice `term` onto the formula if not already on the right.
add_term_if_missing <- function(formula, term) {
  formula <- as_model_formula(formula)

  # one sided lives in [[2L]]; two sided in [[3L]]
  right_side <- if (length(formula) == 2L) formula[[2L]] else formula[[3L]]

  if (!(term %in% all.vars(right_side))) {
    right_side <- call("+", right_side, as.name(term))
  }

  stats::as.formula(call("~", right_side), env = environment(formula))
}

# Locate the tested term among the term labels of a formula
match_tested_term <- function(formula, tested_term) {
  formula <- as_model_formula(formula)

  terms_obj <- stats::terms(formula)
  term_labels <- attr(terms_obj, "term.labels")
  tested_idx <- which(term_labels == tested_term)

  # zero matches -> missing; >1 -> ambiguous label
  if (length(tested_idx) != 1L) {
    stop(
      "tested_term must match exactly one term label. ",
      "Available: ", paste(term_labels, collapse = ", ")
    )
  }

  list(terms_obj = terms_obj, term_labels = term_labels, tested_term_idx = tested_idx)
}

# Build the nuisance only formula (right side without the tested term).
# Collapses to ~ 1 when the tested term was the only right variable.
make_nuisance_formula <- function(formula, tested_term) {
  term_info <- match_tested_term(formula, tested_term)
  nuisance_terms <- setdiff(term_info$term_labels, tested_term)

  right_side <- if (length(nuisance_terms) == 0L) "1" else paste(nuisance_terms, collapse = " + ")
  stats::as.formula(paste("~", right_side))
}

# Make the nuisance only design matrix M from the nuisance formula
make_nuisance_matrix <- function(metadata, formula, tested_term) {
  stats::model.matrix(make_nuisance_formula(formula, tested_term), data = metadata)
}

# Attach a synthetic response column to the left of the formula
make_response_formula <- function(formula, response_name = ".response") {
  formula <- as_model_formula(formula)
  right_side <- if (length(formula) == 2L) formula[[2L]] else formula[[3L]]

  stats::as.formula(
    call("~", as.name(response_name), right_side),
    env = environment(formula)
  )
}

# Partition a model.matrix into tested vs nuisance column index sets via
# the `assign` attribute (which term each column originated from)
split_tested_and_nuisance_columns <- function(model_matrix, tested_term_idx) {
  assign_vec <- attr(model_matrix, "assign")
  tested_cols <- which(assign_vec == tested_term_idx)

  list(
    tested_cols = tested_cols,
    nuisance_cols = setdiff(seq_len(ncol(model_matrix)), tested_cols)
  )
}

# TRUE if the tested term has at least two distinct non NA values in data.
# Returns TRUE when tested_term is not a column (variation comes via model.matrix)
tested_term_has_variation <- function(data, tested_term) {
  if (!tested_term %in% names(data)) return(TRUE)

  v <- data[[tested_term]]

  if (is.factor(v) || is.character(v) || is.logical(v)) {
    return(length(unique(v[!is.na(v)])) >= 2L)
  }

  if (is.numeric(v)) {
    return(length(unique(v[is.finite(v)])) >= 2L)
  }

  TRUE
}


# -----------------------------------------------------------------------------
# Model statistic helpers
# -----------------------------------------------------------------------------

# Wald stat on the tested coefficient block.
# 1 df -> signed t-style theta_hat / SE; multi df -> quadratic form b' V^-1 b
wald_stat <- function(coef_estimates, vcov_matrix) {
  coef_estimates <- as.numeric(coef_estimates)
  vcov_matrix <- as.matrix(vcov_matrix)

  if (length(coef_estimates) == 1L) {
    se <- sqrt(vcov_matrix[1, 1])

    # degenerate design -> neutral stat + status flag, caller drops the taxon
    if (!is.finite(se) || se <= 0) {
      return(list(stat = 0, se = NA_real_, df = 1L, status = "se_invalid"))
    }

    return(list(
      stat = unname(coef_estimates[1] / se),
      se = unname(se),
      df = 1L, status = "ok"
    ))
  }

  # qr.solve fallback for computationally singular Sigma_hat
  vcov_inverse <- tryCatch(solve(vcov_matrix), error = function(e) qr.solve(vcov_matrix))

  wald_value <- as.numeric(t(coef_estimates) %*% vcov_inverse %*% coef_estimates)

  list(
    stat = wald_value,
    se = NA_real_,  # no single SE in the multi df case
    df = length(coef_estimates),
    status = "ok"
  )
}

# -----------------------------------------------------------------------------
# Permutation group helpers
# -----------------------------------------------------------------------------

# TRUE iff every tested variable is constant within each cluster.
# Gate for cluster level permutation: moving a whole cluster shuffles part of
# the exposure if it varies within, breaking the null.
detect_cluster_level_test <- function(metadata, tested_vars, cluster_vector) {
  for (tested_var in tested_vars) {
    values_by_cluster <- split(metadata[[tested_var]], cluster_vector)

    # NAs dropped before uniqueness check (one observed level + some NA = constant)
    constant_within_cluster <- vapply(
      values_by_cluster,
      function(v) length(unique(v[!is.na(v)])) <= 1L,
      logical(1)
    )

    if (!all(constant_within_cluster)) return(FALSE)
  }

  TRUE
}

# Resolve and validate the tested variables vector. Defaults to tested_term.
check_tested_vars <- function(metadata, tested_vars, tested_term = NULL) {
  if (is.null(tested_vars)) {
    if (is.null(tested_term) || !nzchar(tested_term)) {
      stop("Supply tested_vars or tested_term.")
    }
    tested_vars <- tested_term
  }

  if (!all(tested_vars %in% names(metadata))) {
    stop(
      "Some tested_vars are not columns in metadata: ",
      paste(setdiff(tested_vars, names(metadata)), collapse = ", ")
    )
  }

  tested_vars
}

# Validate a grouping vector: length matches metadata, at least 2 distinct levels.
validate_grouping_vector <- function(x, name, n_expected) {
  if (length(x) != n_expected) stop(name, " vector must have length nrow(metadata).")
  if (nlevels(factor(x)) < 2L) stop(name, " must have at least 2 distinct non-NA levels.")
  x
}

# Resolve cluster_id (metadata column name or vector). NULL when omitted.
extract_cluster_vector <- function(metadata, cluster_id = NULL) {
  if (is.null(cluster_id)) return(NULL)

  cluster_vector <- if (is.character(cluster_id) && length(cluster_id) == 1L) {
    if (!cluster_id %in% names(metadata)) stop("cluster_id not found in metadata.")
    metadata[[cluster_id]]
  } else cluster_id

  validate_grouping_vector(cluster_vector, name = "cluster_id", n_expected = nrow(metadata))
}

# Resolve strata (metadata column name or vector). NULL when omitted.
extract_strata_vector <- function(metadata, strata = NULL) {
  if (is.null(strata)) return(NULL)

  strata_vector <- if (is.character(strata) && length(strata) == 1L) {
    if (!strata %in% names(metadata)) stop("strata not found in metadata.")
    metadata[[strata]]
  } else strata

  validate_grouping_vector(strata_vector, name = "strata", n_expected = nrow(metadata))
}
