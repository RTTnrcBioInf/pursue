# -----------------------------------------------------------------------------
# Coercion and input helpers
# -----------------------------------------------------------------------------

# Coerce count-like inputs to a numeric matrix
as_numeric_matrix <- function(x, name = "x") {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  if (!is.matrix(x)) {
    stop(name, " must be a matrix or coercible data.frame.")
  }

  storage.mode(x) <- "double"
  x
}

# Coerce metadata-like inputs to a plain data frame
as_plain_data_frame <- function(x) {
  if (!is.data.frame(x)) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
  }

  x
}

# Accept formulas supplied as formula objects or character strings
as_model_formula <- function(formula, name = "formula") {
  if (is.character(formula)) {
    formula <- stats::as.formula(formula)
  }

  if (!inherits(formula, "formula")) {
    stop(name, " must be a formula or coercible character string.")
  }

  formula
}

# Filter taxa by explicit selection, prevalence, and total count
filter_taxa <- function(
    otu,
    taxa_keep = NULL,
    min_prevalence = NULL,
    min_total_count = NULL,
    verbose = FALSE
) {
  otu <- as_numeric_matrix(otu, "otu")
  keep_taxa <- rep(TRUE, ncol(otu))

  # Apply an explicit logical mask or taxon-index selection first.
  if (!is.null(taxa_keep)) {
    if (is.logical(taxa_keep)) {
      if (length(taxa_keep) != ncol(otu)) {
        stop("Logical taxa_keep must have length ncol(otu).")
      }

      keep_taxa <- keep_taxa & taxa_keep
    } else {
      taxa_idx <- as.integer(unique(taxa_keep))

      if (any(taxa_idx < 1L | taxa_idx > ncol(otu))) {
        stop("taxa_keep contains invalid indices.")
      }

      selected_taxa <- rep(FALSE, ncol(otu))
      selected_taxa[taxa_idx] <- TRUE
      keep_taxa <- keep_taxa & selected_taxa
    }
  }

  # Apply prevalence and library-wide abundance filters on the remaining taxa.
  if (!is.null(min_prevalence)) {
    taxa_prevalence <- colMeans(otu > 0)
    keep_taxa <- keep_taxa & (taxa_prevalence >= min_prevalence)
  }

  if (!is.null(min_total_count)) {
    keep_taxa <- keep_taxa & (colSums(otu) >= min_total_count)
  }

  kept_idx <- which(keep_taxa)
  filtered_otu <- otu[, kept_idx, drop = FALSE]

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

# Add log library size to metadata for depth adjustment
add_log_depth <- function(
    otu,
    metadata,
    depth_var_name = "log_depth"
) {
  otu <- as_numeric_matrix(otu, "otu")
  metadata <- as_plain_data_frame(metadata)

  if (nrow(otu) != nrow(metadata)) {
    stop("nrow(otu) must equal nrow(metadata).")
  }

  sample_depth <- rowSums(otu)

  if (any(!is.finite(sample_depth)) || any(sample_depth < 0)) {
    stop("Sample depths must be finite and non-negative.")
  }

  log_depth <- log(sample_depth)

  if (any(!is.finite(log_depth))) {
    stop("Non-finite log depth values encountered. Consider removing 0 depth samples.")
  }

  metadata[[depth_var_name]] <- log_depth
  metadata
}


# -----------------------------------------------------------------------------
# Formula and design-matrix helpers
# -----------------------------------------------------------------------------

# Add a covariate term to the right-hand side if missing. Used for log(depth)
add_term_if_missing <- function(formula, term) {
  formula <- as_model_formula(formula)

  if (!is.character(term) || length(term) != 1L || !nzchar(term)) {
    stop("term must be a non-empty single character string.")
  }

  right_side <- if (length(formula) == 2L) {
    formula[[2L]]
  } else {
    formula[[3L]]
  }

  rhs_vars <- all.vars(right_side)

  if (!(term %in% rhs_vars)) {
    right_side <- call("+", right_side, as.name(term))
  }

  stats::as.formula(
    call("~", right_side),
    env = environment(formula)
  )
}

# Match the user-supplied tested term to exactly one formula term
match_tested_term <- function(formula, tested_term) {
  formula <- as_model_formula(formula)

  terms_obj <- stats::terms(formula)
  term_labels <- attr(terms_obj, "term.labels")
  tested_idx <- which(term_labels == tested_term)

  if (length(tested_idx) != 1L) {
    stop("tested_term must match exactly one term label.")
  }

  list(
    terms_obj = terms_obj,
    term_labels = term_labels,
    tested_term_idx = tested_idx
  )
}

# Build the nuisance-only formula by removing the tested term
make_nuisance_formula <- function(formula, tested_term) {
  term_info <- match_tested_term(formula, tested_term)
  nuisance_terms <- setdiff(term_info$term_labels, tested_term)

  right_side <- if (length(nuisance_terms) == 0L) {
    "1"
  } else {
    paste(nuisance_terms, collapse = " + ")
  }

  stats::as.formula(paste("~", right_side))
}

# Build the nuisance-only design matrix for residualization or null fitting
make_nuisance_matrix <- function(metadata, formula, tested_term) {
  nuisance_formula <- make_nuisance_formula(formula, tested_term)
  stats::model.matrix(nuisance_formula, data = metadata)
}

# Attach the standard response column name to a model formula
make_response_formula <- function(formula, response_name = ".response") {
  formula <- as_model_formula(formula)

  right_side <- if (length(formula) == 2L) {
    formula[[2L]]
  } else {
    formula[[3L]]
  }

  stats::as.formula(
    call("~", as.name(response_name), right_side),
    env = environment(formula)
  )
}

# Map design-matrix columns to tested and nuisance terms
split_tested_and_nuisance_columns <- function(model_matrix, tested_term_idx) {
  assign_vec <- attr(model_matrix, "assign")
  tested_cols <- which(assign_vec == tested_term_idx)
  nuisance_cols <- setdiff(seq_len(ncol(model_matrix)), tested_cols)

  list(
    tested_cols = tested_cols,
    nuisance_cols = nuisance_cols
  )
}

# Check whether the tested term can support a non-degenerate model fit
tested_term_has_variation <- function(data, tested_term) {
  if (!tested_term %in% names(data)) {
    return(TRUE)
  }

  tested_values <- data[[tested_term]]

  if (is.factor(tested_values) || is.character(tested_values) || is.logical(tested_values)) {
    return(length(unique(tested_values[!is.na(tested_values)])) >= 2L)
  }

  if (is.numeric(tested_values)) {
    return(length(unique(tested_values[is.finite(tested_values)])) >= 2L)
  }

  TRUE
}


# -----------------------------------------------------------------------------
# Model-statistic and fit-result helpers
# -----------------------------------------------------------------------------

# Compute single- or multi-df Wald statistics from coefficient estimates
wald_stat <- function(coef_estimates, vcov_matrix) {
  coef_estimates <- as.numeric(coef_estimates)
  vcov_matrix <- as.matrix(vcov_matrix)

  # Use the signed Wald statistic for scalar tested terms.
  if (length(coef_estimates) == 1L) {
    se <- sqrt(vcov_matrix[1, 1])

    if (!is.finite(se) || se <= 0) {
      return(list(
        stat = 0,
        se = NA_real_,
        df = 1L,
        stat_type = "wald",
        status = "se_invalid"
      ))
    }

    return(list(
      stat = unname(coef_estimates[1] / se),
      se = unname(se),
      df = 1L,
      stat_type = "wald",
      status = "ok"
    ))
  }

  # Use the quadratic Wald statistic when the tested term spans multiple columns.
  vcov_inverse <- tryCatch(
    solve(vcov_matrix),
    error = function(e) qr.solve(vcov_matrix)
  )

  wald_value <- as.numeric(
    t(coef_estimates) %*% vcov_inverse %*% coef_estimates
  )

  list(
    stat = wald_value,
    se = NA_real_,
    df = length(coef_estimates),
    stat_type = "wald",
    status = "ok"
  )
}

# Capture model warnings without printing them during batch fitting
capture_fit_warnings <- function(expr) {
  fit_warnings <- character(0)

  value <- withCallingHandlers(
    expr,
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  list(value = value, warnings = unique(fit_warnings))
}

# Locate coefficient rows generated by scalar or factor tested terms
find_tested_coef_rows <- function(coef_names, tested_term) {
  if (!is.character(tested_term) || length(tested_term) != 1L || !nzchar(tested_term)) {
    stop("tested_term must be a non-empty single character string.")
  }

  tested_rows <- which(coef_names == tested_term)
  if (length(tested_rows) > 0L) {
    return(tested_rows)
  }

  tested_rows <- grep(paste0("^", tested_term), coef_names)
  if (length(tested_rows) > 0L) {
    return(tested_rows)
  }

  grep(paste0("^`?", tested_term, "`?"), coef_names)
}

# Extract tested-term estimates and Wald statistics from a fitted model
extract_test_result <- function(fit, tested_term) {
  fit_summary <- tryCatch(summary(fit), error = function(e) NULL)
  coef_matrix <- tryCatch(fit_summary$coefficients, error = function(e) NULL)

  if (is.null(coef_matrix)) {
    return(list(
      tested_rows = integer(0),
      estimate = NA_real_,
      se = NA_real_,
      stat = NA_real_,
      stat_type = "unavailable",
      df = NA_real_,
      status = "coef_extraction_failed"
    ))
  }

  # Match either a scalar coefficient or all dummy-coded columns for a factor term.
  coef_names <- rownames(coef_matrix)
  tested_rows <- find_tested_coef_rows(coef_names, tested_term)

  if (length(tested_rows) == 0L) {
    return(list(
      tested_rows = integer(0),
      estimate = NA_real_,
      se = NA_real_,
      stat = NA_real_,
      stat_type = "term_not_found",
      df = NA_real_,
      status = "tested_term_not_found"
    ))
  }

  vcov_matrix <- tryCatch(stats::vcov(fit), error = function(e) NULL)

  if (is.null(vcov_matrix)) {
    return(list(
      tested_rows = tested_rows,
      estimate = NA_real_,
      se = NA_real_,
      stat = NA_real_,
      stat_type = "vcov_unavailable",
      df = NA_real_,
      status = "vcov_unavailable"
    ))
  }

  # Reduce the fitted model to the tested coefficient block used by downstream p-values.
  coef_estimates <- stats::coef(fit)[tested_rows]
  tested_vcov <- as.matrix(vcov_matrix[tested_rows, tested_rows, drop = FALSE])
  wald <- wald_stat(coef_estimates = coef_estimates, vcov_matrix = tested_vcov)

  list(
    tested_rows = tested_rows,
    estimate = if (length(coef_estimates) == 1L) {
      unname(coef_estimates[1])
    } else {
      NA_real_
    },
    se = wald$se,
    stat = wald$stat,
    stat_type = wald$stat_type,
    df = if (length(coef_estimates) == 1L) {
      residual_df <- stats::df.residual(fit)
      if (is.finite(residual_df)) as.numeric(residual_df) else NA_real_
    } else {
      wald$df
    },
    status = wald$status
  )
}

# Count tested-term groups when the term is discrete or binary
tested_term_group_counts <- function(data, tested_term) {
  if (tested_term %in% names(data)) {
    tested_values <- data[[tested_term]]

    if (is.factor(tested_values) || is.character(tested_values) || is.logical(tested_values)) {
      return(table(tested_values))
    }

    if (is.numeric(tested_values)) {
      unique_values <- unique(tested_values[is.finite(tested_values)])

      if (length(unique_values) == 2L) {
        return(table(tested_values))
      }
    }
  }

  NULL
}

# Standardize fit outputs across successful and neutral taxon-level fits
make_fit_result <- function(
    n_obs,
    formula_used,
    group_counts,
    estimate = NA_real_,
    se = NA_real_,
    stat = NA_real_,
    stat_type = "unavailable",
    df = NA_real_,
    model_type = "lm",
    converged = FALSE,
    warnings = character(0),
    status = "unavailable",
    tested_rows = integer(0),
    model = NULL
) {
  list(
    estimate = estimate,
    se = se,
    stat = stat,
    stat_type = stat_type,
    df = df,
    n_obs = n_obs,
    model_type = model_type,
    converged = converged,
    warnings = warnings,
    status = status,
    tested_rows = tested_rows,
    formula_used = formula_used,
    group_counts = group_counts,
    model = model
  )
}


# -----------------------------------------------------------------------------
# Metadata alignment and permutation-group helpers
# -----------------------------------------------------------------------------

# Align metadata rows to a response matrix, using row names when available
align_metadata_to_matrix <- function(metadata, response_matrix) {
  metadata <- as_plain_data_frame(metadata)

  if (nrow(metadata) != nrow(response_matrix)) {
    stop("nrow(metadata) must equal nrow(Y).")
  }

  matrix_rows <- rownames(response_matrix)
  metadata_rows <- rownames(metadata)

  if (!is.null(matrix_rows) && !is.null(metadata_rows)) {
    if (!identical(matrix_rows, metadata_rows)) {
      if (setequal(matrix_rows, metadata_rows)) {
        metadata <- metadata[matrix_rows, , drop = FALSE]
      } else {
        stop("Row names of metadata and abundance_matrix are both present but do not match.")
      }
    }
  }

  metadata
}

# Check whether tested variables are constant within each cluster
detect_cluster_level_test <- function(metadata, tested_vars, cluster_id) {
  if (is.null(cluster_id)) {
    return(FALSE)
  }

  if (is.character(cluster_id) && length(cluster_id) == 1L) {
    if (!cluster_id %in% names(metadata)) {
      stop("cluster_id not found in metadata.")
    }

    cluster_vector <- metadata[[cluster_id]]
  } else {
    cluster_vector <- cluster_id
  }

  # A cluster-level test requires every tested variable to be constant within clusters.
  is_cluster_level <- TRUE

  for (tested_var in tested_vars) {
    if (!tested_var %in% names(metadata)) {
      stop("tested_var not found in metadata: ", tested_var)
    }

    values_by_cluster <- split(metadata[[tested_var]], cluster_vector)

    constant_within_cluster <- vapply(
      values_by_cluster,
      function(values) {
        values <- values[!is.na(values)]
        length(unique(values)) <= 1L
      },
      logical(1)
    )

    is_cluster_level <- is_cluster_level && all(constant_within_cluster)
  }

  is_cluster_level
}

# Validate or infer the metadata variables used as tested terms
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

# Require permutation grouping variables to be discrete and aligned to metadata
validate_grouping_vector <- function(x, name, n_expected) {
  if (length(x) != n_expected) {
    stop(name, " vector must have length nrow(metadata).")
  }

  keep <- !is.na(x)

  if (!any(keep)) {
    stop(name, " cannot be all NA.")
  }

  if (inherits(x, c("Date", "POSIXct", "POSIXt", "difftime"))) {
    stop(name, " must be a grouping variable with discrete labels, not a continuous time scale.")
  }

  # Numeric grouping vectors are allowed only when they behave like integer IDs.
  if (is.numeric(x)) {
    x_keep <- x[keep]

    if (any(!is.finite(x_keep))) {
      stop(name, " contains non-finite numeric values.")
    }

    if (any(abs(x_keep - round(x_keep)) > 1e-8)) {
      stop(
        name,
        " must be a grouping variable with discrete labels, not continuous numeric values. ",
        "Use factor/character labels or integer-like IDs."
      )
    }
  }

  x
}

# Extract a cluster vector supplied directly or by metadata column name
extract_cluster_vector <- function(metadata, cluster_id = NULL) {
  if (is.null(cluster_id)) {
    return(NULL)
  }

  cluster_vector <- if (is.character(cluster_id) && length(cluster_id) == 1L) {
    if (!cluster_id %in% names(metadata)) {
      stop("cluster_id not found in metadata.")
    }

    metadata[[cluster_id]]
  } else {
    cluster_id
  }

  validate_grouping_vector(
    cluster_vector,
    name = "cluster_id",
    n_expected = nrow(metadata)
  )
}

# Extract a strata vector supplied directly or by metadata column name
extract_strata_vector <- function(metadata, strata = NULL) {
  if (is.null(strata)) {
    return(NULL)
  }

  strata_vector <- if (is.character(strata) && length(strata) == 1L) {
    if (!strata %in% names(metadata)) {
      stop("strata not found in metadata.")
    }

    metadata[[strata]]
  } else {
    strata
  }

  validate_grouping_vector(
    strata_vector,
    name = "strata",
    n_expected = nrow(metadata)
  )
}
