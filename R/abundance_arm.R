# -----------------------------------------------------------------------------
# Taxon-level abundance fitting
# -----------------------------------------------------------------------------

# Fit one reference-normalized abundance response and extract the tested-term statistic
fit_abundance_single_taxon <- function(
    y_abund,
    metadata,
    formula,
    tested_term,
    min_nonzero_n = 3L,
    return_model = FALSE,
    verbose = FALSE
) {
  if (!is.numeric(y_abund)) {
    stop("y_abund must be numeric.")
  }

  if (length(y_abund) != nrow(metadata)) {
    stop("length(y_abund) must equal nrow(metadata).")
  }


  metadata <- as_plain_data_frame(metadata)

  # Attach the taxon-specific response to the shared metadata frame.
  fit_data <- metadata
  fit_data$.response <- as.numeric(y_abund)

  # Convert the user model into the internal response-bearing formula.
  formula_full <- make_response_formula(
    formula = formula,
    response_name = ".response"
  )

  # Check that all formula variables are present before complete-case filtering.
  needed_vars <- unique(c(".response", all.vars(formula_full)))
  missing_vars <- setdiff(needed_vars, names(fit_data))

  if (length(missing_vars) > 0L) {
    stop(
      "Missing variables in metadata/formula: ",
      paste(missing_vars, collapse = ", ")
    )
  }

  # Drop rows that cannot be used by this taxon/formula pair.
  complete_rows <- stats::complete.cases(fit_data[, needed_vars, drop = FALSE])
  fit_data <- fit_data[complete_rows, , drop = FALSE]
  fit_data <- droplevels(fit_data)

  # Record support diagnostics after taxon-specific missingness is applied.
  n_obs <- nrow(fit_data)
  group_counts <- tested_term_group_counts(fit_data, tested_term)

  # Return neutral fits when the tested term collapses after complete-case filtering.
  if (!tested_term_has_variation(fit_data, tested_term)) {
    return(make_fit_result(
      n_obs = n_obs,
      formula_used = formula_full,
      group_counts = group_counts,
      estimate = 0,
      se = Inf,
      stat = 0,
      stat_type = "insufficient_tested_term_variation",
      df = 1L,
      model_type = "lm",
      converged = TRUE,
      warnings = character(0),
      status = "insufficient_tested_term_variation",
      tested_rows = integer(0),
      model = NULL
    ))
  }

  # Mark sparse responses as failed rather than manufacturing a test statistic.
  if (n_obs < min_nonzero_n) {
    return(make_fit_result(
      n_obs = n_obs,
      formula_used = formula_full,
      group_counts = group_counts,
      estimate = NA_real_,
      se = NA_real_,
      stat = NA_real_,
      stat_type = "unavailable",
      df = NA_real_,
      model_type = "lm",
      converged = FALSE,
      warnings = "too_few_nonzero_samples",
      status = "too_few_nonzero_samples",
      tested_rows = integer(0),
      model = NULL
    ))
  }

  # Treat constant responses as valid neutral tests with no usable effect variation.
  if (length(unique(fit_data$.response)) < 2L) {
    return(make_fit_result(
      n_obs = n_obs,
      formula_used = formula_full,
      group_counts = group_counts,
      estimate = 0,
      se = Inf,
      stat = 0,
      stat_type = "constant_response",
      df = 1L,
      model_type = "lm",
      converged = TRUE,
      warnings = character(0),
      status = "constant_response",
      tested_rows = integer(0),
      model = NULL
    ))
  }

  # Fit the linear model while storing warnings as diagnostics instead of printing them.
  fit_capture <- capture_fit_warnings(
    stats::lm(
      formula = formula_full,
      data = fit_data
    )
  )

  fit <- fit_capture$value
  fit_warnings <- fit_capture$warnings

  # Extract the scalar or multi-column tested-term statistic from the fitted model.
  test_result <- extract_test_result(fit, tested_term = tested_term)
  status <- if (test_result$status == "ok") "ok" else test_result$status

  # Standardize the taxon-level output used by both observed fits and null machinery.
  out <- make_fit_result(
    n_obs = n_obs,
    formula_used = formula_full,
    group_counts = group_counts,
    estimate = test_result$estimate,
    se = test_result$se,
    stat = test_result$stat,
    stat_type = test_result$stat_type,
    df = test_result$df,
    model_type = "lm",
    converged = TRUE,
    warnings = fit_warnings,
    status = status,
    tested_rows = test_result$tested_rows,
    model = if (return_model) fit else NULL
  )

  if (verbose) {
    message(
      "Abundance fit complete: status = ", out$status,
      ", stat_type = ", out$stat_type,
      ", n_obs = ", out$n_obs
    )
  }

  out
}

# Fit the abundance model independently across selected taxa
fit_abundance_arm <- function(
    abundance_matrix,
    metadata,
    formula,
    tested_term,
    taxa_subset = NULL,
    min_nonzero_n = 3L,
    return_models = FALSE,
    verbose = FALSE,
    progress_label = "Abundance arm"
) {
  abundance_matrix <- as_numeric_matrix(abundance_matrix, "abundance_matrix")
  metadata <- align_metadata_to_matrix(metadata, abundance_matrix)

  # Resolve the taxa to fit from NULL, a logical mask, or integer indices.
  if (is.null(taxa_subset)) {
    taxa_subset <- seq_len(ncol(abundance_matrix))
  } else {
    if (is.logical(taxa_subset)) {
      if (length(taxa_subset) != ncol(abundance_matrix)) {
        stop("Logical taxa_subset must have length ncol(abundance_matrix).")
      }

      taxa_subset <- which(taxa_subset)
    }

    taxa_subset <- as.integer(unique(taxa_subset))

    if (any(taxa_subset < 1L | taxa_subset > ncol(abundance_matrix))) {
      stop("taxa_subset contains invalid indices.")
    }
  }

  # Use explicit taxon names when present, otherwise generate stable display labels.
  taxa_names <- colnames(abundance_matrix)

  if (is.null(taxa_names)) {
    taxa_names <- paste0("Taxon_", seq_len(ncol(abundance_matrix)))
  }

  taxon_fits <- vector("list", length(taxa_subset))
  names(taxon_fits) <- taxa_names[taxa_subset]

  # Fit each taxon independently with shared metadata and formula settings.
  for (taxon_pos in seq_along(taxa_subset)) {
    taxon_index <- taxa_subset[taxon_pos]
    taxon_response <- abundance_matrix[, taxon_index]

    taxon_fit <- fit_abundance_single_taxon(
      y_abund = taxon_response,
      metadata = metadata,
      formula = formula,
      tested_term = tested_term,
      min_nonzero_n = min_nonzero_n,
      return_model = return_models,
      verbose = FALSE
    )

    taxon_fits[[taxon_pos]] <- taxon_fit

    if (
      verbose &&
      (
        taxon_pos %% max(1L, floor(length(taxa_subset) / 10L)) == 0L ||
        taxon_pos == length(taxa_subset)
      )
    ) {
      message(progress_label, " progress: ", taxon_pos, "/", length(taxa_subset))
    }
  }

  # Flatten taxon-level fit objects into the component result table.
  result_table <- data.frame(
    taxon_index = taxa_subset,
    taxon_name = taxa_names[taxa_subset],
    estimate = vapply(taxon_fits, function(x) x$estimate, numeric(1)),
    se = vapply(taxon_fits, function(x) x$se, numeric(1)),
    stat = vapply(taxon_fits, function(x) x$stat, numeric(1)),
    stat_type = vapply(taxon_fits, function(x) x$stat_type, character(1)),
    df = vapply(
      taxon_fits,
      function(x) {
        if (is.null(x$df) || length(x$df) == 0L || !is.finite(x$df)) {
          return(NA_real_)
        }

        as.numeric(x$df)
      },
      numeric(1)
    ),
    n_obs = vapply(taxon_fits, function(x) x$n_obs, integer(1)),
    model_type = vapply(taxon_fits, function(x) x$model_type, character(1)),
    converged = vapply(taxon_fits, function(x) x$converged, logical(1)),
    status = vapply(taxon_fits, function(x) x$status, character(1)),
    stringsAsFactors = FALSE
  )

  # Keep full models only when requested; otherwise keep lightweight summaries.
  out <- list(
    result_table = result_table,
    fits = if (return_models) taxon_fits else NULL,
    fit_summaries = if (!return_models) taxon_fits else NULL,
    tested_term = tested_term,
    formula = formula,
    min_nonzero_n = min_nonzero_n,
    n_taxa_fit = length(taxa_subset)
  )

  class(out) <- c("pursue_abundance_arm", "list")
  out
}


# -----------------------------------------------------------------------------
# Component result harmonization
# -----------------------------------------------------------------------------

# Convert prevalence or abundance fits into a shared long-format component table
component_fit_to_common_table <- function(
    fit_obj,
    component = c("auto", "prevalence", "abundance"),
    neutral_stat = 0
) {
  component <- match.arg(component)

  # Infer the component label from the class when the caller does not provide it.
  if (component == "auto") {
    if (inherits(fit_obj, "pursue_prevalence_arm")) {
      component <- "prevalence"
    } else if (inherits(fit_obj, "pursue_abundance_arm")) {
      component <- "abundance"
    } else {
      stop("Could not infer component from fit_obj class.")
    }
  }

  if (!is.list(fit_obj) || is.null(fit_obj$result_table)) {
    stop("fit_obj must contain a result_table.")
  }

  # Require the common fit columns needed by downstream union-level orchestration.
  component_table <- fit_obj$result_table
  needed_cols <- c(
    "taxon_index", "taxon_name", "estimate", "se", "stat",
    "stat_type", "df", "n_obs", "model_type", "converged", "status"
  )
  missing_cols <- setdiff(needed_cols, names(component_table))

  if (length(missing_cols) > 0L) {
    stop(
      "fit_obj$result_table is missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  # Only successful finite tests contribute their observed statistic to permutation layers.
  usable <- component_table$status == "ok" &
    is.finite(component_table$stat) &
    !is.na(component_table$estimate)

  out <- data.frame(
    taxon_index = component_table$taxon_index,
    taxon_name = component_table$taxon_name,
    component = component,
    estimate = component_table$estimate,
    se = component_table$se,
    stat = component_table$stat,
    stat_type = component_table$stat_type,
    df = component_table$df,
    n_obs = component_table$n_obs,
    model_type = component_table$model_type,
    converged = component_table$converged,
    status = component_table$status,
    usable_for_testing = usable,
    stat_for_permutation = ifelse(usable, component_table$stat, neutral_stat),
    reason_unusable = ifelse(usable, NA_character_, component_table$status),
    stringsAsFactors = FALSE
  )

  # Preserve singular-fit diagnostics when supplied by a component implementation.
  if ("singular" %in% names(component_table)) {
    out$singular <- component_table$singular
  } else {
    out$singular <- NA
  }

  out
}

# Merge prevalence and abundance component summaries into one taxon-level wide table
combine_component_tables <- function(
    prevalence_fit = NULL,
    abundance_fit = NULL,
    neutral_stat = 0,
    sort_by_taxon = TRUE
) {
  if (is.null(prevalence_fit) && is.null(abundance_fit)) {
    stop("At least one of prevalence_fit or abundance_fit must be supplied.")
  }

  # Standardize each supplied component before merging them by taxon.
  prevalence_table <- if (!is.null(prevalence_fit)) {
    component_fit_to_common_table(
      prevalence_fit,
      component = "prevalence",
      neutral_stat = neutral_stat
    )
  } else {
    NULL
  }

  abundance_table <- if (!is.null(abundance_fit)) {
    component_fit_to_common_table(
      abundance_fit,
      component = "abundance",
      neutral_stat = neutral_stat
    )
  } else {
    NULL
  }

  # Define the merged taxon universe across whatever components were fit.
  all_taxon_indices <- unique(c(
    if (!is.null(prevalence_table)) prevalence_table$taxon_index else integer(0),
    if (!is.null(abundance_table)) abundance_table$taxon_index else integer(0)
  ))

  taxon_name_map <- c(
    if (!is.null(prevalence_table)) {
      stats::setNames(prevalence_table$taxon_name, prevalence_table$taxon_index)
    } else {
      NULL
    },
    if (!is.null(abundance_table)) {
      stats::setNames(abundance_table$taxon_name, abundance_table$taxon_index)
    } else {
      NULL
    }
  )
  taxon_names <- unname(taxon_name_map[as.character(all_taxon_indices)])

  out <- data.frame(
    taxon_index = all_taxon_indices,
    taxon_name = taxon_names,
    stringsAsFactors = FALSE
  )

  # Prefix component columns before widening to avoid prevalence/abundance name collisions.
  prefix_component_columns <- function(component_table, prefix) {
    value_cols <- setdiff(
      names(component_table),
      c("taxon_index", "taxon_name", "component")
    )

    component_table <- component_table[
      ,
      c("taxon_index", "taxon_name", value_cols),
      drop = FALSE
    ]
    names(component_table)[match(value_cols, names(component_table))] <- paste0(
      prefix,
      value_cols
    )

    component_table
  }

  # Add prevalence summaries when the prevalence arm was fit.
  if (!is.null(prevalence_table)) {
    prevalence_wide <- prefix_component_columns(prevalence_table, "prev_")
    out <- merge(
      out,
      prevalence_wide,
      by = c("taxon_index", "taxon_name"),
      all.x = TRUE,
      sort = FALSE
    )
  }

  # Add abundance summaries when the abundance arm was fit.
  if (!is.null(abundance_table)) {
    abundance_wide <- prefix_component_columns(abundance_table, "abund_")
    out <- merge(
      out,
      abundance_wide,
      by = c("taxon_index", "taxon_name"),
      all.x = TRUE,
      sort = FALSE
    )
  }

  fill_missing_values <- function(df, col_name, value) {
    if (!col_name %in% names(df)) {
      return(df)
    }

    missing_rows <- is.na(df[[col_name]])
    df[[col_name]][missing_rows] <- value
    df
  }

  # Components missing for a taxon receive neutral statistics and inactive flags.
  neutral_cols <- c("prev_stat_for_permutation", "abund_stat_for_permutation")

  for (col_name in neutral_cols) {
    out <- fill_missing_values(out, col_name, neutral_stat)
  }

  logical_false_cols <- c(
    "prev_usable_for_testing", "abund_usable_for_testing",
    "prev_converged", "abund_converged"
  )

  for (col_name in logical_false_cols) {
    out <- fill_missing_values(out, col_name, FALSE)
  }

  if (sort_by_taxon) {
    out <- out[order(out$taxon_index), , drop = FALSE]
    rownames(out) <- NULL
  }

  out
}
