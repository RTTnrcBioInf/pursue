# -----------------------------------------------------------------------------
# Component inference orchestration
# -----------------------------------------------------------------------------

# Convert a component FL cache into the long table shape used by the empirical
# p-value
component_cache_to_common_table <- function(
    component_cache,
    component = c("prevalence", "abundance")
) {
  if (!inherits(component_cache, "pursue_component_permutation_cache")) {
    stop("component_cache must come from prepare_component_permutation_cache().")
  }

  component <- match.arg(component)

  data.frame(
    taxon_index = seq_len(component_cache$n_taxa),
    taxon_name = component_cache$taxon_name,
    component = component,
    estimate = component_cache$observed_estimate,
    se = component_cache$observed_se,
    stat = component_cache$observed_stat,
    df = component_cache$observed_df,
    n_obs = component_cache$n_obs,
    converged = component_cache$usable_for_testing,
    status = component_cache$status,
    usable_for_testing = component_cache$usable_for_testing,
    stat_for_permutation = component_cache$stat_for_permutation,
    reason_unusable = component_cache$reason_unusable,
    stringsAsFactors = FALSE
  )
}

# Wide format merge of prevalence + abundance cache-derived summaries.
# Each value column gets a prev_ or abund_ prefix so the two components can
# live on the same row without colliding.
combine_component_cache_tables <- function(
    prevalence_cache,
    abundance_cache,
    sort_by_taxon = TRUE
) {
  prevalence_table <- component_cache_to_common_table(
    component_cache = prevalence_cache,
    component = "prevalence"
  )
  abundance_table <- component_cache_to_common_table(
    component_cache = abundance_cache,
    component = "abundance"
  )

  if (!identical(prevalence_table$taxon_index, abundance_table$taxon_index)) {
    stop("prevalence_cache and abundance_cache must cover the same taxa in the same order.")
  }

  value_cols <- setdiff(names(prevalence_table), c("taxon_index", "taxon_name", "component"))

  prev_block <- prevalence_table[, value_cols, drop = FALSE]
  names(prev_block) <- paste0("prev_", value_cols)

  abund_block <- abundance_table[, value_cols, drop = FALSE]
  names(abund_block) <- paste0("abund_", value_cols)

  out <- cbind(
    prevalence_table[, c("taxon_index", "taxon_name"), drop = FALSE],
    prev_block,
    abund_block
  )

  if (sort_by_taxon) {
    out <- out[order(out$taxon_index), , drop = FALSE]
    rownames(out) <- NULL
  }

  out
}


# Permutation orchestration
run_component_permutations <- function(
    prevalence_matrix,
    abundance_matrix,
    metadata,
    formula,
    tested_term,
    prevalence_formula = NULL,
    tested_vars = NULL,
    n_perm = 199,
    cluster_id = NULL,
    strata = NULL,
    cluster_exposure_mode = c("sample", "cluster"),
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

  cluster_exposure_mode <- match.arg(cluster_exposure_mode)

  prevalence_formula_use <- if (is.null(prevalence_formula)) formula else prevalence_formula

  prevalence_fit_args <- list(min_nonzero_n = 1L, resid_winsorize = FALSE)

  abundance_fit_args <- modifyList(
    list(
      min_nonzero_n = 3L,
      resid_winsorize = FALSE,
      resid_winsor_quantile = 0.05,
      resid_winsor_min_n = 20L
    ),
    abundance_fit_args
  )

  perm_design <- build_permutation_design(
    metadata = metadata,
    tested_vars = tested_vars,
    tested_term = tested_term,
    cluster_id = cluster_id,
    strata = strata,
    cluster_exposure_mode = cluster_exposure_mode,
    verbose = FALSE
  )
  tested_vars <- perm_design$tested_vars

  perm_list <- generate_permutation_indices(
    perm_design = perm_design,
    n_perm = n_perm,
    seed = seed,
    verbose = FALSE
  )

  prev_cache <- prepare_component_permutation_cache(
    response_matrix = prevalence_matrix,
    metadata = metadata,
    formula = prevalence_formula_use,
    tested_term = tested_term,
    tested_vars = tested_vars,
    fit_args = prevalence_fit_args,
    neutral_stat = neutral_stat,
    perm_design = perm_design,
    verbose = verbose,
    component = "prevalence",
    cache_label = "Prevalence"
  )

  abund_cache <- prepare_component_permutation_cache(
    response_matrix = abundance_matrix,
    metadata = metadata,
    formula = formula,
    tested_term = tested_term,
    tested_vars = tested_vars,
    fit_args = abundance_fit_args,
    neutral_stat = neutral_stat,
    perm_design = perm_design,
    verbose = verbose,
    component = "abundance",
    cache_label = "Abundance"
  )

  observed_component_table <- combine_component_cache_tables(
    prevalence_cache = prev_cache,
    abundance_cache = abund_cache
  )

  prev_null <- compute_component_permutation_null(
    component_cache = prev_cache,
    perm_list = perm_list,
    perm_design = perm_design,
    neutral_stat = neutral_stat,
    parallelize_permutations = parallelize_permutations,
    permutation_n_cores = permutation_n_cores,
    verbose = verbose,
    null_label = "Prevalence"
  )

  if (verbose) message(strrep("-", 50))

  abund_null <- compute_component_permutation_null(
    component_cache = abund_cache,
    perm_list = perm_list,
    perm_design = perm_design,
    neutral_stat = neutral_stat,
    parallelize_permutations = parallelize_permutations,
    permutation_n_cores = permutation_n_cores,
    verbose = verbose,
    null_label = "Abundance"
  )

  if (verbose) message(strrep("-", 50))

  prev_perm_stats <- prev_null$stat_mat
  abund_perm_stats <- abund_null$stat_mat

  colnames(prev_perm_stats) <- paste0("perm_", seq_len(n_perm))
  colnames(abund_perm_stats) <- paste0("perm_", seq_len(n_perm))
  rownames(prev_perm_stats) <- observed_component_table$taxon_name
  rownames(abund_perm_stats) <- observed_component_table$taxon_name

  out <- list(
    observed_component_table = observed_component_table,
    observed_prev_stat = prev_cache$stat_for_permutation,
    observed_abund_stat = abund_cache$stat_for_permutation,
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
    abund_null = abund_null,
    prev_perm_status = prev_null$status_mat,
    abund_perm_status = abund_null$status_mat
  )
  class(out) <- c("pursue_component_permutations", "list")

  if (verbose) {
    message("Component resampling complete.")
    message("  Taxa: ", nrow(observed_component_table))
    message("  Permutations: ", n_perm)
    message("  Prev usable: ", sum(observed_component_table$prev_usable_for_testing, na.rm = TRUE))
    message("  Abund usable: ", sum(observed_component_table$abund_usable_for_testing, na.rm = TRUE))
    message("  Permutation unit level: ", perm_design$unit_level)
  }

  out
}

# -----------------------------------------------------------------------------
# Empirical p-values and union layer
# -----------------------------------------------------------------------------

# Empirical p per taxon, per component.
# Abs value |W|. only perm slots with status "ok" count toward the
# denominator.
# Verbose path shows the failed slot fraction across components and warns
# about taxa whose valid perm count fell below min_valid_perm_warn.
compute_empirical_component_pvalues <- function(
    perm_obj,
    min_valid_perm_warn = 200L,
    verbose = FALSE
) {
  if (!inherits(perm_obj, "pursue_component_permutations")) {
    stop("perm_obj must come from run_component_permutations().")
  }

  observed_table   <- perm_obj$observed_component_table
  n_taxa           <- nrow(observed_table)
  n_perm           <- perm_obj$n_perm

  prev_stat_mat    <- perm_obj$prev_perm_stats
  abund_stat_mat   <- perm_obj$abund_perm_stats
  prev_status_mat  <- perm_obj$prev_perm_status
  abund_status_mat <- perm_obj$abund_perm_status

  prev_p_values    <- rep(NA_real_, n_taxa)
  abund_p_values   <- rep(NA_real_, n_taxa)
  prev_n_valid     <- rep(NA_integer_, n_taxa)
  abund_n_valid    <- rep(NA_integer_, n_taxa)

  # Per taxon loop. The two components are computed independently in the same
  # loop body rather than factored into a helper.
  for (taxon_pos in seq_len(n_taxa)) {
    observed_prev  <- perm_obj$observed_prev_stat[taxon_pos]
    observed_abund <- perm_obj$observed_abund_stat[taxon_pos]

    # prevalence component: keep only perm slots whose status is "ok"
    prev_ok          <- prev_status_mat[taxon_pos, ] == "ok"
    prev_null_valid  <- prev_stat_mat[taxon_pos, prev_ok]
    n_valid_prev     <- length(prev_null_valid)

    prev_n_valid[taxon_pos]  <- n_valid_prev
    # p_jc = (1 + sum_b indicator(|W_jc_b| >= |W_jc_obs|)) / (1 + B_jc)
    prev_p_values[taxon_pos] <- if (n_valid_prev == 0L) {
      1                          # no valid perms; cannot reject so treat as p = 1
    } else {
      (1 + sum(prev_null_valid >= observed_prev, na.rm = TRUE)) / (1 + n_valid_prev)
    }

    # abundance component: identical machinery on the abundance null matrix
    abund_ok         <- abund_status_mat[taxon_pos, ] == "ok"
    abund_null_valid <- abund_stat_mat[taxon_pos, abund_ok]
    n_valid_abund    <- length(abund_null_valid)

    abund_n_valid[taxon_pos]  <- n_valid_abund
    abund_p_values[taxon_pos] <- if (n_valid_abund == 0L) {
      1
    } else {
      (1 + sum(abund_null_valid >= observed_abund, na.rm = TRUE)) / (1 + n_valid_abund)
    }
  }

  # Unusable observed taxa still got a null matrix computed for them (the
  # cache laid out a neutral_stat row, so the per perm |W| stays at neutral_stat).
  prev_p_values[!observed_table$prev_usable_for_testing]   <- 1
  abund_p_values[!observed_table$abund_usable_for_testing] <- 1

  if (verbose) {
    prev_n_failed  <- n_perm - prev_n_valid
    abund_n_failed <- n_perm - abund_n_valid

    total_prev_failed  <- sum(prev_n_failed,  na.rm = TRUE)
    total_abund_failed <- sum(abund_n_failed, na.rm = TRUE)
    total_cells        <- n_taxa * n_perm

    message("Empirical p-value computation complete.")
    message(
      "  Prevalence component:  ", total_prev_failed, " / ", total_cells,
      " permutation slots (perm x feature) failed or were neutralized (",
      round(100 * total_prev_failed / total_cells, 1), "%)"
    )
    message(
      "  Abundance component:   ", total_abund_failed, " / ", total_cells,
      " permutation slots (perm x feature) failed or were neutralized (",
      round(100 * total_abund_failed / total_cells, 1), "%)"
    )

    prev_low  <- which(prev_n_valid  < min_valid_perm_warn &
                         observed_table$prev_usable_for_testing)
    abund_low <- which(abund_n_valid < min_valid_perm_warn &
                         observed_table$abund_usable_for_testing)

    if (length(prev_low) > 0L) {
      message("  WARNING: ", length(prev_low),
              " taxon/taxa had fewer than ", min_valid_perm_warn,
              " valid prevalence permutations:")
      for (i in prev_low) {
        message("    ", observed_table$taxon_name[i],
                " (", prev_n_valid[i], " / ", n_perm, " valid)")
      }
    }

    if (length(abund_low) > 0L) {
      message("  WARNING: ", length(abund_low),
              " taxon/taxa had fewer than ", min_valid_perm_warn,
              " valid abundance permutations:")
      for (i in abund_low) {
        message("    ", observed_table$taxon_name[i],
                " (", abund_n_valid[i], " / ", n_perm, " valid)")
      }
    }

    if (total_prev_failed == 0L && total_abund_failed == 0L) {
      message("  No failed permutation slots. All denominators are ", n_perm + 1L, ".")
    }
  }

  out <- observed_table
  out$prev_empirical_p   <- prev_p_values
  out$abund_empirical_p  <- abund_p_values
  out$prev_n_valid_perm  <- prev_n_valid
  out$abund_n_valid_perm <- abund_n_valid
  out
}


# ACAT combine the per-component p values to get one union p value per taxon
combine_observed_components <- function(
    component_p_table,
    acat_weights = c(0.5, 0.5)
) {
  # require the columns the ACAT step depends on (per-component p values plus
  # per-component usable flags).
  needed <- c(
    "taxon_index", "taxon_name",
    "prev_empirical_p", "abund_empirical_p",
    "prev_usable_for_testing", "abund_usable_for_testing"
  )
  missing_cols <- setdiff(needed, names(component_p_table))
  if (length(missing_cols) > 0L) {
    stop("component_p_table is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  if (
    !is.numeric(acat_weights) || length(acat_weights) != 2L ||
      any(!is.finite(acat_weights)) || any(acat_weights < 0)
  ) {
    stop("acat_weights must be two finite non-negative numeric values.")
  }

  n_taxa <- nrow(component_p_table)
  union_p <- rep(NA_real_, n_taxa)

  # Per taxon loop. Each taxon is combined independently, so we just look
  # up the two p values and the two usable flags
  for (i in seq_len(n_taxa)) {
    p_vec  <- c(component_p_table$prev_empirical_p[i],
                component_p_table$abund_empirical_p[i])
    active <- c(component_p_table$prev_usable_for_testing[i],
                component_p_table$abund_usable_for_testing[i])

    # both components unusable: there is no signal to combine, return 1
    if (!any(active)) {
      union_p[i] <- 1
      next
    }

    # drop unusable components and redistribute the weight over the rest
    w <- acat_weights[active]

    if (sum(w) <= 0) w <- rep(1, sum(active))

    w <- w / sum(w)
    union_p[i] <- ACAT::ACAT(p_vec[active], weights = w)
  }

  out <- component_p_table
  out$union_p_value <- union_p
  out
}


# -----------------------------------------------------------------------------
# Independent filtering and final union layer
# -----------------------------------------------------------------------------

#  ACAT combine the per-component empirical p values,
# BH adjust, optionally apply independent filtering, and emit the per taxon
# rejection flags at the defined q_alpha.
# Independent filtering drops the bottom independent_filter_quantile of
# taxa before BH
compute_union_layer <- function(
    component_p_table,
    acat_weights = c(0.5, 0.5),
    independent_filtering = FALSE,
    q_alpha = 0.05,
    independent_filter_quantile = 0.20,
    independent_filter_stat = NULL,
    independent_filter_filtered_q_value = NA_real_
) {
  # ACAT combine the two components into p_union per taxon
  union_table <- combine_observed_components(
    component_p_table = component_p_table,
    acat_weights = acat_weights
  )

  # alias kept for backwards compatibility; some plotting code expects this
  union_table$union_combined_stat_obs <- union_table$union_p_value

  # Raw BH over ALL modeled taxa (the layer without independent filtering)
  union_q_raw <- stats::p.adjust(union_table$union_p_value, method = "BH")

  union_table$union_q_value_raw_bh <- union_q_raw
  union_table$if_enabled           <- isTRUE(independent_filtering)
  # IF diagnostic columns: initialized empty here, overwritten below when IF runs
  union_table$if_filter_stat       <- NA_real_
  union_table$if_keep              <- NA
  union_table$if_threshold         <- NA_real_

  union_q_final <- union_q_raw

  if (isTRUE(independent_filtering)) {

    filter_stat <- as.numeric(independent_filter_stat)

    # threshold = independent_filter_quantile of F_j
    cutoff <- as.numeric(stats::quantile(
      filter_stat[is.finite(filter_stat)],
      probs  = independent_filter_quantile,
      names  = FALSE,
      na.rm  = TRUE
    ))

    # taxa at or above the threshold survive to the final BH set
    keep <- is.finite(filter_stat) & (filter_stat >= cutoff)

    # Filtered out taxa get a placeholder q value (NA by default)
    union_q_final <- rep(independent_filter_filtered_q_value, nrow(union_table))
    if (any(keep)) {
      # BH applied only to the kept taxa
      union_q_final[keep] <- stats::p.adjust(union_table$union_p_value[keep], method = "BH")
    }

    union_table$if_filter_stat <- filter_stat
    union_table$if_keep        <- keep
    union_table$if_threshold   <- cutoff
  }

  union_table$union_q_value <- union_q_final
  union_table$q_alpha       <- q_alpha

  # Two rejection columns are returned:
  # rejected_bh_raw : BH over all modeled taxa
  # rejected_bh     : (raw BH without IF, IF restricted BH with IF)

  union_table$rejected_bh_raw <- is.finite(union_q_raw)   & !is.na(union_q_raw)   & union_q_raw   <= q_alpha
  union_table$rejected_bh     <- is.finite(union_q_final) & !is.na(union_q_final) & union_q_final <= q_alpha

  list(
    result_table                = union_table,
    independent_filtering       = isTRUE(independent_filtering),
    q_alpha                     = q_alpha,
    independent_filter_quantile = independent_filter_quantile
  )
}
