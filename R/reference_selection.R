# -----------------------------------------------------------------------------
# Reference selector orchestrator
# -----------------------------------------------------------------------------

# Pipeline: score -> bootstrap rerank -> smallest stable prefix
#           -> optional hostage-prevention expansion.
select_reference_bootstrap_prefix <- function(
    X,
    metadata,
    formula,
    tested_term,
    p_max = 1000L,
    pseudocount = 1,
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
    return_ratio_matrix = FALSE,
    seed = 1,
    verbose = FALSE
) {
  # Caller: candidate scoring on the full data
  score_obj <- score_reference_candidates(
    X = X,
    metadata = metadata,
    formula = formula,
    tested_term = tested_term,
    p_max = p_max,
    pseudocount = pseudocount,
    return_ratio_matrix = return_ratio_matrix
  )

  # Caller: bootstrap rerank
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

  # Caller: bootstrap inclusion rates per prefix size
  stability_table <- compute_prefix_stability(
    full_sorted_candidate_taxa = score_obj$sorted_candidate_taxa,
    bootstrap_rankings = bootstrap_rankings
  )

  # Caller: reference-mass adequacy per prefix size
  adequacy_table <- compute_prefix_adequacy(
    X = X,
    full_sorted_candidate_taxa = score_obj$sorted_candidate_taxa,
    bulk_count_quantile = bulk_count_quantile,
    bulk_rel_quantile = bulk_rel_quantile
  )

  # Caller: smallest prefix that clears stability_threshold
  selection_obj <- find_stable_prefix(
    prefix_stability = stability_table,
    prefix_adequacy = adequacy_table,
    stability_threshold = bootstrap_stability_threshold
  )

  selection_table <- selection_obj$selection_table

  if (hostage_prevention) {
    # Caller: hostage-prevention expansion (only when mass conditions fail at k0)
    hostage_obj <- extend_fragile_reference_prefix(
      selection_table = selection_table,
      sorted_scores = score_obj$sorted_scores,
      n_samples = nrow(X),
      k0 = selection_obj$selected_k,
      X = X,
      full_sorted_candidate_taxa = score_obj$sorted_candidate_taxa,
      require_no_zero_mass = require_no_zero_mass,
      min_bulk_count = min_bulk_count,
      min_bulk_rel = min_bulk_rel,
      hostage_max_size_inflation = hostage_max_size_inflation,
      hostage_max_score_inflation = hostage_max_score_inflation,
      hostage_max_stability_drop = hostage_max_stability_drop
    )

    # Warn when the function had to exceed the expansion budget
    # specifically to rescue samples with zero reference mass.
    if (isTRUE(hostage_obj$override_invoked)) {
      warning(
        "Reference expansion exceeded collateral budget to satisfy zero-mass requirement. ",
        "Final prefix size: ", hostage_obj$selected_k,
        " (base: ", selection_obj$selected_k,
        ", size inflation: ", round(100 * hostage_obj$final_size_inflation, 1), "%",
        ", score inflation: ", round(100 * hostage_obj$final_score_inflation, 1), "%",
        ", stability drop: ", round(100 * hostage_obj$final_stability_drop, 1), "%). ",
        "Sample(s) requiring expansion: ",
        paste(hostage_obj$override_samples, collapse = ", "),
        ". Consider removing these samples or inspecting their library composition.",
        call. = FALSE
      )
    }

    # Warn when the function had to exceed the expansion budget but still failed to satisfy all adequacy conditions
    if (identical(hostage_obj$hostage_reason, "budget_exhausted_quality_compromise")) {
      failed_conds <- character(0)
      if (!hostage_obj$final_cond2_ok) failed_conds <- c(failed_conds, "condition 2 (bulk count quantile)")
      if (!hostage_obj$final_cond3_ok) failed_conds <- c(failed_conds, "condition 3 (bulk relative mass quantile)")

      compromise_sample_text <- format_reference_mass_problem_samples(hostage_obj$compromise_samples)

      # Warn when expansion fails to rescue conditions 2 and 3 within expansion budget
      warning(
        "Reference expansion exhausted collateral budget without satisfying all adequacy ",
        "conditions. Failed: ", paste(failed_conds, collapse = "; "),
        ". Final prefix size: ", hostage_obj$selected_k,
        if (nzchar(compromise_sample_text)) {
          paste0(". Samples still below failed mass thresholds: ", compromise_sample_text)
        } else "",
        ". Log-ratios will be defined but may have elevated noise.",
        call. = FALSE
      )
    }

    selected_k <- hostage_obj$selected_k
  } else {
    hostage_obj <- list(   # if hostage prevention is off, fill in a placeholder object
      selected_k = selection_obj$selected_k,
      hostage_base_k = selection_obj$selected_k,
      hostage_reason = "not_used",
      override_invoked = FALSE,
      override_samples = character(0),
      expansion_trigger_samples = make_empty_reference_mass_problem_table(),
      compromise_samples = make_empty_reference_mass_problem_table(),
      final_cond1_ok = NA,
      final_cond2_ok = NA,
      final_cond3_ok = NA
    )
    selected_k <- selection_obj$selected_k # return selection without mass checks
  }

  selected_references <- score_obj$sorted_candidate_taxa[seq_len(selected_k)] # final selected reference indices
  selected_reference_names <- colnames(X)[selected_references]                # before group cleanup

  mass_diag <- summarize_reference_mass(
    X = X, ref_idx = selected_references, include_zero_count = FALSE
  )

  # collect all refsel outputs + diagnostics
  out <- list(
    X = X,
    selected_references = selected_references,
    selected_reference_names = selected_reference_names,
    selected_cutpoint = selected_k,
    selected_MinAbundance = unname(min(rowSums(X[, selected_references, drop = FALSE]))),

    scores = score_obj$scores_full,
    sorted_candidate_taxa = score_obj$sorted_candidate_taxa,
    sorted_candidate_names = score_obj$sorted_candidate_names,
    sorted_scores = score_obj$sorted_scores,

    ratio_matrix = score_obj$ratio_matrix,
    reference_mass = mass_diag$reference_mass,
    reference_mass_min = mass_diag$reference_mass_min,
    reference_mass_median = mass_diag$reference_mass_median,
    reference_mass_max = mass_diag$reference_mass_max,

    bootstrap_B = bootstrap_B,
    bootstrap_stability_threshold = bootstrap_stability_threshold,
    require_no_zero_mass = require_no_zero_mass,
    min_bulk_count = min_bulk_count,
    bulk_count_quantile = bulk_count_quantile,
    min_bulk_rel = min_bulk_rel,
    bulk_rel_quantile = bulk_rel_quantile,
    bootstrap_hostage_prevention = hostage_prevention,
    bootstrap_override_invoked = hostage_obj$override_invoked,
    bootstrap_override_samples = hostage_obj$override_samples,
    bootstrap_expansion_trigger_samples = hostage_obj$expansion_trigger_samples,
    bootstrap_compromise_samples = hostage_obj$compromise_samples,
    bootstrap_final_cond1_ok = hostage_obj$final_cond1_ok,
    bootstrap_final_cond2_ok = hostage_obj$final_cond2_ok,
    bootstrap_final_cond3_ok = hostage_obj$final_cond3_ok,
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

      trigger_sample_text <- format_reference_mass_problem_samples(hostage_obj$expansion_trigger_samples)
      if (nzchar(trigger_sample_text)) {
        message("  Reference expansion trigger samples: ", trigger_sample_text)
      }
    }
  }

  out
}

# -----------------------------------------------------------------------------
# Function definitions
# -----------------------------------------------------------------------------

# all taxa are used except when p_max is defined
select_reference_candidates <- function(X, p_max = 1000L) {
  candidate_idx <- seq_len(ncol(X))

  if (length(candidate_idx) > p_max) {
    abundance_order <- order(colSums(X[, candidate_idx, drop = FALSE]), decreasing = TRUE) # when pmax is on,
    candidate_idx <- candidate_idx[abundance_order[seq_len(p_max)]]                        # select by total abundance
  }

  candidate_idx
}

# Score candidates by nuisance-adjusted pairwise log-ratio stability:
# pairwise log-ratios per (j, k), OLS-residualize against the nuisance
# design, score each taxon by the median residual SD across its partners.
score_reference_candidates <- function(
    X,
    metadata,
    formula,
    tested_term,
    p_max = 1000L,
    pseudocount = 1,
    return_ratio_matrix = TRUE
) {
  nuisance_matrix <- make_nuisance_matrix(                # M : nuisance design
    metadata = metadata, formula = formula, tested_term = tested_term
  )

  candidate_idx <- select_reference_candidates(X, p_max)
  candidate_counts <- X[, candidate_idx, drop = FALSE]

  candidate_table <- t(candidate_counts) + pseudocount   # X_ij + c
  log_candidate_counts <- log(candidate_table)           # log(X_ij + c)

  n_candidates <- nrow(log_candidate_counts)
  n_samples <- ncol(log_candidate_counts)

  # (M^T M)^-1, with qr.solve as a fallback when M is singular (determinant = 0) (safety guard)
  xtx_inverse <- tryCatch(
    solve(crossprod(nuisance_matrix)),
    error = function(e) qr.solve(crossprod(nuisance_matrix))
  )

  projection_solver <- xtx_inverse %*% t(nuisance_matrix)  # (M^T M)^-1 M^T : applied to y_jk gives beta_hat_jk
  residual_df <- n_samples - ncol(nuisance_matrix)         # n - q

  # number of nuisance variables should not exceed sample size
  if (residual_df <= 0L) {
    stop("Not enough samples to residualize pairwise log-ratios against nuisance design.")
  }

  # first initialize an empty matrix to hold pairwise residual SDs
  pairwise_sd <- matrix(
    NA_real_,
    nrow = n_candidates, ncol = n_candidates,
    dimnames = list(colnames(candidate_counts), colnames(candidate_counts))
  )

  # iterate through candidates, with -1 because we compare each candidate to those ranked after it
  for (candidate_pos in seq_len(n_candidates - 1L)) {
    comparison_block <- log_candidate_counts[(candidate_pos + 1L):n_candidates, , drop = FALSE]

    # y_i_jk = log(X_ij + c) - log(X_ik + c) for every k ranked after j
    log_ratio_block <- matrix(
      log_candidate_counts[candidate_pos, ],
      nrow = nrow(comparison_block),
      ncol = ncol(comparison_block),
      byrow = TRUE
    ) - comparison_block

    log_ratio_block <- t(log_ratio_block)                                       # samples x pairs

    coef_estimates <- projection_solver %*% log_ratio_block                     # beta_hat_jk = (M^T M)^-1 M^T y_jk
    residual_matrix <- log_ratio_block - nuisance_matrix %*% coef_estimates     # r_jk = y_jk - M beta_hat_jk
    residual_sd <- sqrt(colSums(residual_matrix^2) / residual_df)               # sigma_jk = sqrt(sum_i r_i_jk^2 / (n - q))

    pairwise_sd[candidate_pos, (candidate_pos + 1L):n_candidates] <- residual_sd # fill in the upper triangle of the pairwise SD matrix
    pairwise_sd[(candidate_pos + 1L):n_candidates, candidate_pos] <- residual_sd # mirror to the lower triangle
  }                                                                              # mirror is needed for apply() so that each taxon has the same amount of rows
  # the median is equivalent to unique pairs, because each vector entry is essentially doubled
  diag(pairwise_sd) <- NA_real_ # ignore self comparison

  # per taxon median sigma_jk across pair partners
  median_pairwise_score <- apply(pairwise_sd, 1L, median, na.rm = TRUE) # apply median over rows (taxa)
  score_order <- order(median_pairwise_score, na.last = TRUE)

  sorted_candidate_taxa <- candidate_idx[score_order]
  sorted_candidate_names <- colnames(X)[sorted_candidate_taxa]
  sorted_scores <- median_pairwise_score[score_order]

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

# Rerank candidates on bootstrap resamples for downstream stability scoring
# Rows of (X, metadata) are resampled jointly
bootstrap_reference_rankings <- function(
    X,
    metadata,
    formula,
    tested_term,
    p_max = 1000L,
    pseudocount = 1,
    n_bootstrap = 100L,
    seed = 1
) {
  if (!is.null(seed)) set.seed(seed)

  n_samples <- nrow(X)
  bootstrap_rankings <- vector("list", n_bootstrap)

  for (boot_pos in seq_len(n_bootstrap)) {
    sample_rows <- sample.int(n_samples, size = n_samples, replace = TRUE) # get sampled sample indices

    boot_counts <- X[sample_rows, , drop = FALSE]  # boot otu table
    boot_metadata <- metadata[sample_rows, , drop = FALSE] # boot metadata

    # unique row names per boot (model.matrix breaks on duplicates)
    rownames(boot_counts) <- paste0("boot_", boot_pos, "_", seq_len(nrow(boot_counts)))
    rownames(boot_metadata) <- rownames(boot_counts)

    boot_score <- score_reference_candidates( # repeat scoring with new otu and metadata tables
      X = boot_counts, metadata = boot_metadata,
      formula = formula, tested_term = tested_term,
      p_max = p_max, pseudocount = pseudocount,
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
# Prefix feasibility
# -----------------------------------------------------------------------------

# Bootstrap stability per prefix size
# Computed as the mean over taxa in R_k of the per taxon inclusion rate
# across bootstraps (algebraically equal to the per bootstrap overlap form).
compute_prefix_stability <- function(
    full_sorted_candidate_taxa,
    bootstrap_rankings
) {
  n_candidates <- length(full_sorted_candidate_taxa)
  n_bootstrap  <- length(bootstrap_rankings)
  stability    <- numeric(n_candidates)

  for (prefix_size in seq_len(n_candidates)) {                       # k
    full_prefix <- full_sorted_candidate_taxa[seq_len(prefix_size)]  # R_k
    inclusion_rate <- numeric(prefix_size)

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

    stability[prefix_size] <- mean(inclusion_rate)  # S_k = (1/B) sum_b O_k^(b)
  }

  data.frame(
    k = seq_len(n_candidates),
    stability = stability,
    n_boot = n_bootstrap,
    stringsAsFactors = FALSE
  )
}

# Per prefix size: reference mass diagnostics used by the three adequacy
# conditions
compute_prefix_adequacy <- function(
    X,
    full_sorted_candidate_taxa,
    bulk_count_quantile = 0.05,
    bulk_rel_quantile = 0.10
) {
  n_candidates <- length(full_sorted_candidate_taxa)

  sample_depth <- rowSums(X)                        # sum_{j=1}^m X_ij = N_i total library size
  adequacy_rows <- vector("list", n_candidates)

  for (prefix_size in seq_len(n_candidates)) {      # k
    reference_indices <- full_sorted_candidate_taxa[seq_len(prefix_size)]

    reference_mass <- rowSums(X[, reference_indices, drop = FALSE])   # sum_{r in R_k} X_ir = M_i(k)   sample-wise reference mass
    relative_reference_mass <- reference_mass / pmax(sample_depth, 1) # M_i(k) / N_i # relative reference mass

    adequacy_rows[[prefix_size]] <- data.frame(
      k = prefix_size,
      n_zero_ref_mass = sum(reference_mass == 0),   # cond. 1 input
      bulk_count_q = as.numeric(stats::quantile(reference_mass, probs = bulk_count_quantile, names = FALSE, na.rm = TRUE)),  # cond. 2 input
      bulk_rel_q   = as.numeric(stats::quantile(relative_reference_mass, probs = bulk_rel_quantile, names = FALSE, na.rm = TRUE)),  # cond. 3 input
      ref_mass_min = min(reference_mass),
      ref_mass_median = unname(stats::median(reference_mass)),
      ref_mass_max = max(reference_mass),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, adequacy_rows)
}

# Evaluate the three adequacy conditions for one prefix row
evaluate_adequacy_conditions <- function(
    adequacy_row,
    require_no_zero_mass,
    min_bulk_count,
    min_bulk_rel
) {
  # cond. 1 : no sample has M_i(k) = 0
  cond1_ok <- if (isTRUE(require_no_zero_mass)) adequacy_row$n_zero_ref_mass == 0L else TRUE
  # cond. 2 : 5th percentile of M_i(k)       >= min_bulk_count
  cond2_ok <- is.finite(adequacy_row$bulk_count_q) && (adequacy_row$bulk_count_q >= min_bulk_count)
  # cond. 3 : 5th percentile of M_i(k) / N_i >= min_bulk_rel
  cond3_ok <- is.finite(adequacy_row$bulk_rel_q)   && (adequacy_row$bulk_rel_q   >= min_bulk_rel)

  list(
    cond1_zero_mass_ok = cond1_ok,
    cond2_bulk_count_ok = cond2_ok,
    cond3_bulk_rel_ok = cond3_ok,
    all_ok = cond1_ok && cond2_ok && cond3_ok
  )
}

# Smallest prefix that clears the stability threshold. Mass adequacy is
# handled separately downstream by extend_fragile_reference_prefix()
# Called from the orchestrator select_reference_bootstrap_prefix() once the
# stability + adequacy tables are built
find_stable_prefix <- function(
    prefix_stability,
    prefix_adequacy,
    stability_threshold = 0.80
) {
  # Join the two diagnostic tables by prefix size; sort = TRUE keeps k ascending
  selection_table <- merge(prefix_stability, prefix_adequacy, by = "k", sort = TRUE)

  # Feasible = S_k >= S_min (NaN stabilities are not feasible)
  feasible <- is.finite(selection_table$stability) &
    (selection_table$stability >= stability_threshold)

  feasible_rows <- which(feasible)
  if (length(feasible_rows) == 0L) {
    stop("No prefix satisfied the bootstrap stability threshold.")
  }

  selected_row <- feasible_rows[1L]   # smallest k that clears the threshold

  list(
    selected_k = selection_table$k[selected_row],
    selection_table = selection_table,     # full table, extend_fragile_reference_prefix uses it for k_prime > k_star rows
    selected_row = selection_table[selected_row, , drop = FALSE]
  )
}

# constructor for the empty problem sample table (zero row template).
# used as the no-problem-samples return from diagnose_prefix_mass_problem_samples()
# and as the placeholder for `expansion_trigger_samples` and `compromise_samples
# in extend_fragile_reference_prefix() when those branches dont fire.
make_empty_reference_mass_problem_table <- function() {
  data.frame(
    sample_index = integer(0),
    sample_name = character(0),
    k = integer(0),
    reason = character(0),
    reference_mass = numeric(0),
    sample_depth = numeric(0),
    relative_reference_mass = numeric(0),
    threshold = numeric(0),
    stringsAsFactors = FALSE
  )
}

# Per sample reasons for failed adequacy at prefix k
# Quantile conditions report the failing low mass tail (not a unique
# causal subset, its the region, not the cause).
# Called twice from extend_fragile_reference_prefix() once for
# expansion_trigger_samples (at k0) and once for compromise_samples (at the
# final selected_k).
diagnose_prefix_mass_problem_samples <- function(
    X,
    full_sorted_candidate_taxa,
    k,
    conditions,
    require_no_zero_mass = TRUE,
    min_bulk_count = 10,
    min_bulk_rel = 0.01
) {
  reference_indices <- full_sorted_candidate_taxa[seq_len(k)]

  reference_mass <- rowSums(X[, reference_indices, drop = FALSE])      # M_i(k)
  sample_depth <- rowSums(X)                                            # N_i
  relative_reference_mass <- reference_mass / pmax(sample_depth, 1)     # M_i(k) / N_i  (pmax avoids 0/0)

  sample_names <- rownames(X)

  # one row per (failed condition, sample). NULL when the
  # failed sample set is empty so we dont pollute the rbind with zero row frames
  build_rows <- function(sample_pos, reason, threshold) {
    if (length(sample_pos) == 0L) return(NULL)
    data.frame(
      sample_index = sample_pos,
      sample_name = sample_names[sample_pos],
      k = k,
      reason = reason,
      reference_mass = unname(reference_mass[sample_pos]),
      sample_depth = unname(sample_depth[sample_pos]),
      relative_reference_mass = unname(relative_reference_mass[sample_pos]),
      threshold = threshold,
      stringsAsFactors = FALSE
    )
  }

  pieces <- list()

  # cond 1 failure: samples with M_i(k) = 0 (undefined log-ratios downstream)
  if (isTRUE(require_no_zero_mass) && !isTRUE(conditions$cond1_zero_mass_ok)) {
    pieces <- c(pieces, list(
      build_rows(which(reference_mass == 0), "zero_reference_mass", 0)
    ))
  }

  # cond 2 failure: samples below the absolute count threshold
  # (note: this reports the failing TAIL, not just the defining sample)
  if (!isTRUE(conditions$cond2_bulk_count_ok)) {
    pieces <- c(pieces, list(
      build_rows(which(reference_mass < min_bulk_count),
                 "low_absolute_reference_mass", min_bulk_count)
    ))
  }

  # cond 3 failure: samples below the relative reference mass threshold
  if (!isTRUE(conditions$cond3_bulk_rel_ok)) {
    pieces <- c(pieces, list(
      build_rows(which(relative_reference_mass < min_bulk_rel),
                 "low_relative_reference_mass", min_bulk_rel)
    ))
  }

  pieces <- Filter(Negate(is.null), pieces)
  # all conditions OK (or all failing sets empty) -> return the empty template
  # so callers can always rbind / nrow() without special casing
  if (length(pieces) == 0L) return(make_empty_reference_mass_problem_table())

  do.call(rbind, pieces)
}

# Format a small summary string for samples that failed reference mass adequacy
# Used by select_reference_bootstrap_prefix() for verbose messages and inside
# the budget exhausted warning
format_reference_mass_problem_samples <- function(problem_table, max_samples = 12L) {
  # empty table -> empty string (callers check nzchar() before printing)
  if (is.null(problem_table) || nrow(problem_table) == 0L) return("")

  # group reason rows by sample so each sample appears once with all its reasons
  split_reasons <- split(problem_table$reason, problem_table$sample_name)

  sample_strings <- vapply(
    names(split_reasons),
    function(sample_name) {
      reasons <- paste(unique(split_reasons[[sample_name]]), collapse = ", ")
      paste0(sample_name, " [", reasons, "]")
    },
    character(1L)
  )

  # cap output length. long lists become "first 12 ... + N more"
  n_total <- length(sample_strings)
  if (n_total > max_samples) {
    sample_strings <- c(
      sample_strings[seq_len(max_samples)],
      paste0("... ", n_total - max_samples, " more")
    )
  }

  paste(sample_strings, collapse = "; ")
}


# Walk the prefix forward when k0 still leaves small mass.
# Called from select_reference_bootstrap_prefix() right after find_stable_prefix()
# resolves the smallest feasible k0;
#
#  1. within_budget : grow past k0 inside size/score/stability caps.
#                     stop when all conditions clear

#  2. override : if budget runs out and cond1 (no zero-mass) still fails
#                AND require_no_zero_mass = TRUE, keep growing
#                past budget specifically for cond1, then warn.

#  3. hard fail : if we reach the end of the candidate list still failing
#                 cond1, stop() - those samples cant be
#                 supported by any subset and need to be dropped.
extend_fragile_reference_prefix <- function(
    selection_table,
    sorted_scores,
    n_samples,
    k0,
    X,
    full_sorted_candidate_taxa,
    require_no_zero_mass = TRUE,
    min_bulk_count = 10,
    min_bulk_rel = 0.01,
    hostage_max_size_inflation = 0.50,
    hostage_max_score_inflation = 0.25,
    hostage_max_stability_drop = 0.05
) {
  # locate the k0 row in the merged selection table (built by find_stable_prefix)
  base_row <- selection_table[selection_table$k == k0, , drop = FALSE]
  if (nrow(base_row) != 1L) stop("k0 not found in selection_table.")

  base_stability <- base_row$stability         # S_{k_star} (denominator for stability drop)
  base_score <- sorted_scores[k0]              # sigma_{k_star} (denominator for score inflation)
  if (!is.finite(base_score) || base_score <= 0) base_score <- 1e-8   # avoid divide by zero on perfect scores

  # adequacy + per-sample diagnostics at k0 (used as the base reference for everything below)
  base_cond <- evaluate_adequacy_conditions(
    adequacy_row = base_row,
    require_no_zero_mass = require_no_zero_mass,
    min_bulk_count = min_bulk_count,
    min_bulk_rel = min_bulk_rel
  )

  # which samples fail at k0, this is what gets returned as expansion_trigger_samples
  base_problem_samples <- diagnose_prefix_mass_problem_samples(
    X = X,
    full_sorted_candidate_taxa = full_sorted_candidate_taxa,
    k = k0,
    conditions = base_cond,
    require_no_zero_mass = require_no_zero_mass,
    min_bulk_count = min_bulk_count,
    min_bulk_rel = min_bulk_rel
  )

  # if k0 already passes all three conditions, no expansion needed
  if (base_cond$all_ok) {
    return(list(
      selected_k = k0,
      hostage_base_k = k0,
      hostage_reason = "all_conditions_ok_at_bulk_cutoff",
      override_invoked = FALSE,
      override_samples = character(0),
      expansion_trigger_samples = make_empty_reference_mass_problem_table(),
      compromise_samples = make_empty_reference_mass_problem_table(),
      final_size_inflation = 0,
      final_score_inflation = 0,
      final_stability_drop = 0,
      final_cond1_ok = base_cond$cond1_zero_mass_ok,
      final_cond2_ok = base_cond$cond2_bulk_count_ok,
      final_cond3_ok = base_cond$cond3_bulk_rel_ok
    ))
  }

  # walk forward through k > k0 candidates in order
  later_prefix_sizes <- selection_table$k[selection_table$k > k0]
  selected_k <- k0                              # current best k, overwritten when we accept an expansion
  best_cond <- base_cond                        # adequacy state at selected_k
  selection_reason <- "no_expansion_attempted"  # overwritten by the loop exit branch
  override_invoked <- FALSE
  mode <- "within_budget"                       # state: within_budget -> override_for_condition_1

  # last-accepted collateral metrics (for the callers warning output)
  final_size_inflation <- 0
  final_score_inflation <- 0
  final_stability_drop <- 0

  for (prefix_size in later_prefix_sizes) {
    candidate_row <- selection_table[selection_table$k == prefix_size, , drop = FALSE]

    cand_cond <- evaluate_adequacy_conditions(
      adequacy_row = candidate_row,
      require_no_zero_mass = require_no_zero_mass,
      min_bulk_count = min_bulk_count,
      min_bulk_rel = min_bulk_rel
    )

    # Expansion budget for k_prime vs base k_star
    size_inflation  <- (prefix_size - k0) / max(k0, 1) # size inflation  = (k_prime - k_star) / k_star
    score_inflation <- (sorted_scores[prefix_size] - base_score) / base_score # score inflation = (sigma_k_prime - sigma_k_star) / sigma_k_star
    stability_drop  <- base_stability - candidate_row$stability # stability drop  = S_{k_star} - S_{k_prime}

    within_budget <- (size_inflation  <= hostage_max_size_inflation) &&
      (score_inflation <= hostage_max_score_inflation) &&
      (stability_drop  <= hostage_max_stability_drop)

    if (identical(mode, "within_budget")) {
      if (within_budget) {
        # accept this prefix, keep going if conditions arent all clear yet
        selected_k <- prefix_size
        best_cond <- cand_cond
        final_size_inflation <- size_inflation
        final_score_inflation <- score_inflation
        final_stability_drop <- stability_drop

        if (cand_cond$all_ok) {
          selection_reason <- "all_conditions_rescued_within_budget"
          break
        }
      } else {
        # budget exhausted - decide whether to enter override
        if (!best_cond$cond1_zero_mass_ok && require_no_zero_mass) {
          mode <- "override_for_condition_1"
          override_invoked <- TRUE
          # reevaluate this prefix under override rules below
        } else {
          selection_reason <- if (best_cond$all_ok) {
            "all_conditions_rescued_within_budget"
          } else {
            "budget_exhausted_quality_compromise"
          }
          break
        }
      }
    }

    if (identical(mode, "override_for_condition_1")) {
      # override: keep expanding past budget until cond1 clears
      selected_k <- prefix_size
      best_cond <- cand_cond
      final_size_inflation <- size_inflation
      final_score_inflation <- score_inflation
      final_stability_drop <- stability_drop

      if (cand_cond$cond1_zero_mass_ok) {
        selection_reason <- "condition_1_rescued_via_override"
        break
      }
    }
  }

  # loop exited without setting a specific reason. Happens if we ran out
  # of later_prefix_sizes (candidates exhausted) before any break fired
  if (selection_reason == "no_expansion_attempted") {
    selection_reason <- if (identical(mode, "override_for_condition_1")) {
      "candidates_exhausted_in_override"
    } else if (best_cond$all_ok) {
      "all_conditions_rescued_within_budget"
    } else {
      "budget_exhausted_quality_compromise"
    }
  }

  # hard fail: cond1 still bad after exhausting candidates.
  # Some sample has zero counts across every candidate reference taxon
  if (require_no_zero_mass && !best_cond$cond1_zero_mass_ok) {
    final_reference_indices <- full_sorted_candidate_taxa[seq_len(selected_k)]
    final_reference_mass <- rowSums(X[, final_reference_indices, drop = FALSE])
    zero_mass_samples <- rownames(X)[final_reference_mass == 0]

    stop(
      "Reference selection failed: condition 1 (no zero-mass samples) could not be ",
      "satisfied even at maximum expansion. ",
      length(zero_mass_samples), " sample(s) have zero reference mass at the largest ",
      "tested prefix size (k = ", selected_k, "): ",
      paste(zero_mass_samples, collapse = ", "), ". ",
      "Remove these samples and re-run."
    )
  }

  # which samples specifically needed override? compare k0 mass vs final mass:
  # those that had M_i(k0) = 0 but now have M_i(selected_k) > 0 are the rescued ones
  override_samples <- character(0)
  if (override_invoked) {
    base_reference_indices <- full_sorted_candidate_taxa[seq_len(k0)]
    base_reference_mass <- rowSums(X[, base_reference_indices, drop = FALSE])
    selected_reference_indices <- full_sorted_candidate_taxa[seq_len(selected_k)]
    selected_reference_mass <- rowSums(X[, selected_reference_indices, drop = FALSE])

    override_samples <- rownames(X)[(base_reference_mass == 0) & (selected_reference_mass > 0)]
  }

  # per sample diagnostics at the FINAL selected_k (for compromise reporting)
  final_problem_samples <- diagnose_prefix_mass_problem_samples(
    X = X,
    full_sorted_candidate_taxa = full_sorted_candidate_taxa,
    k = selected_k,
    conditions = best_cond,
    require_no_zero_mass = require_no_zero_mass,
    min_bulk_count = min_bulk_count,
    min_bulk_rel = min_bulk_rel
  )

  # only report compromise_samples when we actually compromised on quality
  # otherwise leave the empty template so the caller knows nothing went wrong
  compromise_samples <- if (identical(selection_reason, "budget_exhausted_quality_compromise")) {
    final_problem_samples
  } else {
    make_empty_reference_mass_problem_table()
  }

  list(
    selected_k = selected_k,
    hostage_base_k = k0,
    hostage_reason = selection_reason,
    override_invoked = override_invoked,
    override_samples = override_samples,
    expansion_trigger_samples = base_problem_samples,
    compromise_samples = compromise_samples,
    final_size_inflation = final_size_inflation,
    final_score_inflation = final_score_inflation,
    final_stability_drop = final_stability_drop,
    final_cond1_ok = best_cond$cond1_zero_mass_ok,
    final_cond2_ok = best_cond$cond2_bulk_count_ok,
    final_cond3_ok = best_cond$cond3_bulk_rel_ok
  )
}


# -----------------------------------------------------------------------------
# Group-aware reference cleanup
# -----------------------------------------------------------------------------

# Group-aware reference cleanup: drop the top `frac` references by group-shift
# score, with a floor of `min_ref` retained. WIP, binary tested_term only.
cleanup_reference_group_shift <- function(
    X,
    ref_idx,
    metadata,
    tested_term,
    frac = 0.15,
    min_ref = 10L,
    pseudocount = 1
) {
  ref_idx <- as.integer(unique(ref_idx))   # tolerate caller passing repeats

  # not enough references to drop anything safely -> return untouched with NA scores
  if (length(ref_idx) <= min_ref) {
    return(list(
      selected_references = ref_idx,
      removed_references = integer(0),
      group_cleanup_scores = setNames(rep(NA_real_, length(ref_idx)), colnames(X)[ref_idx])
    ))
  }

  # delegate to score_reference_group_shift (returns D_r per reference taxon)
  scores <- score_reference_group_shift(
    X = X, ref_idx = ref_idx, metadata = metadata,
    tested_term = tested_term, pseudocount = pseudocount
  )

  # how many to drop: floor(N * frac), but never go below min_ref retained
  n_drop <- min(floor(length(ref_idx) * frac), length(ref_idx) - min_ref)

  # frac too small to drop anyone -> return untouched (still keep the scores for diagnostics)
  if (n_drop <= 0L) {
    return(list(
      selected_references = ref_idx,
      removed_references = integer(0),
      group_cleanup_scores = scores
    ))
  }

  # decreasing order: largest D_r first, drop the top n_drop
  score_order <- order(scores, decreasing = TRUE, na.last = TRUE)
  drop_idx <- ref_idx[score_order[seq_len(n_drop)]]
  keep_idx <- setdiff(ref_idx, drop_idx)

  list(
    selected_references = keep_idx,
    removed_references = drop_idx,
    group_cleanup_scores = scores             # scores for ALL refs, not just the kept ones
  )
}


# Per reference taxon r: log ratio of r vs the rest of the selected set,
# then absolute group mean difference
score_reference_group_shift <- function(
    X,
    ref_idx,
    metadata,
    tested_term,
    pseudocount = 1
) {

  group_factor <- metadata[[tested_term]]
  group_factor <- droplevels(as.factor(group_factor))

  if (nlevels(group_factor) != 2L) {
    stop("group-aware reference cleanup expects a binary tested_term.")
  }

  group_levels <- levels(group_factor)

  group1_rows <- which(group_factor == group_levels[1L])
  group2_rows <- which(group_factor == group_levels[2L])

  scores <- numeric(length(ref_idx))
  names(scores) <- colnames(X)[ref_idx]

  for (ref_pos in seq_along(ref_idx)) {
    reference_index <- ref_idx[ref_pos]
    other_references <- setdiff(ref_idx, reference_index)

    if (length(other_references) < 1L) {
      scores[ref_pos] <- 0
      next
    }

    other_reference_mass <- rowSums(X[, other_references, drop = FALSE])    # sum_{q in R, q != r} X_iq
    log_reference_ratio <- log(
      (X[, reference_index] + pseudocount) /
        (other_reference_mass + pseudocount)
    )                                                                       # G_ir = log((X_ir + c) / (sum_{q in R, q != r} X_iq + c))

    # Gbar_rg = (1 / n_g) sum_{i : T_i = g} G_ir, then D_r = | Gbar_r1 - Gbar_r0 |
    scores[ref_pos] <- abs(
      mean(log_reference_ratio[group1_rows], na.rm = TRUE) -
        mean(log_reference_ratio[group2_rows], na.rm = TRUE)
    )
  }

  scores
}
