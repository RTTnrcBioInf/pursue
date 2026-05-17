library(GUniFrac)
library(parallel)

# Obtaining stool and vaginal reference OTU tables and model parameters
#
# stool_ref <- data(stool.otu.tab)
# vaginal_ref <- data(vaginal.otu.tab)
# 
# stool_paras   <- GUniFrac:::EstPara(stool_ref)
# vaginal_paras <- GUniFrac:::EstPara(vaginal_ref)

sf_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) 0 else stats::sd(x)
}

sf_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else mean(x)
}

sf_median <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else stats::median(x)
}

prefilter_features_raw <- function(count_mat,
                                   min_prevalence = 0.05,
                                   min_total_count = 10) {
  keep_prev  <- rowMeans(count_mat > 0) >= min_prevalence
  keep_total <- rowSums(count_mat) >= min_total_count
  keep <- keep_prev & keep_total
  count_mat[keep, , drop = FALSE]
}

calc_run_metrics <- function(qvals, truth_logical, alpha = 0.05) {
  qvals <- as.numeric(qvals)
  qvals[!is.finite(qvals)] <- 1
  truth_logical <- as.logical(truth_logical)
  
  pred <- qvals < alpha
  
  tp <- sum(pred & truth_logical)
  fp <- sum(pred & !truth_logical)
  fn <- sum(!pred & truth_logical)
  tn <- sum(!pred & !truth_logical)
  
  precision <- if ((tp + fp) == 0) 0 else tp / (tp + fp)
  recall    <- if ((tp + fn) == 0) 0 else tp / (tp + fn)
  f1        <- if ((precision + recall) == 0) 0 else {
    2 * precision * recall / (precision + recall)
  }
  fpr <- if ((fp + tn) == 0) NA_real_ else fp / (fp + tn)
  fdr <- if ((tp + fp) == 0) 0 else fp / (tp + fp)
  
  data.frame(
    tp = tp,
    fp = fp,
    fn = fn,
    tn = tn,
    fpr = fpr,
    fdr = fdr,
    precision = precision,
    recall = recall,
    f1 = f1,
    stringsAsFactors = FALSE
  )
}

make_all_regimes <- function(grp.ratio = 1) {
  null_grid <- expand.grid(
    diversity = c("stool", "vaginal"),
    nSam = c(50, 100, 200),
    nOTU = c(50, 500),
    stringsAsFactors = FALSE
  )
  
  null_grid$regime_family   <- "global_null"
  null_grid$regime_name     <- paste0(
    "global_null__", null_grid$diversity,
    "__n", null_grid$nSam,
    "__m", null_grid$nOTU
  )
  null_grid$is_null         <- TRUE
  null_grid$diff.otu.pct    <- 0
  null_grid$diff.otu.direct <- "balanced"
  null_grid$diff.otu.mode   <- "abundant"
  null_grid$depth.conf.factor <- 0
  null_grid$confounder.type <- "none"
  
  null_grid <- null_grid[, c(
    "regime_family", "regime_name", "diversity", "nSam", "nOTU",
    "is_null", "diff.otu.pct", "diff.otu.direct", "diff.otu.mode",
    "depth.conf.factor", "confounder.type"
  )]
  
  nonnull <- data.frame(
    regime_family = c(
      "baseline_nonnull",
      "unbalanced_shift",
      "low_diversity",
      "sample_size",
      "sample_size",
      "small_taxa_number",
      "depth_confounding",
      "depth_confounding",
      "signal_density",
      "signal_density",
      "diff_otu_rarity",
      "diff_otu_rarity",
      "confounder_adjustment"
    ),
    regime_name = c(
      "baseline_nonnull",
      "unbalanced_shift",
      "low_diversity",
      "sample_size__n50",
      "sample_size__n200",
      "small_taxa_number",
      "depth_confounding__4x",
      "depth_confounding__9x",
      "signal_density__5pct",
      "signal_density__20pct",
      "diff_otu_rarity__rare",
      "diff_otu_rarity__mix",
      "confounder_adjustment"
    ),
    diversity = c(
      "stool","stool","vaginal","stool","stool","stool",
      "stool","stool","stool","stool","stool","stool","stool"
    ),
    nSam = c(
      100,100,100,50,200,100,100,100,100,100,100,100,100
    ),
    nOTU = c(
      500,500,500,500,500,50,500,500,500,500,500,500,500
    ),
    is_null = FALSE,
    diff.otu.pct = c(
      0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.05,0.20,0.10,0.10,0.10
    ),
    diff.otu.direct = c(
      "balanced","unbalanced","balanced","balanced","balanced","balanced",
      "balanced","balanced","balanced","balanced","balanced","balanced","balanced"
    ),
    diff.otu.mode = c(
      "abundant","abundant","abundant","abundant","abundant","abundant",
      "abundant","abundant","abundant","abundant","rare","mix","abundant"
    ),
    depth.conf.factor = c(
      0,0,0,0,0,0,0.69,1.1,0,0,0,0,0
    ),
    confounder.type = c(
      "none","none","none","none","none","none",
      "none","none","none","none","none","none","continuous"
    ),
    stringsAsFactors = FALSE
  )
  
  regimes <- rbind(null_grid, nonnull)
  
  regimes$regime_id <- sprintf("R%02d", seq_len(nrow(regimes)))
  regimes$grp.ratio <- grp.ratio
  regimes <- regimes[, c(
    "regime_id", "regime_family", "regime_name",
    "is_null", "diversity", "nSam", "nOTU",
    "diff.otu.pct", "diff.otu.direct", "diff.otu.mode",
    "depth.conf.factor",
    "confounder.type", "grp.ratio"
  )]
  
  regimes
}

default_sim_background <- function() {
  list(
    covariate.type = "binary",
    grp.ratio = 1,
    
    covariate.eff.mean = 1,
    covariate.eff.sd   = 0,
    
    error.sd = 0,
    
    depth.mu    = 10000,
    depth.theta = 5,
    
    conf.cov.cor         = 0.6,
    conf.diff.otu.pct    = 0,
    conf.nondiff.otu.pct = 0.1,
    confounder.eff.mean  = 1,
    confounder.eff.sd    = 0
  )
}

build_sim_args_from_regime <- function(regime_row, sim_background) {
  conf_type <- regime_row$confounder.type
  
  list(
    nSam = regime_row$nSam,
    nOTU = regime_row$nOTU,
    
    diff.otu.pct    = regime_row$diff.otu.pct,
    diff.otu.direct = regime_row$diff.otu.direct,
    diff.otu.mode   = regime_row$diff.otu.mode,
    
    covariate.type     = sim_background$covariate.type,
    grp.ratio          = regime_row$grp.ratio,
    covariate.eff.mean = sim_background$covariate.eff.mean,
    covariate.eff.sd   = sim_background$covariate.eff.sd,
    
    error.sd          = sim_background$error.sd,
    depth.mu          = sim_background$depth.mu,
    depth.theta       = sim_background$depth.theta,
    depth.conf.factor = regime_row$depth.conf.factor,
    
    confounder.type      = conf_type,
    conf.cov.cor         = if (conf_type == "none") 0 else sim_background$conf.cov.cor,
    conf.diff.otu.pct    = if (conf_type == "none") 0 else sim_background$conf.diff.otu.pct,
    conf.nondiff.otu.pct = if (conf_type == "none") 0 else sim_background$conf.nondiff.otu.pct,
    confounder.eff.mean  = if (conf_type == "none") 0 else sim_background$confounder.eff.mean,
    confounder.eff.sd    = if (conf_type == "none") 0 else sim_background$confounder.eff.sd
  )
}

default_meta_builder <- function(sim_obj) {
  otu.sim <- sim_obj$otu.tab.sim
  
  x <- sim_obj$covariate
  if (is.matrix(x) || is.data.frame(x)) x <- x[, 1]
  x <- as.numeric(x)
  
  meta.sim <- data.frame(
    group = factor(
      ifelse(x == 1, "case", "control"),
      levels = c("control", "case")
    ),
    row.names = colnames(otu.sim),
    stringsAsFactors = FALSE
  )
  
  conf <- sim_obj$confounder
  if (!is.null(conf)) {
    if (is.matrix(conf) || is.data.frame(conf)) {
      meta.sim$confounder <- as.numeric(conf[, 1])
    } else {
      meta.sim$confounder <- as.numeric(conf)
    }
  }
  
  as.data.frame(meta.sim)
}

simulate_regime_dataset <- function(regime_row,
                                    ref_otu_list,
                                    model_paras_list,
                                    sim_background,
                                    seed,
                                    meta_builder = default_meta_builder) {
  set.seed(seed)
  
  ref_name <- regime_row$diversity
  ref_otu <- ref_otu_list[[ref_name]]
  model.paras <- model_paras_list[[ref_name]]
  
  sim_args <- build_sim_args_from_regime(regime_row, sim_background)
  
  sim_obj <- do.call(
    GUniFrac::SimulateMSeq,
    c(
      list(
        ref.otu.tab = ref_otu,
        model.paras = model.paras
      ),
      sim_args
    )
  )
  
  otu.sim_full <- sim_obj$otu.tab.sim
  all_taxa <- rownames(otu.sim_full)
  
  truth_names <- character(0)
  if (!is.null(sim_obj$diff.otu.ind) && length(sim_obj$diff.otu.ind) > 0) {
    truth_names <- sim_obj$otu.names[sim_obj$diff.otu.ind]
  }
  
  truth_logical_full <- setNames(rep(FALSE, length(all_taxa)), all_taxa)
  truth_logical_full[truth_names] <- TRUE
  
  meta_full <- as.data.frame(meta_builder(sim_obj))
  meta_full <- meta_full[colnames(otu.sim_full), , drop = FALSE]
  
  list(
    sim_obj = sim_obj,
    sim_args = sim_args,
    otu_full = otu.sim_full,
    meta_full = meta_full,
    all_taxa = all_taxa,
    truth_logical_full = truth_logical_full
  )
}

align_method_results_to_full_taxa <- function(all_taxa,
                                              truth_logical_full,
                                              retained_taxa,
                                              feature_stats,
                                              alpha = 0.05,
                                              iter = NULL,
                                              seed = NULL,
                                              missing_taxa_policy = "negative",
                                              exclude_from_truth_comparison = character(0)) {
  retained_taxa <- unique(as.character(retained_taxa))
  exclude_from_truth_comparison <- unique(as.character(exclude_from_truth_comparison))
  exclude_from_truth_comparison <- exclude_from_truth_comparison[
    !is.na(match(exclude_from_truth_comparison, all_taxa))
  ]
  
  full_res <- data.frame(
    taxon_name = all_taxa,
    retained_after_prefilter = all_taxa %in% retained_taxa,
    returned_by_method = all_taxa %in% feature_stats$taxon_name,
    excluded_by_method_ineligibility = all_taxa %in% exclude_from_truth_comparison,
    eligible_for_truth_comparison =
      (all_taxa %in% retained_taxa) &
      !(all_taxa %in% exclude_from_truth_comparison),
    truth = as.integer(truth_logical_full[all_taxa]),
    stringsAsFactors = FALSE
  )
  
  idx <- match(full_res$taxon_name, feature_stats$taxon_name)
  
  extra_cols <- setdiff(colnames(feature_stats), "taxon_name")
  for (nm in extra_cols) {
    full_res[[nm]] <- feature_stats[[nm]][idx]
  }
  
  full_res$pval <- as.numeric(full_res$pval)
  full_res$qval <- as.numeric(full_res$qval)
  
  eligible <- full_res$eligible_for_truth_comparison
  missing_and_eligible <- eligible & !full_res$returned_by_method
  
  if (missing_taxa_policy == "negative") {
    full_res$pval[missing_and_eligible] <- 1
    full_res$qval[missing_and_eligible] <- 1
  } else if (missing_taxa_policy == "exclude") {
    full_res$eligible_for_truth_comparison[missing_and_eligible] <- FALSE
    full_res$excluded_by_method_ineligibility[missing_and_eligible] <- TRUE
  }
  
  eligible <- full_res$eligible_for_truth_comparison
  
  full_res$pval[eligible & (is.na(full_res$pval) | !is.finite(full_res$pval))] <- 1
  full_res$qval[eligible & (is.na(full_res$qval) | !is.finite(full_res$qval))] <- 1
  
  full_res$pval[!eligible] <- NA_real_
  full_res$qval[!eligible] <- NA_real_
  
  full_res$called <- ifelse(
    eligible,
    as.integer(full_res$qval < alpha),
    NA_integer_
  )
  
  if (!is.null(iter)) full_res$iter <- iter
  if (!is.null(seed)) full_res$seed <- seed
  
  lead_cols <- c(
    "iter", "seed",
    "taxon_name",
    "retained_after_prefilter",
    "returned_by_method",
    "excluded_by_method_ineligibility",
    "eligible_for_truth_comparison",
    "truth",
    "pval", "qval", "called"
  )
  lead_cols <- lead_cols[lead_cols %in% colnames(full_res)]
  other_cols <- setdiff(colnames(full_res), lead_cols)
  full_res <- full_res[, c(lead_cols, other_cols), drop = FALSE]
  
  full_res
}

make_regime_summary <- function(run_metrics) {
  metric_cols <- c("tp", "fp", "fn", "tn", "fpr", "fdr", "precision", "recall", "f1")
  split_runs <- split(run_metrics, run_metrics$regime_id)
  
  out <- lapply(split_runs, function(df) {
    row <- df[1, c(
      "regime_id", "regime_family", "regime_name", "is_null",
      "diversity", "nSam", "nOTU",
      "diff.otu.pct", "diff.otu.direct", "diff.otu.mode",
      "depth.conf.factor", "confounder.type"
    ), drop = FALSE]
    
    row$n_runs <- nrow(df)
    row$mean_n_taxa_total           <- sf_mean(df$n_taxa_total)
    row$mean_n_taxa_after_prefilter <- sf_mean(df$n_taxa_after_prefilter)
    row$mean_n_true_signal          <- sf_mean(df$n_true_signal)
    row$mean_runtime_sec            <- sf_mean(df$runtime_sec)
    row$median_runtime_sec          <- sf_median(df$runtime_sec)
    row$sd_runtime_sec              <- sf_sd(df$runtime_sec)
    
    for (m in metric_cols) {
      x <- df[[m]]
      row[[paste0("mean_", m)]]   <- sf_mean(x)
      row[[paste0("median_", m)]] <- sf_median(x)
      row[[paste0("sd_", m)]]     <- sf_sd(x)
    }
    
    row
  })
  
  do.call(rbind, out)
}

save_regime_rds_files <- function(task_results,
                                  regimes,
                                  method_name,
                                  min_prevalence_filter,
                                  min_total_count_filter,
                                  raw_dir) {
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  
  split_res <- split(
    task_results,
    vapply(task_results, `[[`, character(1), "regime_id")
  )
  
  for (rid in names(split_res)) {
    regime_row <- regimes[regimes$regime_id == rid, , drop = FALSE]
    sub <- split_res[[rid]]
    
    run_metrics_regime <- do.call(rbind, lapply(sub, `[[`, "run_metrics"))
    
    feature_list <- Filter(
      Negate(is.null),
      lapply(sub, `[[`, "feature_results")
    )
    
    feature_results_regime <- if (length(feature_list) == 0) {
      data.frame()
    } else {
      do.call(rbind, feature_list)
    }
    
    sim_args_regime <- sub[[1]]$sim_args
    
    obj <- list(
      regime_info = regime_row,
      method_name = method_name,
      sim_args = sim_args_regime,
      prefilter = list(
        min_prevalence_filter = min_prevalence_filter,
        min_total_count_filter = min_total_count_filter
      ),
      run_metrics = run_metrics_regime,
      feature_results = feature_results_regime
    )
    
    saveRDS(
      obj,
      file = file.path(raw_dir, paste0("regime_", rid, ".rds")),
      compress = "xz"
    )
  }
}

run_one_benchmark_task <- function(task_row,
                                   regimes,
                                   ref_otu_list,
                                   model_paras_list,
                                   sim_background,
                                   adapter,
                                   method_args,
                                   alpha = 0.05,
                                   tested_metadata = "group",
                                   min_prevalence_filter = 0.05,
                                   min_total_count_filter = 10,
                                   meta_builder = default_meta_builder) {
  regime_row <- regimes[match(task_row$regime_id, regimes$regime_id), , drop = FALSE]
  regime_row <- regime_row[1, , drop = FALSE]
  
  run_prefix <- regime_row[, c(
    "regime_id", "regime_family", "regime_name", "is_null",
    "diversity", "nSam", "nOTU",
    "diff.otu.pct", "diff.otu.direct", "diff.otu.mode",
    "depth.conf.factor",
    "confounder.type"
  ), drop = FALSE]
  run_prefix$iter <- task_row$iter
  run_prefix$seed <- task_row$seed
  
  sim_dat <- simulate_regime_dataset(
    regime_row = regime_row,
    ref_otu_list = ref_otu_list,
    model_paras_list = model_paras_list,
    sim_background = sim_background,
    seed = task_row$seed,
    meta_builder = meta_builder
  )
  
  otu_filt <- prefilter_features_raw(
    count_mat = sim_dat$otu_full,
    min_prevalence = min_prevalence_filter,
    min_total_count = min_total_count_filter
  )
  
  meta_filt <- sim_dat$meta_full[colnames(otu_filt), , drop = FALSE]
  
  runtime_start <- proc.time()[["elapsed"]]
  adapter_out <- adapter$run_method(
    otu_full = sim_dat$otu_full,
    otu_filt = otu_filt,
    meta_full = sim_dat$meta_full,
    meta_filt = meta_filt,
    regime_row = regime_row,
    alpha = alpha,
    tested_metadata = tested_metadata,
    method_args = method_args
  )
  runtime_sec <- proc.time()[["elapsed"]] - runtime_start
  
  feature_results <- align_method_results_to_full_taxa(
    all_taxa = sim_dat$all_taxa,
    truth_logical_full = sim_dat$truth_logical_full,
    retained_taxa = rownames(otu_filt),
    feature_stats = adapter_out$feature_stats,
    alpha = alpha,
    iter = task_row$iter,
    seed = task_row$seed,
    missing_taxa_policy = if (!is.null(adapter_out$missing_taxa_policy)) {
      adapter_out$missing_taxa_policy
    } else {
      "negative"
    },
    exclude_from_truth_comparison = if (!is.null(adapter_out$exclude_from_truth_comparison)) {
      adapter_out$exclude_from_truth_comparison
    } else {
      character(0)
    }
  )
  
  eligible <- as.logical(feature_results$eligible_for_truth_comparison)
  
  run_metrics <- calc_run_metrics(
    qvals = feature_results$qval[eligible],
    truth_logical = feature_results$truth[eligible],
    alpha = alpha
  )
  
  run_metrics <- cbind(
    run_prefix,
    data.frame(
      n_taxa_total = length(sim_dat$all_taxa),
      n_taxa_after_prefilter = nrow(otu_filt),
      n_taxa_returned_by_method = sum(feature_results$returned_by_method, na.rm = TRUE),
      n_taxa_excluded_method = sum(feature_results$excluded_by_method_ineligibility, na.rm = TRUE),
      n_taxa_evaluated = sum(eligible),
      n_true_signal = sum(feature_results$truth[eligible]),
      n_true_signal_total = sum(sim_dat$truth_logical_full),
      runtime_sec = runtime_sec,
      stringsAsFactors = FALSE
    ),
    run_metrics
  )
  
  list(
    regime_id = regime_row$regime_id,
    sim_args = sim_dat$sim_args,
    run_metrics = run_metrics,
    feature_results = feature_results
  )
}

run_method_benchmark <- function(
    adapter,
    ref_otu_list,
    model_paras_list,
    n_times = 100,
    seeds = NULL,
    parallelize = TRUE,
    n_cores = max(1L, parallel::detectCores() - 1L),
    alpha = 0.05,
    tested_metadata = "group",
    
    min_prevalence_filter = 0.05,
    min_total_count_filter = 10,
    
    sim_background = default_sim_background(),
    method_args = list(),
    
    regimes = NULL,
    meta_builder = default_meta_builder,
    
    out_dir = NULL
) {
  if (is.null(seeds)) {
    seeds <- seq_len(n_times)
  }
  if (is.null(regimes)) {
    regimes <- make_all_regimes(grp.ratio = sim_background$grp.ratio)
  }
  
  if (is.null(out_dir)) {
    out_dir <- paste0(tolower(gsub("[^A-Za-z0-9]+", "_", adapter$method_name)), "_25reg_benchmark_output")
  }
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  raw_dir <- file.path(out_dir, "raw")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  
  tasks <- do.call(
    rbind,
    lapply(seq_len(nrow(regimes)), function(i) {
      data.frame(
        regime_id = regimes$regime_id[i],
        iter = seq_len(n_times),
        seed = seeds,
        stringsAsFactors = FALSE
      )
    })
  )
  one_task <- function(idx) {
    run_one_benchmark_task(
      task_row = tasks[idx, , drop = FALSE],
      regimes = regimes,
      ref_otu_list = ref_otu_list,
      model_paras_list = model_paras_list,
      sim_background = sim_background,
      adapter = adapter,
      method_args = method_args,
      alpha = alpha,
      tested_metadata = tested_metadata,
      min_prevalence_filter = min_prevalence_filter,
      min_total_count_filter = min_total_count_filter,
      meta_builder = meta_builder
    )
  }
  
  task_results <- if (!parallelize || nrow(tasks) == 1L) {
    lapply(seq_len(nrow(tasks)), one_task)
  } else {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    
    pkgs_to_load <- unique(c("GUniFrac", adapter$required_packages))
    invisible(parallel::clusterCall(cl, function(pkgs) {
      for (p in pkgs) {
        suppressPackageStartupMessages(library(p, character.only = TRUE))
      }
      NULL
    }, pkgs_to_load))
    
    parallel::clusterExport(
      cl,
      varlist = c(
        "tasks",
        "regimes",
        "ref_otu_list",
        "model_paras_list",
        "sim_background",
        "adapter",
        "method_args",
        "alpha",
        "tested_metadata",
        "min_prevalence_filter",
        "min_total_count_filter",
        "meta_builder",
        "prefilter_features_raw",
        "calc_run_metrics",
        "build_sim_args_from_regime",
        "simulate_regime_dataset",
        "align_method_results_to_full_taxa",
        "run_one_benchmark_task",
        "one_task"
      ),
      envir = environment()
    )
    
    parallel::parLapply(cl, seq_len(nrow(tasks)), one_task)
  }
  
  run_metrics <- do.call(rbind, lapply(task_results, `[[`, "run_metrics"))
  regime_summary <- make_regime_summary(run_metrics)
  
  write.csv(
    regime_summary,
    file.path(out_dir, "regime_summary.csv"),
    row.names = FALSE
  )
  
  save_regime_rds_files(
    task_results = task_results,
    regimes = regimes,
    method_name = adapter$method_name,
    min_prevalence_filter = min_prevalence_filter,
    min_total_count_filter = min_total_count_filter,
    raw_dir = raw_dir
  )
  
  invisible(list(
    regimes = regimes,
    run_metrics = run_metrics,
    regime_summary = regime_summary,
    out_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  ))
}