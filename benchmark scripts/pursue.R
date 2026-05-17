
suppressPackageStartupMessages({
  suppressMessages(suppressWarnings(invisible(capture.output(
    source("benchmark_engine.R"),
    file = nullfile()
  ))))
  library(PURSUE)
})

run_pursue_method <- function(otu_full,
                              otu_filt,
                              meta_full,
                              meta_filt,
                              regime_row,
                              alpha,
                              tested_metadata,
                              method_args) {
  otu_pursue <- as.matrix(t(otu_filt))
  storage.mode(otu_pursue) <- "numeric"
  
  meta_pursue <- as.data.frame(meta_filt)
  meta_pursue <- meta_pursue[rownames(otu_pursue), , drop = FALSE]
  meta_pursue$group <- factor(as.character(meta_pursue$group), levels = c("control", "case"))
  
  rhs_terms <- tested_metadata
  if (regime_row$confounder.type != "none") {
    rhs_terms <- c(rhs_terms, "confounder")
  }
  
  pursue_args <- method_args
  pursue_args$formula <- stats::as.formula(paste("~", paste(rhs_terms, collapse = " + ")))
  pursue_args$tested_term <- tested_metadata
  pursue_args$tested_vars <- tested_metadata
  pursue_args$min_prevalence_filter <- NULL
  pursue_args$min_total_count_filter <- NULL
  
  out <- NULL
  suppressMessages(suppressWarnings(invisible(capture.output(
    out <- do.call(
      PURSUE::run_pursue,
      c(
        list(
          otu = otu_pursue,
          meta = meta_pursue
        ),
        pursue_args
      )
    ),
    file = nullfile()
  ))))
  
  tab <- out$results
  
  feature_stats <- data.frame(
    taxon_name = tab$taxon_name,
    pval = as.numeric(tab$union_p_value),
    qval = as.numeric(tab$union_q_value),
    stringsAsFactors = FALSE
  )
  
  feature_stats$pval[is.na(feature_stats$pval) | !is.finite(feature_stats$pval)] <- 1
  feature_stats$qval[is.na(feature_stats$qval) | !is.finite(feature_stats$qval)] <- 1
  
  list(
    feature_stats = feature_stats,
    analysis_formula = paste(deparse(pursue_args$formula), collapse = ""),
    missing_taxa_policy = "negative",
    exclude_from_truth_comparison = character(0)
  )
}

pursue_adapter <- list(
  method_name = "PURSUE",
  required_packages = c("PURSUE"),
  run_method = run_pursue_method
)

main <- function() {
  ref_otu_list <- list(
    stool = read.csv("stool_otu.csv", row.names = 1, check.names = FALSE),
    vaginal = read.csv("vaginal_otu.csv", row.names = 1, check.names = FALSE)
  )
  
  model_paras_list <- list(
    stool = readRDS("stool_model_paras.rds"),
    vaginal = readRDS("vaginal_model_paras.rds")
  )
  
  suppressMessages(suppressWarnings(invisible(capture.output(
    run_method_benchmark(
      adapter = pursue_adapter,
      ref_otu_list = ref_otu_list,
      model_paras_list = model_paras_list,
      n_times = 100,
      seeds = 1:100,
      parallelize = FALSE,
      n_cores = 48,
      alpha = 0.05,
      tested_metadata = "group",
      min_prevalence_filter = 0.05,
      min_total_count_filter = 10,
      sim_background = default_sim_background(),
      method_args = list(
        formula = ~ group,
        tested_term = "group",
        tested_vars = "group",
        n_perm = 19999,
        min_prevalence_filter = NULL,
        min_total_count_filter = NULL,
        reference_pseudocount = 1,
        bootstrap_B = 200L,
        bootstrap_stability_threshold = 0.80,
        hostage_prevention = TRUE,
        hostage_max_score_inflation = 0.25,
        hostage_max_stability_drop = 0.05,
        reference_group_cleanup = TRUE,
        reference_group_cleanup_frac = 0.15,
        reference_group_cleanup_min_ref = 10L,
        reference_group_cleanup_pseudocount = 1,
        depth_adjust = TRUE,
        prevalence_expected_rarefied_eps = 1e-12,
        parallelize_permutations = FALSE,
        abundance_min_nonzero_n = 3L,
        abund_resid_winsorize = TRUE,
        abund_resid_winsor_quantile = 0.05,
        abund_resid_winsor_min_n = 20L,
        union_acat_weights = c(0.5, 0.5),
        independent_filtering = TRUE,
        independent_filter_quantile = 0.2,
        keep_diagnostics = FALSE,
        verbose = FALSE
      ),
      out_dir = "pursue_25regtest"
    ),
    file = nullfile()
  ))))
}

if (sys.nframe() == 0) {
  main()
}
