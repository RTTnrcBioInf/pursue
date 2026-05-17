
suppressPackageStartupMessages({
  suppressMessages(suppressWarnings(invisible(capture.output(
    source("benchmark_engine.R"),
    file = nullfile()
  ))))
  library(corncob)
})

run_corncob_method <- function(otu_full,
                               otu_filt,
                               meta_full,
                               meta_filt,
                               regime_row,
                               alpha,
                               tested_metadata,
                               method_args) {
  meta_cc <- as.data.frame(meta_filt)
  otu_cc <- as.matrix(otu_filt)
  meta_cc <- meta_cc[colnames(otu_cc), , drop = FALSE]
  meta_cc$group <- factor(meta_cc$group, levels = c("control", "case"))
  
  if (!identical(rownames(meta_cc), colnames(otu_cc))) {
    stop("corncob metadata/count sample alignment failed before fitting.")
  }
  
  if (regime_row$confounder.type == "none") {
    formula <- ~ group
    formula_null <- ~ 1
    phi.formula <- ~ group
    phi.formula_null <- ~ group
    analysis_formula <- "abundance: ~ group vs ~ 1 ; dispersion fixed at ~ group"
  } else {
    formula <- ~ group + confounder
    formula_null <- ~ confounder
    phi.formula <- ~ group + confounder
    phi.formula_null <- ~ group + confounder
    analysis_formula <- "abundance: ~ group + confounder vs ~ confounder ; dispersion fixed at ~ group + confounder"
  }
  
  out_cc <- NULL
  invisible(capture.output({
    suppressWarnings(suppressMessages({
      out_cc <- do.call(
        corncob::differentialTest,
        c(
          list(
            formula = formula,
            phi.formula = phi.formula,
            formula_null = formula_null,
            phi.formula_null = phi.formula_null,
            data = otu_cc,
            sample_data = meta_cc,
            taxa_are_rows = TRUE
          ),
          method_args
        )
      )
    }))
  }, file = nullfile()))
  
  p_vec <- out_cc[["p"]]
  q_vec <- out_cc[["p_fdr"]]
  taxa_names <- names(q_vec)
  
  # corncob reports two separate "discriminant taxa" sets.
  # discriminant_taxa_DA = taxa where at least one abundance-model covariate
  # is perfectly discriminant, i.e. the covariate structure separates the taxon
  # so cleanly that the differential-abundance fit is not a normal comparable
  # test result.
  # discriminant_taxa_DV = the same idea, but for the dispersion/variability
  # model.
  #
  # The benchmark only needs one exclusion vector, so the DA and DV lists are
  # merged. These taxa are still flagged in feature_stats, and are also passed
  # to exclude_from_truth_comparison so they do not get scored as ordinary
  # positives/negatives in the simulation truth comparison.
  discriminant_taxa <- unique(c(
    as.character(out_cc[["discriminant_taxa_DA"]]),
    as.character(out_cc[["discriminant_taxa_DV"]])
  ))
  
  feature_stats <- data.frame(
    taxon_name = taxa_names,
    pval = suppressWarnings(as.numeric(p_vec)),
    qval = suppressWarnings(as.numeric(q_vec)),
    called_native = as.integer(taxa_names %in% out_cc[["significant_taxa"]]),
    flagged_discriminant_native = as.integer(taxa_names %in% discriminant_taxa),
    stringsAsFactors = FALSE
  )
  
  list(
    feature_stats = feature_stats,
    analysis_formula = analysis_formula,
    missing_taxa_policy = "negative",
    exclude_from_truth_comparison = discriminant_taxa
  )
}

corncob_adapter <- list(
  method_name = "corncob",
  required_packages = c("corncob"),
  run_method = run_corncob_method
)

main <- function() {
  res <- NULL
  invisible(capture.output({
    suppressWarnings(suppressMessages({
      res <- run_method_benchmark(
        adapter = corncob_adapter,
        ref_otu_list = list(
          stool = read.csv("stool_otu.csv", row.names = 1, check.names = FALSE),
          vaginal = read.csv("vaginal_otu.csv", row.names = 1, check.names = FALSE)
        ),
        model_paras_list = list(
          stool = readRDS("stool_model_paras.rds"),
          vaginal = readRDS("vaginal_model_paras.rds")
        ),
        n_times = 2,
        seeds = 1:2,
        parallelize = TRUE,
        n_cores = 2,
        alpha = 0.05,
        tested_metadata = "group",
        min_prevalence_filter = 0.05,
        min_total_count_filter = 10,
        sim_background = default_sim_background(),
        method_args = list(
          test = "Wald",
          boot = FALSE,
          fdr_cutoff = 0.05,
          filter_discriminant = TRUE,
          verbose = FALSE
        ),
        out_dir = "corncob_25reg_benchmark_output"
      )
    }))
  }, file = nullfile()))
  
  invisible(res)
}

if (sys.nframe() == 0) {
  invisible(main())
}
