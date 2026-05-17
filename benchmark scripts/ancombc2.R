
suppressPackageStartupMessages({
  suppressMessages(suppressWarnings(invisible(capture.output(
    source("benchmark_engine.R"),
    file = nullfile()
  ))))
  library(ANCOMBC)
})

run_ancombc2_method <- function(otu_full,
                                otu_filt,
                                meta_full,
                                meta_filt,
                                regime_row,
                                alpha,
                                tested_metadata,
                                method_args) {
  formula <- if (regime_row$confounder.type == "none") {
    "group"
  } else {
    "group + confounder"
  }
  
  meta_filt <- as.data.frame(meta_filt)
  meta_filt <- meta_filt[colnames(otu_filt), , drop = FALSE]
  meta_filt$group <- factor(meta_filt$group, levels = c("control", "case"))
  
  
  out_ancom <- NULL
  invisible(capture.output({
    suppressWarnings(suppressMessages({
      out_ancom <- do.call(
        ANCOMBC::ancombc2,
        c(
          list(
            data = otu_filt,
            taxa_are_rows = TRUE,
            meta_data = meta_filt,
            fix_formula = formula
          ),
          method_args
        )
      )
    }))
  }, file = nullfile()))
  
  res_df <- out_ancom$res
  
  called_robust <- as.logical(res_df[["diff_robust_groupcase"]])
  p_raw <- suppressWarnings(as.numeric(res_df[["p_groupcase"]]))
  q_raw <- suppressWarnings(as.numeric(res_df[["q_groupcase"]]))
  
  feature_stats <- data.frame(
    taxon_name = as.character(res_df[["taxon"]]),
    pval = ifelse(called_robust & !is.na(p_raw), p_raw, 1),
    qval = ifelse(called_robust & !is.na(q_raw), q_raw, 1),
    pval_raw_ancombc2 = p_raw,
    qval_raw_ancombc2 = q_raw,
    called_ancombc2_raw = as.logical(res_df[["diff_groupcase"]]),
    passed_ss_ancombc2 = as.logical(res_df[["passed_ss_groupcase"]]),
    called_ancombc2_robust = called_robust,
    called_ancom = called_robust,
    lfc_ref = suppressWarnings(as.numeric(res_df[["lfc_groupcase"]])),
    stringsAsFactors = FALSE
  )
  
  list(
    feature_stats = feature_stats,
    analysis_formula = paste0(formula, "; decision = diff_robust_groupcase")
  )
}

ancombc2_adapter <- list(
  method_name = "ANCOM-BC2-SS",
  required_packages = c("ANCOMBC"),
  run_method = run_ancombc2_method
)

main <- function() {
  res <- NULL
  invisible(capture.output({
    suppressWarnings(suppressMessages({
      res <- run_method_benchmark(
        adapter = ancombc2_adapter,
        ref_otu_list = list(
          stool = read.csv("stool_otu.csv", row.names = 1, check.names = FALSE),
          vaginal = read.csv("vaginal_otu.csv", row.names = 1, check.names = FALSE)
        ),
        model_paras_list = list(
          stool = readRDS("stool_model_paras.rds"),
          vaginal = readRDS("vaginal_model_paras.rds")
        ),
        n_times = 100,
        seeds = 1:100,
        parallelize = TRUE,
        n_cores = 50,
        alpha = 0.05,
        tested_metadata = "group",
        min_prevalence_filter = 0.05,
        min_total_count_filter = 10,
        sim_background = default_sim_background(),
        method_args = list(
          rand_formula = NULL,
          p_adj_method = "BH",
          pseudo = 0,
          pseudo_sens = TRUE,
          prv_cut = 0,
          lib_cut = 0,
          s0_perc = 0.05,
          group = "group",
          struc_zero = FALSE,
          neg_lb = FALSE,
          alpha = 0.05,
          n_cl = 1,
          verbose = FALSE,
          global = FALSE,
          pairwise = FALSE,
          dunnet = FALSE,
          trend = FALSE
        ),
        out_dir = "ancombc2_25reg_benchmark_output_2"
      )
    }))
  }, file = nullfile()))
  
  invisible(res)
}

if (sys.nframe() == 0) {
  main()
}
