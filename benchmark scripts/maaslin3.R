suppressPackageStartupMessages({
  suppressMessages(suppressWarnings(invisible(capture.output(
    source("benchmark_engine.R"),
    file = nullfile()
  ))))
  library(maaslin3)
})

run_maaslin3_method <- function(otu_full,
                                otu_filt,
                                meta_full,
                                meta_filt,
                                regime_row,
                                alpha,
                                tested_metadata,
                                method_args) {
  formula <- if (regime_row$confounder.type == "none") {
    "~ group + log_depth"
  } else {
    "~ group + confounder + log_depth"
  }

  input_data <- as.data.frame(t(otu_filt))
  input_metadata <- as.data.frame(meta_filt)
  input_metadata <- input_metadata[rownames(input_data), , drop = FALSE]
  input_metadata$group <- factor(as.character(input_metadata$group), levels = c("control", "case"))

  sample_depth <- colSums(otu_full)
  sample_depth <- sample_depth[rownames(input_metadata)]
  input_metadata$log_depth <- log(pmax(as.numeric(sample_depth), 1))

  tmp_out <- file.path(
    tempdir(),
    paste0(
      "maaslin3_", regime_row$regime_id, "_",
      Sys.getpid(), "_", as.integer(stats::runif(1, 1, 1e9))
    )
  )
  dir.create(tmp_out, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_out, recursive = TRUE, force = TRUE), add = TRUE)

  invisible(utils::capture.output(
    suppressMessages(suppressWarnings(
      do.call(
        maaslin3::maaslin3,
        c(
          list(
            input_data = input_data,
            input_metadata = input_metadata,
            output = tmp_out,
            formula = formula,
            normalization = "TSS",
            transform = "LOG",
            augment = TRUE,
            standardize = TRUE,
            median_comparison_abundance = TRUE,
            median_comparison_prevalence = FALSE,
            plot_summary_plot = FALSE,
            plot_associations = FALSE,
            save_models = FALSE,
            max_significance = 1,
            min_abundance = 0,
            min_prevalence = 0,
            max_prevalence = 1.01,
            min_variance = 0,
            verbosity = "ERROR"
          )
        )
      )
    )),
    type = "output"
  ))

  res <- read.delim(file.path(tmp_out, "all_results.tsv"), sep = "\t", check.names = FALSE)

  res <- res[
    res$metadata == "group" &
      res$value == "case" &
      res$name == "groupcase",
    ,
    drop = FALSE
  ]

  p_df <- suppressWarnings(aggregate(
    as.numeric(res[["pval_joint"]]),
    by = list(taxon_name = res$feature),
    FUN = function(x) min(x, na.rm = TRUE)
  ))
  names(p_df) <- c("taxon_name", "pval")

  q_df <- suppressWarnings(aggregate(
    as.numeric(res[["qval_joint"]]),
    by = list(taxon_name = res$feature),
    FUN = function(x) min(x, na.rm = TRUE)
  ))
  names(q_df) <- c("taxon_name", "qval")

  feature_stats <- merge(p_df, q_df, by = "taxon_name", all = TRUE, sort = FALSE)

  list(
    feature_stats = feature_stats,
    analysis_formula = formula
  )
}

maaslin3_adapter <- list(
  method_name = "MaAsLin3",
  required_packages = c("maaslin3"),
  run_method = run_maaslin3_method
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

  invisible(utils::capture.output(
    suppressMessages(suppressWarnings(
      run_method_benchmark(
        adapter = maaslin3_adapter,
        ref_otu_list = ref_otu_list,
        model_paras_list = model_paras_list,
        n_times = 100,
        seeds = 1:100,
        parallelize = TRUE,
        n_cores = 48,
        alpha = 0.05,
        tested_metadata = "group",
        min_prevalence_filter = 0.05,
        min_total_count_filter = 10,
        sim_background = default_sim_background(),
        method_args = list(),
        out_dir = "maaslin3_25reg_benchmark_output"
      )
    )),
    type = "output"
  ))
}

if (sys.nframe() == 0) {
  main()
}
