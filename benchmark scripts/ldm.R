
suppressPackageStartupMessages({
  suppressMessages(suppressWarnings(invisible(capture.output(
    source("benchmark_engine.R"),
    file = nullfile()
  ))))
  library(LDM)
})

run_ldm_method <- function(otu_full,
                           otu_filt,
                           meta_full,
                           meta_filt,
                           regime_row,
                           alpha,
                           tested_metadata,
                           method_args) {
  y <- t(otu_filt)
  storage.mode(y) <- "numeric"

  meta_ldm <- as.data.frame(meta_filt)
  meta_ldm$group <- factor(meta_ldm$group, levels = c("control", "case"))

  assign("y", y, envir = .GlobalEnv)
  on.exit(rm(y, envir = .GlobalEnv), add = TRUE)

  formula_obj <- if (as.character(regime_row$confounder.type) == "none") {
    y ~ group
  } else {
    y | confounder ~ group
  }

  out_ldm <- suppressWarnings(suppressMessages({
    fit <- NULL
    utils::capture.output(
      fit <- do.call(
        LDM::ldm,
        c(
          list(
            formula = formula_obj,
            data = meta_ldm
          ),
          method_args
        )
      )
    )
    fit
  }))

  p_vec <- as.numeric(as.matrix(out_ldm[["p.otu.omni"]])["cov1", ])
  q_vec <- as.numeric(as.matrix(out_ldm[["q.otu.omni"]])["cov1", ])

  feature_stats <- data.frame(
    taxon_name = colnames(y),
    pval = p_vec,
    qval = q_vec,
    stringsAsFactors = FALSE
  )

  detected_omni <- out_ldm[["detected.otu.omni"]][["cov1"]]
  feature_stats$called_native <- as.integer(feature_stats$taxon_name %in% detected_omni)

  list(
    feature_stats = feature_stats,
    analysis_formula = paste(deparse(formula_obj), collapse = "")
  )
}

ldm_adapter <- list(
  method_name = "LDM",
  required_packages = c("LDM"),
  run_method = run_ldm_method
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

  res <- suppressWarnings(suppressMessages({
    benchmark <- NULL
    utils::capture.output(
      benchmark <- run_method_benchmark(
        adapter = ldm_adapter,
        ref_otu_list = ref_otu_list,
        model_paras_list = model_paras_list,
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
          dist.method = "bray",
          n.perm.max = 20000,
          seed = 1,
          fdr.nominal = 0.05,
          scale.otu.table = TRUE,
          center.otu.table = TRUE,
          freq.scale.only = FALSE,
          binary = FALSE,
          n.rarefy = 0,
          test.omni3 = FALSE,
          comp.anal = FALSE,
          n.cores = 1,
          verbose = FALSE
        ),
        out_dir = "ldm_25reg_benchmark_output"
      )
    )
    benchmark
  }))

  invisible(res)
}

if (sys.nframe() == 0) {
  main()
}
