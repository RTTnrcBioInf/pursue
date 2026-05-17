
suppressPackageStartupMessages({
  suppressMessages(suppressWarnings(invisible(capture.output(
    source("benchmark_engine.R"),
    file = nullfile()
  ))))
  library(GUniFrac)
})

run_zicoseq_method <- function(otu_full,
                               otu_filt,
                               meta_full,
                               meta_filt,
                               regime_row,
                               alpha,
                               tested_metadata,
                               method_args) {
  meta_zico <- as.data.frame(meta_filt)
  feature_zico <- as.matrix(otu_filt)
  meta_zico <- meta_zico[colnames(feature_zico), , drop = FALSE]
  meta_zico$group <- factor(as.character(meta_zico$group), levels = c("control", "case"))

  call_args <- c(
    list(
      meta.dat = meta_zico,
      feature.dat = feature_zico,
      grp.name = tested_metadata
    ),
    method_args
  )

  if (regime_row$confounder.type != "none") {
    call_args$adj.name <- "confounder"
    analysis_formula <- "ZicoSeq: group + confounder"
  } else {
    call_args$adj.name <- NULL
    analysis_formula <- "ZicoSeq: group"
  }

  out_zico <- suppressWarnings(suppressMessages(capture.output(
    zico_fit <- do.call(GUniFrac::ZicoSeq, call_args)
  )))
  out_zico <- zico_fit

  result_taxa <- names(out_zico[["p.adj.fdr"]])

  feature_stats <- data.frame(
    taxon_name = result_taxa,
    pval = as.numeric(out_zico[["p.raw"]]),
    qval = as.numeric(out_zico[["p.adj.fdr"]]),
    R2_Func1 = as.numeric(out_zico[["R2"]][, 1]),
    stringsAsFactors = FALSE
  )

  list(
    feature_stats = feature_stats,
    analysis_formula = analysis_formula
  )
}

zicoseq_adapter <- list(
  method_name = "ZicoSeq",
  required_packages = c("GUniFrac"),
  run_method = run_zicoseq_method
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

  res <- suppressWarnings(suppressMessages(capture.output(
    benchmark_res <- run_method_benchmark(
      adapter = zicoseq_adapter,
      ref_otu_list = ref_otu_list,
      model_paras_list = model_paras_list,
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
        feature.dat.type = "count",
        prev.filter = 0,
        mean.abund.filter = 0,
        max.abund.filter = 0,
        min.prop = 0,
        is.winsor = TRUE,
        outlier.pct = 0.03,
        winsor.end = "top",
        is.post.sample = TRUE,
        post.sample.no = 25,
        link.func = list(function(x) x^0.5),
        stats.combine.func = max,
        perm.no = 1999,
        strata = NULL,
        ref.pct = 0.5,
        stage.no = 6,
        excl.pct = 0.2,
        is.fwer = FALSE,
        return.feature.dat = TRUE,
        verbose = FALSE
      ),
      out_dir = "zicoseq_25reg_benchmark_output"
    )
  )))

  invisible(benchmark_res)
}

if (sys.nframe() == 0) {
  main()
}
