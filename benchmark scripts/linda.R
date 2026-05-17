
suppressPackageStartupMessages({
  suppressMessages(suppressWarnings(invisible(capture.output(
    source("benchmark_engine.R"),
    file = nullfile()
  ))))
  library(MicrobiomeStat)
})

run_linda_method <- function(otu_full,
                             otu_filt,
                             meta_full,
                             meta_filt,
                             regime_row,
                             alpha,
                             tested_metadata,
                             method_args) {
  formula_chr <- if (regime_row$confounder.type == "none") {
    "~ group"
  } else {
    "~ group + confounder"
  }
  
  otu_mat <- as.matrix(otu_filt)
  storage.mode(otu_mat) <- "numeric"
  
  meta_linda <- as.data.frame(meta_filt)
  meta_linda$group <- factor(as.character(meta_linda$group), levels = c("control", "case"))
  
  out_linda <- local({
    result <- NULL
    suppressMessages(suppressWarnings(invisible(capture.output(
      result <- do.call(
        MicrobiomeStat::linda,
        c(
          list(
            feature.dat = otu_mat,
            meta.dat = meta_linda,
            formula = formula_chr
          ),
          method_args
        )
      ),
      file = nullfile()
    ))))
    result
  })
  
  res_tab <- out_linda$output[["groupcase"]]
  
  feature_stats <- data.frame(
    taxon_name = rownames(res_tab),
    pval = res_tab[["pvalue"]],
    qval = res_tab[["padj"]],
    called_linda = as.integer(res_tab[["reject"]]),
    baseMean = res_tab[["baseMean"]],
    log2FoldChange = res_tab[["log2FoldChange"]],
    lfcSE = res_tab[["lfcSE"]],
    stat = res_tab[["stat"]],
    df = res_tab[["df"]],
    stringsAsFactors = FALSE
  )
  rownames(feature_stats) <- NULL
  
  list(
    feature_stats = feature_stats,
    analysis_formula = formula_chr
  )
}

linda_adapter <- list(
  method_name = "LinDA",
  required_packages = c("MicrobiomeStat"),
  run_method = run_linda_method
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
      adapter = linda_adapter,
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
        is.winsor = TRUE,
        outlier.pct = 0.03,
        adaptive = TRUE,
        p.adj.method = "BH",
        alpha = 0.05,
        n.cores = 1,
        verbose = FALSE
      ),
      out_dir = "linda_25reg_benchmark_output"
    ),
    file = nullfile()
  ))))
}

if (sys.nframe() == 0) {
  main()
}
