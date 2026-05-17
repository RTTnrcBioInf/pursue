
suppressPackageStartupMessages({
  suppressMessages(suppressWarnings(invisible(capture.output(
    source("benchmark_engine.R"),
    file = nullfile()
  ))))
  library(LOCOM)
})

run_locom_method <- function(otu_full,
                             otu_filt,
                             meta_full,
                             meta_filt,
                             regime_row,
                             alpha,
                             tested_metadata,
                             method_args) {
  otu_locom <- t(otu_filt)
  storage.mode(otu_locom) <- "numeric"
  
  meta_locom <- as.data.frame(meta_filt)
  meta_locom$group <- factor(as.character(meta_locom$group), levels = c("control", "case"))
  
  call_args <- c(
    list(
      otu.table = otu_locom,
      Y = ifelse(meta_locom$group == "case", 1, 0),
      seed = sample.int(.Machine$integer.max, 1)
    ),
    method_args
  )
  
  analysis_formula <- "LOCOM: Y = case indicator"
  
  if (regime_row$confounder.type != "none") {
    call_args$C <- meta_locom$confounder
    analysis_formula <- "LOCOM: Y = case indicator; C = confounder"
  }
  
  out_locom <- local({
    result <- NULL
    suppressMessages(suppressWarnings(invisible(capture.output(
      result <- do.call(LOCOM::locom, call_args),
      file = nullfile()
    ))))
    result
  })
  
  result_taxa <- colnames(out_locom[["p.otu"]])
  p_vec <- as.numeric(out_locom[["p.otu"]][1, result_taxa])
  q_vec <- as.numeric(out_locom[["q.otu"]][1, result_taxa])
  
  feature_stats <- data.frame(
    taxon_name = result_taxa,
    pval = p_vec,
    qval = q_vec,
    called_native = as.integer(result_taxa %in% out_locom[["detected.otu"]]),
    effect_size = as.numeric(out_locom[["effect.size"]])[seq_along(result_taxa)],
    stringsAsFactors = FALSE
  )
  rownames(feature_stats) <- NULL
  
  list(
    feature_stats = feature_stats,
    analysis_formula = analysis_formula
  )
}

locom_adapter <- list(
  method_name = "LOCOM",
  required_packages = c("LOCOM"),
  run_method = run_locom_method
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
      adapter = locom_adapter,
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
      method_args = list(
        n.cores = 1,
        n.perm.max = 1000,
        filter.thresh = 0,
        fdr.nominal = 0.05,
        verbose = FALSE
      ),
      out_dir = "locom_25reg_benchmark_output"
    ),
    file = nullfile()
  ))))
}

if (sys.nframe() == 0) {
  main()
}
