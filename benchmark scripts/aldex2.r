
invisible(capture.output(suppressMessages(suppressWarnings(source("benchmark_engine.R")))))
invisible(capture.output(suppressPackageStartupMessages(library(ALDEx2))))

run_aldex2_method <- function(otu_full,
                              otu_filt,
                              meta_full,
                              meta_filt,
                              regime_row,
                              alpha,
                              tested_metadata,
                              method_args) {
  otu_aldex <- as.data.frame(otu_filt, check.names = FALSE)
  otu_aldex[] <- lapply(otu_aldex, function(x) {
    if (is.factor(x)) x <- as.character(x)
    as.numeric(x)
  })
  otu_aldex <- as.matrix(otu_aldex)
  rownames(otu_aldex) <- rownames(otu_filt)
  meta_aldex <- as.data.frame(meta_filt)[colnames(otu_aldex), , drop = FALSE]
  meta_aldex$group <- factor(meta_aldex$group, levels = c("control", "case"))
  
  if (regime_row$confounder.type == "none") {
    conds <- ifelse(as.character(meta_aldex$group) == "case", "case", "control")
    
    x.clr <- ALDEx2::aldex.clr(
      otu_aldex,
      conds,
      mc.samples = 128,
      denom = "all",
      verbose = FALSE
    )
    
    x.tt <- ALDEx2::aldex.ttest(
      x.clr,
      paired.test = FALSE,
      hist.plot = FALSE,
      verbose = FALSE
    )
    
    stopifnot(all(c("wi.ep", "wi.eBH", "we.ep", "we.eBH") %in% colnames(x.tt)))
    
    feature_stats <- data.frame(
      taxon_name = rownames(x.tt),
      pval = x.tt[["wi.ep"]],
      qval = x.tt[["wi.eBH"]],
      we_pval = x.tt[["we.ep"]],
      we_qval = x.tt[["we.eBH"]],
      stringsAsFactors = FALSE
    )
    
    analysis_formula <- "ALDEx2 Wilcoxon on group"
  } else {
    mm <- model.matrix(~ group + confounder, data = meta_aldex)
    
    x.clr <- ALDEx2::aldex.clr(
      otu_aldex,
      mm,
      mc.samples = 128,
      denom = "all",
      verbose = FALSE
    )
    
    x.glm <- ALDEx2::aldex.glm(
      x.clr,
      verbose = FALSE,
      fdr.method = "BH"
    )
    
    stopifnot("groupcase" %in% colnames(mm))
    
    p_col <- "groupcase:pval"
    q_col <- "groupcase:pval.padj"
    
    if (!all(c(p_col, q_col) %in% colnames(x.glm))) {
      stop(
        "ALDEx2 GLM output did not contain the expected exact group columns. ",
        "Expected columns: ", paste(c(p_col, q_col), collapse = ", "), ". ",
        "Observed columns: ", paste(colnames(x.glm), collapse = ", ")
      )
    }
    
    feature_stats <- data.frame(
      taxon_name = rownames(x.glm),
      pval = x.glm[[p_col]],
      qval = x.glm[[q_col]],
      stringsAsFactors = FALSE
    )
    
    analysis_formula <- "ALDEx2 GLM on group + confounder"
  }
  
  list(feature_stats = feature_stats, analysis_formula = analysis_formula)
}

aldex2_adapter <- list(
  method_name = "ALDEx2",
  required_packages = "ALDEx2",
  run_method = run_aldex2_method
)

main <- function() {
  res <- NULL
  invisible(capture.output(
    res <- suppressMessages(suppressWarnings(run_method_benchmark(
      adapter = aldex2_adapter,
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
      method_args = list(),
      out_dir = "aldex2_25reg_benchmark_output"
    )))
  ))
  invisible(res)
}

if (sys.nframe() == 0) {
  invisible(main())
}
