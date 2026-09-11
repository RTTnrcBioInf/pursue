# -----------------------------------------------------------------------------
# io.R -- the results contract (protocol section 12) and run manifests.
# One gzipped CSV per (axis, simulator, template, regime, replicate) with the
# feature-level long table; one JSON manifest per cell; cell-level metrics as a
# second CSV. Parquet is used instead of CSV when the `arrow` package is present.
# -----------------------------------------------------------------------------

contract_columns <- c("axis", "simulator", "template", "regime_id", "replicate", "seed", "method", "method_version",
                      "arm", "feature", "truth_abs", "truth_rel", "truth_type", "truth_effect", "p", "q",
                      "estimate", "se", "ci_lo", "ci_hi", "status", "runtime_s", "mem_mb")

assemble_contract <- function(cell, method_id, run, truth, version) {
  r <- run$result; tr <- truth[match(r$feature, truth$feature), ]
  te <- if ("truth_lfc2" %in% names(tr)) ifelse(tr$truth_type %in% c("prevalence"), if ("truth_logor" %in% names(tr)) tr$truth_logor else NA, tr$truth_lfc2) else NA
  data.frame(axis = cell$axis, simulator = cell$simulator, template = cell$template, regime_id = cell$regime_id,
             replicate = cell$replicate, seed = cell$seed, method = method_id, method_version = version,
             arm = r$arm, feature = r$feature, truth_abs = tr$truth_abs, truth_rel = tr$truth_rel, truth_type = tr$truth_type,
             truth_effect = te, p = r$p, q = r$q, estimate = r$estimate, se = r$se, ci_lo = r$ci_lo, ci_hi = r$ci_hi,
             status = r$status, runtime_s = run$runtime_s, mem_mb = run$mem_mb, stringsAsFactors = FALSE)[, contract_columns]
}

cell_id <- function(cell) sprintf("%s__%s__%s__%s__r%03d", cell$axis, cell$simulator, cell$template, cell$regime_id, cell$replicate)

write_cell <- function(cell, contract, metrics, extra = list(), out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  id <- cell_id(cell)
  if (requireNamespace("arrow", quietly = TRUE)) {
    arrow::write_parquet(contract, file.path(out_dir, paste0(id, ".features.parquet")))
  } else {
    con <- gzfile(file.path(out_dir, paste0(id, ".features.csv.gz")), "w"); utils::write.csv(contract, con, row.names = FALSE); close(con)
  }
  utils::write.csv(metrics, file.path(out_dir, paste0(id, ".metrics.csv")), row.names = FALSE)
  manifest <- c(list(cell = cell, written = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), r_version = R.version.string,
                     protocol_version = "0.9", git_commit = tryCatch(system("git rev-parse HEAD", intern = TRUE, ignore.stderr = TRUE), error = function(e) NA),
                     hostname = Sys.info()[["nodename"]]), extra)
  writeLines(jsonlite::toJSON(manifest, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null"), file.path(out_dir, paste0(id, ".manifest.json")))
  invisible(id)
}

read_all_metrics <- function(out_dir) {
  fs <- list.files(out_dir, pattern = "\\.metrics\\.csv$", full.names = TRUE, recursive = TRUE)
  do.call(rbind, lapply(fs, utils::read.csv, stringsAsFactors = FALSE))
}
