# -----------------------------------------------------------------------------
# templates.R -- load real 16S count tables used as simulation templates and as
# real-data axes. Every loader returns list(counts = integer matrix FEATURES x SAMPLES,
# meta = data.frame with rownames = sample ids, id = template id).
#
# Registry: benchmarks/templates.tsv with columns
#   id, source, format, path, subset, min_prevalence, min_depth, group_var, notes
# format in {phyloseq_rda, tsv_features_x_samples, csv_samples_x_features}
# subset: an R expression on the metadata evaluated with `with(meta, ...)`, or ""
# -----------------------------------------------------------------------------

read_template_registry <- function(path = file.path(bench_root(), "templates.tsv")) {
  utils::read.delim(path, stringsAsFactors = FALSE, comment.char = "#", quote = "")
}

bench_root <- function() {
  r <- Sys.getenv("PURSUE_BENCH_ROOT", unset = "")
  if (nzchar(r)) return(r)
  # fall back: directory containing this file's grandparent (benchmarks/R/engine -> benchmarks)
  normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), "..", ".."), mustWork = FALSE)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

data_root <- function() {
  r <- Sys.getenv("PURSUE_DATA_ROOT", unset = "")
  if (nzchar(r)) r else file.path(bench_root(), "data")
}

#' Load a template by id (or a registry row)
load_template <- function(id, registry = read_template_registry(), min_prevalence = NULL,
                          min_depth = NULL, max_samples = NULL, seed = 1) {
  row <- if (is.data.frame(id)) id else registry[registry$id == id, , drop = FALSE]
  if (nrow(row) != 1L) stop("template id not found or not unique: ", id)
  path <- file.path(data_root(), row$path)
  if (!file.exists(path)) stop("template file missing: ", path, "\n  run hpc/download_templates.sh first")
  obj <- switch(row$format,
    phyloseq_rda = load_phyloseq_rda(path),
    tsv_features_x_samples = load_tsv_fxs(path, row),
    csv_samples_x_features = load_csv_sxf(path, row),
    stop("unknown template format: ", row$format))
  counts <- obj$counts; meta <- obj$meta
  # subset
  if (!is.na(row$subset) && nzchar(row$subset)) {
    keep <- eval(parse(text = row$subset), envir = meta)
    keep[is.na(keep)] <- FALSE
    counts <- counts[, keep, drop = FALSE]; meta <- meta[keep, , drop = FALSE]
  }
  # depth and prevalence filters
  md <- if (!is.null(min_depth)) min_depth else if (!is.na(row$min_depth)) row$min_depth else 1000
  keep_s <- colSums(counts) >= md
  counts <- counts[, keep_s, drop = FALSE]; meta <- meta[keep_s, , drop = FALSE]
  mp <- if (!is.null(min_prevalence)) min_prevalence else if (!is.na(row$min_prevalence)) row$min_prevalence else 0.05
  keep_f <- rowMeans(counts > 0) >= mp
  counts <- counts[keep_f, , drop = FALSE]
  if (!is.null(max_samples) && ncol(counts) > max_samples) {
    set.seed(seed); s <- sort(sample.int(ncol(counts), max_samples))
    counts <- counts[, s, drop = FALSE]; meta <- meta[s, , drop = FALSE]
  }
  storage.mode(counts) <- "integer"
  # sps_group: the metadata variable SPsimSeq borrows real between-group differences from
  # (see simulate_sps). Falls back to group_var; NA means sps cannot use this template.
  sg <- if (!is.null(row$sps_group) && !is.na(row$sps_group) && nzchar(row$sps_group)) row$sps_group else
        if (is.na(row$group_var)) NA_character_ else row$group_var
  list(counts = counts, meta = meta, id = row$id,
       group_var = if (is.na(row$group_var)) NA else row$group_var, sps_group = sg)
}

load_phyloseq_rda <- function(path) {
  if (!requireNamespace("phyloseq", quietly = TRUE)) stop("phyloseq is required to read ", path)
  e <- new.env(); nm <- load(path, envir = e); ps <- get(nm[1], envir = e)
  ot <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) ot <- t(ot)
  sd <- tryCatch(data.frame(phyloseq::sample_data(ps), check.names = FALSE, stringsAsFactors = FALSE),
                 error = function(e) data.frame(row.names = colnames(ot)))
  sd <- sd[colnames(ot), , drop = FALSE]
  list(counts = ot, meta = sd)
}

load_tsv_fxs <- function(path, row) {
  ct <- as.matrix(utils::read.delim(path, row.names = 1, check.names = FALSE, quote = ""))
  meta_path <- file.path(data_root(), sub("_count_matrix\\.tsv$", "_sample_metadata.tsv", row$path))
  meta <- if (file.exists(meta_path)) utils::read.delim(meta_path, row.names = 1, check.names = FALSE, quote = "", stringsAsFactors = FALSE)
          else data.frame(row.names = colnames(ct))
  common <- intersect(colnames(ct), rownames(meta))
  list(counts = ct[, common, drop = FALSE], meta = meta[common, , drop = FALSE])
}

load_csv_sxf <- function(path, row) {
  ct <- utils::read.csv(path, row.names = 1, check.names = FALSE)
  ct <- t(as.matrix(ct))
  meta_path <- file.path(data_root(), dirname(row$path), "metadata.csv")
  meta <- if (file.exists(meta_path)) utils::read.csv(meta_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
          else data.frame(row.names = colnames(ct))
  common <- intersect(colnames(ct), rownames(meta))
  list(counts = ct[, common, drop = FALSE], meta = meta[common, , drop = FALSE])
}

#' Per-feature template summaries used by the in-house simulator and the realism gate
template_profile <- function(tpl) {
  ct <- tpl$counts; N <- colSums(ct); R <- sweep(ct, 2L, N, "/")
  prev <- rowMeans(ct > 0)
  logr <- log(R); logr[!is.finite(logr)] <- NA
  list(prevalence = prev,
       mean_log_relab = apply(logr, 1L, mean, na.rm = TRUE),
       sd_log_relab = apply(logr, 1L, stats::sd, na.rm = TRUE),
       depth = N, n_features = nrow(ct), n_samples = ncol(ct))
}
