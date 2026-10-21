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
  if (!file.exists(sub("::.*$", "", path))) stop("template file missing: ", path, "\n  run hpc/download_templates.sh first")
  obj <- switch(row$format,
    phyloseq_rda = load_phyloseq_rda(path),
    tsv_features_x_samples = load_tsv_fxs(path, row),
    csv_samples_x_features = load_csv_sxf(path, row),
    microbiomehd_tar = load_microbiomehd_tar(path),          # conf2 formats (protocol section 13, 2026-10-21)
    gems1_txt = load_gems1_txt(path),
    nearing_tar = load_nearing_tar(path),
    biom_morton = load_biom_morton(path),
    vandeputte_bundle = load_vandeputte(path),
    template_rds = load_template_rds(path),
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
  # conf2's real data (protocol 13): one sample per person -- repeat samples are dropped here, the first kept, and logged
  if (!is.null(row$pool) && row$pool %in% c("axisD2", "axisE2")) {
    sc <- names(meta)[tolower(names(meta)) %in% c("subject", "subject_id", "host_subject_id", "participant", "participant_id",
                                                  "patient", "patient_id", "individual", "anonymized_name", "host_id")]
    if (length(sc)) { s <- as.character(meta[[sc[1]]]); d <- duplicated(s) & !is.na(s) & nzchar(s)
      if (any(d)) { message(sprintf("template %s: %d repeat samples of %d people dropped (%s; first kept)", row$id, sum(d), length(unique(s[d])), sc[1]))
        counts <- counts[, !d, drop = FALSE]; meta <- meta[!d, , drop = FALSE] } }
  }
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

# ---- conf2 real-data formats (protocol section 13, 2026-10-21; data fetched by hpc/run_fetch_conf2.sh) ----
.tar_members <- function(tarfile, pattern) {            # extract the members whose path matches `pattern`; return local paths
  all <- utils::untar(tarfile, list = TRUE); hit <- all[grepl(pattern, all)]
  if (!length(hit)) stop("no member matching ", pattern, " in ", tarfile)
  td <- tempfile("tpl"); dir.create(td); utils::untar(tarfile, files = hit, exdir = td); file.path(td, hit)
}
.align_ids <- function(ct, meta) {                      # match count-table and metadata sample ids; normalise if few match
  norm <- function(x) gsub("[.]", "-", sub("^X(?=[0-9])", "", x, perl = TRUE))
  a <- colnames(ct); b <- rownames(meta)
  if (mean(a %in% b) < 0.5) { a2 <- norm(a); b2 <- norm(b); ib <- match(a2, b2); keep <- !is.na(ib)
    ct <- ct[, keep, drop = FALSE]; meta <- meta[ib[keep], , drop = FALSE]; rownames(meta) <- colnames(ct)
  } else { common <- intersect(a, b); ct <- ct[, common, drop = FALSE]; meta <- meta[common, , drop = FALSE] }
  list(counts = ct, meta = meta)
}
# MicrobiomeHD (Duvallet et al. 2017): RDP-assigned de novo OTU table, collapsed to genus by summing the OTUs that share
# kingdom..genus (as Pelto et al. 2025); metadata.txt with DiseaseState. DiseaseClass = DiseaseState up to the first
# "_" (t1d_alkanani's "T1D_new-onset" joins "T1D").
load_microbiomehd_tar <- function(path) {
  st <- sub("_results[.]tar[.]gz$", "", basename(path))
  f <- .tar_members(path, paste0("(RDP/", st, "[.]otu_table[.]100[.]denovo[.]rdp_assigned|/", st, "[.]metadata[.]txt)$"))
  ot <- utils::read.delim(f[grepl("rdp_assigned$", f)], row.names = 1, check.names = FALSE, quote = "", comment.char = "")
  lin <- vapply(strsplit(rownames(ot), ";", fixed = TRUE), function(z) paste(utils::head(z, 6), collapse = ";"), "")
  ct <- rowsum(as.matrix(ot), lin)
  meta <- utils::read.delim(f[grepl("metadata[.]txt$", f)], row.names = 1, check.names = FALSE, quote = "",
                            stringsAsFactors = FALSE, fileEncoding = "latin1")
  if ("DiseaseState" %in% names(meta)) meta$DiseaseClass <- sub("_.*$", "", meta$DiseaseState)
  .align_ids(ct, meta)
}
# GEMS1 (microbiomeData data-raw): genus-level read counts (EUPATH_0009356_*; the EUPATH_00092xx columns are proportions),
# samples in rows; participant table gives case / control and country. One sample per participant.
load_gems1_txt <- function(path) {
  rd <- function(f) utils::read.delim(file.path(path, f), check.names = FALSE, quote = "", stringsAsFactors = FALSE)
  a <- rd("GEMS_16S_rRNA_V1-V2_assay.txt"); p <- rd("GEMS_Participant.txt")
  g <- grepl("[[]EUPATH_0009356_", colnames(a)); if (!any(g)) stop("GEMS1: no genus-count columns (EUPATH_0009356)")
  ct <- t(as.matrix(a[, g, drop = FALSE])); ct[is.na(ct)] <- 0
  dimnames(ct) <- list(sub("^.*[[]EUPATH_0009356_(.*)[]]$", "\\1", colnames(a)[g]), a$Sample_Id)
  col <- function(pat) { k <- grep(pat, names(p), ignore.case = TRUE); if (!length(k)) stop("GEMS1: no column ", pat); p[[k[1]]] }
  i <- match(a$Participant_Id, p$Participant_Id)
  meta <- data.frame(row.names = a$Sample_Id, participant = a$Participant_Id, case_control = col("^Case or control")[i],
                     Country = col("^Country")[i], age_months = col("^Age [(]months")[i], stringsAsFactors = FALSE)
  list(counts = ct, meta = meta)
}
# Nearing et al. 2022 (figshare 14531724): "<archive>::<Study>"; the study's genus table (taxa x samples, an optional
# "# Constructed from biom file" first line) and its two-column metadata (sample, grouping -> `comparison`).
load_nearing_tar <- function(path) {
  tf <- sub("::.*$", "", path); ds <- sub("^.*::", "", path)
  f <- .tar_members(tf, paste0("Studies/", ds, "/", ds, "_(genus_table|meta)[.]tsv$"))
  gt <- f[grepl("_genus_table[.]tsv$", f)][1]; mf <- f[grepl("_meta[.]tsv$", f)][1]
  sk <- as.integer(grepl("Constructed from biom file", readLines(gt, n = 1L)))
  ct <- as.matrix(utils::read.delim(gt, skip = sk, row.names = 1, check.names = FALSE, quote = "", comment.char = ""))
  meta <- utils::read.delim(mf, row.names = 1, check.names = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE)
  names(meta)[1] <- "comparison"
  .align_ids(ct, meta)
}
# Morton et al. 2023 ASD cohorts: biom 2.1 (HDF5) table + the matching metadata file in the same folder.
load_biom_morton <- function(path) {
  if (!requireNamespace("rhdf5", quietly = TRUE)) stop("rhdf5 is required to read ", path)
  rd <- function(n) as.vector(rhdf5::h5read(path, n))
  obs <- rd("observation/ids"); smp <- rd("sample/ids"); x <- rd("observation/matrix/data")
  j <- rd("observation/matrix/indices"); ip <- rd("observation/matrix/indptr")
  ct <- matrix(0, length(obs), length(smp), dimnames = list(obs, smp)); ct[cbind(rep(seq_along(obs), diff(ip)), j + 1L)] <- x
  mf <- file.path(dirname(path), if (grepl("^sibling", basename(path))) "sibling_match_metadata.txt" else "age_sex_match_metadata.txt")
  meta <- utils::read.delim(mf, check.names = FALSE, quote = "", stringsAsFactors = FALSE)
  rownames(meta) <- meta$sampleid
  .align_ids(ct, meta)
}
# Vandeputte et al. 2017 disease cohort: deblurred sOTU table (DC01-DC95, Brill et al.'s bundle of the paper's data) and
# flow-cytometry cell counts by serial sample. Labels as in that bundle's analysis code: DC01-DC29 Crohn's disease,
# DC30-DC95 healthy (29 vs 66, as in the paper).
load_vandeputte <- function(path) {
  e <- new.env(); load(file.path(path, "gut_otu_table.RData"), envir = e); ot <- as.matrix(get("otu_table", envir = e))
  id <- sub("[._]R1$", "", colnames(ot)); serial <- suppressWarnings(as.integer(sub("^DC", "", id)))   # DC01_R1 .. DC29_R1, DC30.R1 ..
  o <- order(serial); ot <- ot[, o, drop = FALSE]; id <- id[o]; serial <- serial[o]
  flow <- utils::read.csv(file.path(path, "FlowCytometryCounts.csv"))$Average_Cell_Count
  if (anyNA(serial) || length(flow) != 95L || any(serial < 1 | serial > 95)) stop("Vandeputte bundle: unexpected sample ids or flow table")
  meta <- data.frame(row.names = colnames(ot), sample = id, group = ifelse(serial <= 29L, "CD", "healthy"),
                     cell_count = flow[serial], stringsAsFactors = FALSE)
  list(counts = ot, meta = meta)
}
# A saved template list(counts, meta, subject_var, time_var) (dev/tools/fetch_longitudinal.R): adds `subject` and
# `visit` (order of the time variable within subject).
load_template_rds <- function(path) {
  x <- readRDS(path); m <- x$meta
  if (!is.null(x$subject_var) && !is.na(x$subject_var)) {
    m$subject <- as.character(m[[x$subject_var]])
    tv <- if (!is.null(x$time_var) && !is.na(x$time_var)) m[[x$time_var]] else seq_len(nrow(m))
    tn <- tryCatch(as.numeric(as.Date(as.character(tv))), error = function(e) rep(NA_real_, length(tv)))
    if (all(is.na(tn))) tn <- suppressWarnings(as.numeric(tv)); if (all(is.na(tn))) tn <- seq_len(nrow(m))
    m$visit <- stats::ave(tn, m$subject, FUN = function(r) rank(r, ties.method = "first"))
  }
  list(counts = as.matrix(x$counts), meta = m)
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
