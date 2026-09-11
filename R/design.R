# -----------------------------------------------------------------------------
# design.R -- model formula handling shared by both arms
# -----------------------------------------------------------------------------

#' Build the design used by both arms
#'
#' Expands `formula` on `meta`, identifies the columns that belong to `tested_term`,
#' and returns everything the arms need. The abundance arm appends a centred
#' log-depth column; the presence arm uses the depth as an offset instead.
#'
#' @param formula RHS-only formula, e.g. `~ group + age`.
#' @param meta sample metadata (data.frame), rows aligned to the count table.
#' @param tested_term single term label present in `formula` (e.g. `"group"`).
#' @param depth numeric vector of per-sample library sizes.
#' @return list with `X` (model matrix, intercept first), `tested_cols` (integer
#'   indices of the tested term's columns in `X`), `nuisance_cols`, `log_depth`
#'   (centred), `term_labels`, and `coef_names`.
#' @keywords internal
build_design <- function(formula, meta, tested_term, depth) {
  if (!inherits(formula, "formula")) stop("`formula` must be a formula, e.g. ~ group + age.")
  if (length(formula) == 3L) formula <- formula[-2L]           # drop any LHS
  tt <- stats::terms(formula, data = meta)
  labels <- attr(tt, "term.labels")
  if (length(labels) == 0L) stop("`formula` must contain at least the tested term.")
  if (!tested_term %in% labels)
    stop("`tested_term` = '", tested_term, "' is not a term of `formula` (terms: ",
         paste(labels, collapse = ", "), ").")
  needed <- all.vars(formula)
  missing_vars <- setdiff(needed, colnames(meta))
  if (length(missing_vars)) stop("Variables not in `meta`: ", paste(missing_vars, collapse = ", "))
  mf <- stats::model.frame(tt, data = meta, na.action = stats::na.pass)
  X  <- stats::model.matrix(tt, mf)
  if (anyNA(X)) stop("Metadata used in `formula` contain NA; remove those samples first.")
  assign <- attr(X, "assign")
  tested_idx <- which(labels == tested_term)
  tested_cols <- which(assign == tested_idx)
  if (length(tested_cols) == 0L) stop("The tested term contributes no columns (constant variable?).")
  nuisance_cols <- setdiff(seq_len(ncol(X)), tested_cols)
  if (qr(X)$rank < ncol(X)) stop("Design matrix is rank deficient; check for collinear covariates.")
  ld <- log(depth); ld <- ld - mean(ld)
  list(X = X, tested_cols = tested_cols, nuisance_cols = nuisance_cols,
       log_depth = ld, term_labels = labels, coef_names = colnames(X),
       tested_term = tested_term, n_tested = length(tested_cols))
}

#' Coerce a count table to a numeric samples x features matrix
#' @keywords internal
as_count_matrix <- function(otu) {
  if (is.data.frame(otu)) otu <- as.matrix(otu)
  if (!is.numeric(otu)) stop("`otu` must be numeric (a count matrix or data.frame).")
  storage.mode(otu) <- "double"
  if (any(!is.finite(otu))) stop("`otu` contains non-finite values.")
  if (any(otu < 0)) stop("`otu` contains negative values.")
  if (is.null(colnames(otu))) colnames(otu) <- paste0("feature_", seq_len(ncol(otu)))
  if (is.null(rownames(otu))) rownames(otu) <- paste0("sample_", seq_len(nrow(otu)))
  otu
}

#' Align metadata rows to the count table
#' @keywords internal
align_meta <- function(otu, meta) {
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  if (nrow(meta) != nrow(otu)) stop("`meta` must have one row per sample (row of `otu`).")
  if (!is.null(rownames(meta)) && !identical(rownames(meta), rownames(otu))) {
    if (setequal(rownames(meta), rownames(otu))) meta <- meta[rownames(otu), , drop = FALSE]
    else if (all(rownames(meta) == as.character(seq_len(nrow(meta))))) rownames(meta) <- rownames(otu)
    else stop("Row names of `meta` do not match row names of `otu`.")
  }
  meta
}
