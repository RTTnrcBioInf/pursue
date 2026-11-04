# msq_cache.R -- SimulateMSeq without refitting its template for every cell (2026-11-04). Development drivers only
# (devsuite.R, bank_cells.R); the benchmark engine is not changed.
#
# GUniFrac's EstPara() is called once per msq cell (benchmarks/R/simulators/dispatch.R) and spends nearly all its time in
# dirmult::dirmult(), a Newton fit of a Dirichlet-multinomial to the whole template: about 25 minutes per twinsuk_stool cell
# at single-threaded BLAS, a minute for hmp_tongue. That fit is deterministic (no random numbers) and its input is the
# template, not the cell. The rest of EstPara (Dirichlet draws, size factors) does use random numbers, from the cell's own
# seed, so caching EstPara's result would change every cell; caching only dirmult() does not. This memoises dirmult() on
# disk under <cache_dir>/dirmult/, keyed by the md5 of its serialised input, so the fit runs once per template (per R
# version) and every cell stays bit-identical. hpc/run_srv35.sh checks that against cells banked without the cache.
# Off with PURSUE_MSQ_CACHE=0.
.msq_cache_on <- function(cache_dir) {
  if (identical(Sys.getenv("PURSUE_MSQ_CACHE"), "0") || !requireNamespace("dirmult", quietly = TRUE)) return(invisible(FALSE))
  ns <- asNamespace("dirmult"); orig <- get("dirmult", envir = ns)
  if (isTRUE(attr(orig, "pursue_memo"))) return(invisible(TRUE))
  dir <- file.path(cache_dir, "dirmult"); dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  memo <- function(data, init, initscalar, epsilon = 10^(-4), trace = TRUE, mode) {
    if (!missing(init) || !missing(initscalar) || !missing(mode) || !identical(epsilon, 10^(-4)))   # only EstPara's call is cached
      return(orig(data, init, initscalar, epsilon, trace, mode))
    tf <- tempfile(fileext = ".rds"); on.exit(unlink(tf), add = TRUE); saveRDS(data, tf, compress = FALSE)
    f <- file.path(dir, paste0(unname(tools::md5sum(tf)), ".rds"))
    if (file.exists(f)) { v <- tryCatch(readRDS(f), error = function(e) NULL); if (!is.null(v)) return(v) }
    v <- orig(data, trace = trace); tmp <- paste0(f, ".tmp", Sys.getpid()); saveRDS(v, tmp); file.rename(tmp, f); v
  }
  attr(memo, "pursue_memo") <- TRUE
  utils::assignInNamespace("dirmult", memo, ns = "dirmult")
  invisible(TRUE)
}
# Fit each template once, in parallel, before the workers fork -- otherwise every worker that starts on the same template
# fits it at the same time. The input is built exactly as simulate_msq() and EstPara() build it.
.msq_cache_warm <- function(TPL, cores = 1L) {
  if (!isTRUE(attr(get("dirmult", envir = asNamespace("dirmult")), "pursue_memo"))) return(invisible(NULL))
  t0 <- proc.time()[["elapsed"]]
  invisible(parallel::mclapply(TPL, function(tpl) { ref <- tpl$counts; storage.mode(ref) <- "numeric"
    if (is.null(rownames(ref))) rownames(ref) <- paste0("OTU", seq_len(nrow(ref)))
    dirmult::dirmult(t(ref), trace = FALSE); NULL }, mc.cores = max(1L, min(cores, length(TPL)))))
  cat(sprintf("msq cache: %d template fit(s) ready (%.0f s)\n", length(TPL), proc.time()[["elapsed"]] - t0))
}
