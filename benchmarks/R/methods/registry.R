# -----------------------------------------------------------------------------
# registry.R -- the method roster (protocol section 6). `available` is checked at
# run time from the installed packages; a method that is not installed is recorded
# in the results with status "not_installed" instead of crashing the cell.
# -----------------------------------------------------------------------------

method_registry <- function() {
  data.frame(stringsAsFactors = FALSE, rbind(
    c("wilcoxon_tss",      "Wilcoxon on TSS",                "relative",  "single",   "",                "method_wilcoxon_tss"),
    c("lm_logtss",         "LM on log(TSS+pc)",              "relative",  "single",   "",                "method_lm_logtss"),
    c("limma_logtss",      "limma on log(TSS+pc)",           "relative",  "single",   "limma",           "method_limma_logtss"),
    c("logistic_presence", "Logistic on presence",           "prevalence","single",   "",                "method_logistic_presence"),
    c("pursue",            "PURSUE 0.2",                     "absolute",  "two-part", "PURSUE",          "method_pursue"),
    c("linda",             "LinDA",                          "absolute",  "single",   "MicrobiomeStat",  "method_linda"),
    c("ancombc2",          "ANCOM-BC2",                      "absolute",  "single",   "ANCOMBC",         "method_ancombc2"),
    c("maaslin3",          "MaAsLin 3",                      "absolute",  "two-part", "maaslin3",        "method_maaslin3"),
    c("aldex2",            "ALDEx2",                         "relative",  "single",   "ALDEx2",          "method_aldex2"),
    c("corncob",           "corncob",                        "relative",  "single",   "corncob",         "method_corncob"),
    c("ldm",               "LDM",                            "relative",  "single",   "LDM",             "method_ldm"),
    c("locom",             "LOCOM",                          "relative",  "single",   "LOCOM",           "method_locom"),
    c("locom2",            "LOCOM2",                         "relative",  "single",   "LOCOM2",          "method_locom2"),
    c("zicoseq",           "ZicoSeq",                        "relative",  "single",   "GUniFrac",        "method_zicoseq"),
    c("fastancom",         "fastANCOM",                      "absolute",  "single",   "fastANCOM",       "method_fastancom"),
    c("adapt",             "ADAPT",                          "absolute",  "single",   "ADAPT",           "method_adapt"),
    c("fastemu",           "fastEmu / radEmu",               "absolute",  "single",   "fastEmu",         "method_fastemu"),
    c("pursue03_erdl",     "PURSUE 0.3 cand.: ERD+ERL max",  "absolute",  "single",   "",                "method_pursue03_erdl"),
    c("pursue03_erdc",     "PURSUE 0.3 cand.: ERD centred",  "absolute",  "single",   "",                "method_pursue03_erdc"),
    c("pursue03_erlc",     "PURSUE 0.3 cand.: ERL centred",  "absolute",  "single",   "",                "method_pursue03_erlc"),
    c("pursue03_erdlu",    "PURSUE 0.3 cand.: pairwise max", "absolute",  "single",   "",                "method_pursue03_erdlu"),
    c("pursue03_erdu",     "PURSUE 0.3 cand.: pairwise ERD", "absolute",  "single",   "",                "method_pursue03_erdu"),
    c("pursue03_erlu",     "PURSUE 0.3 cand.: pairwise ERL", "absolute",  "single",   "",                "method_pursue03_erlu")))
}
.mr <- method_registry(); names(.mr) <- c("id", "label", "estimand", "arms", "package", "fn")
method_registry <- function() .mr

method_available <- function(id) {
  r <- .mr[.mr$id == id, ]; if (nrow(r) != 1L) return(FALSE)
  pk <- r$package; nzchar(pk) == FALSE || requireNamespace(pk, quietly = TRUE)
}

run_method <- function(id, counts, meta, formula, tested_term, args = list(), timeout_s = 3600) {
  r <- .mr[.mr$id == id, ]; feats <- rownames(counts)
  if (nrow(r) != 1L) stop("unknown method id: ", id)
  if (!method_available(id)) return(list(result = .empty_result(feats, status = "not_installed"), runtime_s = 0, mem_mb = NA, version = NA))
  fn <- get(r$fn, mode = "function")
  ver <- if (nzchar(r$package)) as.character(utils::packageVersion(r$package)) else as.character(getRversion())
  gc(); t0 <- proc.time()[["elapsed"]]
  res <- tryCatch(R.utils_withTimeout(fn(counts, meta, formula, tested_term, args), timeout_s),
                  error = function(e) { d <- .empty_result(feats, status = paste0("error: ", substr(conditionMessage(e), 1, 120))); d })
  rt <- proc.time()[["elapsed"]] - t0
  mem <- tryCatch(sum(gc()[, 6]) , error = function(e) NA_real_)
  list(result = res, runtime_s = rt, mem_mb = mem, version = ver)
}

R.utils_withTimeout <- function(expr, timeout_s) {
  setTimeLimit(elapsed = timeout_s, transient = TRUE); on.exit(setTimeLimit(elapsed = Inf))
  expr
}
