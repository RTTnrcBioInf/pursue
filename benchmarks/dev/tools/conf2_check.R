#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# conf2_check.R -- pre-launch check of conf2's new real data (protocol 13, 2026-10-21). Loads every axisD2 / axisE2 /
# axisRM template through the engine's own loaders and reports format and design only: dimensions, depths, whether the
# counts are raw integers, the group sizes of the contrasts hpc/run_DE.sh SET=conf2 will run, repeated samples per
# subject, and the taxa two GEMS1 sites share. No method runs and no taxon is compared between groups. The one outcome
# summary is Vandeputte's measured load by group: the paper reports lower load in Crohn's disease, so it checks that the
# labels and the flow-cytometry counts line up. Repeated-measures templates are checked by rm_suite.R --dry.
#   Rscript benchmarks/dev/tools/conf2_check.R            (hpc/run_conf2_check.sh runs it and pushes the report)
# Exit status 1 if any template fails to load or a contrast lacks samples.
# -----------------------------------------------------------------------------
root <- normalizePath(Sys.getenv("PURSUE_BENCH_ROOT", unset = "benchmarks"), mustWork = TRUE)
Sys.setenv(PURSUE_BENCH_ROOT = root)
if (!nzchar(Sys.getenv("PURSUE_DATA_ROOT"))) Sys.setenv(PURSUE_DATA_ROOT = file.path(root, "data"))
source(file.path(root, "R", "engine", "templates.R"))
reg <- read_template_registry(); reg <- reg[reg$pool %in% c("axisD2", "axisE2", "axisRM"), , drop = FALSE]
LOADER <- c(phyloseq_rda = "load_phyloseq_rda", microbiomehd_tar = "load_microbiomehd_tar", gems1_txt = "load_gems1_txt",
            nearing_tar = "load_nearing_tar", biom_morton = "load_biom_morton", vandeputte_bundle = "load_vandeputte",
            template_rds = "load_template_rds")
# the contrasts of hpc/run_DE.sh SET=conf2 (control level first); NULL = the first two levels in sort order
LEVELS <- list(d2_vandeputte = c("healthy", "CD"), e2_cdi_schubert = c("H", "CDI"), e2_par_scheperjans = c("H", "PAR"),
               e2_t1d_alkanani = c("H", "T1D"), e2_asd_dan2020 = c("Control", "ASD"), e2_arctic_freshwater = NULL,
               e2_gems1_bangladesh = c("Control", "Case"), e2_gems1_kenya = c("Control", "Case"),
               e2_gems1_mali = c("Control", "Case"), e2_gems1_gambia = c("Control", "Case"))
cat(sprintf("conf2 data check -- %s\nR %s; rhdf5 %s\n\n", format(Sys.time(), "%F %T"), getRversion(),
            if (requireNamespace("rhdf5", quietly = TRUE)) as.character(utils::packageVersion("rhdf5")) else "MISSING (needed for e2_asd_dan2020)"))
q3 <- function(x) paste(format(stats::quantile(x, c(0, 0.5, 1), names = FALSE), big.mark = ",", scientific = FALSE, trim = TRUE), collapse = " / ")
fail <- character(0); kept <- list()

for (k in seq_len(nrow(reg))) {
  row <- reg[k, , drop = FALSE]; id <- row$id
  cat(sprintf("== %s  [%s, %s]\n", id, row$pool, row$format))
  r <- tryCatch({
    raw <- get(LOADER[[row$format]])(file.path(data_root(), row$path))
    rc <- raw$counts; int <- all(abs(rc - round(rc)) < 1e-8, na.rm = TRUE)
    cat(sprintf("   loader: %d features x %d samples with metadata; counts %s%s%s\n",
                nrow(rc), ncol(rc), if (int) "integer (raw)" else "NOT INTEGER (normalised?)", if (any(rc < 0, na.rm = TRUE)) "; NEGATIVE values" else "",
                if (anyNA(rc)) "; NA present" else ""))
    if (!is.na(row$subset) && nzchar(row$subset)) {
      s <- eval(parse(text = row$subset), envir = raw$meta); cat(sprintf("   subset %s: %d of %d samples\n", row$subset, sum(s %in% TRUE), length(s))) }
    tp <- load_template(id); ct <- tp$counts; N <- colSums(ct)
    cat(sprintf("   template: %d taxa (prevalence >= %s) x %d samples (depth >= %s); depth min / median / max %s%s\n",
                nrow(ct), row$min_prevalence, ncol(ct), row$min_depth, q3(N), if (length(unique(N)) == 1L) "  (all equal: rarefied)" else ""))
    m <- tp$meta; cat("   metadata columns:", paste(utils::head(names(m), 30), collapse = ", "), if (ncol(m) > 30) "..." else "", "\n")
    sc <- grep("^(subject|participant|host_subject_id|anonymized_name|individual|patient|host_id|subject_id)$", names(m), value = TRUE, ignore.case = TRUE)
    for (c in sc) cat(sprintf("   %s: %d distinct for %d samples%s\n", c, length(unique(m[[c]])), nrow(m),
                              if (anyDuplicated(m[[c]])) " -- REPEATED samples per subject" else " (one sample each)"))
    gv <- row$group_var
    if (!is.na(gv) && id %in% names(LEVELS)) {
      v <- as.character(m[[gv]]); tb <- table(v, useNA = "ifany"); names(tb)[is.na(names(tb))] <- "<NA>"; names(tb)[names(tb) == ""] <- "<blank>"
      cat("   ", gv, ": ", paste(names(tb), tb, sep = " = ", collapse = "; "), "\n", sep = "")
      lv <- LEVELS[[id]]; if (is.null(lv)) lv <- sort(unique(v[!is.na(v)]))[1:2]
      lv <- vapply(lv, function(l) { if (l %in% v) return(l); h <- unique(v[!is.na(v) & tolower(v) == tolower(l)]); if (length(h) == 1L) h else l }, "")
      nl <- vapply(lv, function(l) sum(v == l, na.rm = TRUE), 1L)
      ok <- all(nl >= 5L); if (!ok) fail <- c(fail, id)
      cat(sprintf("   contrast %s vs %s: %d vs %d samples%s%s\n", lv[2], lv[1], nl[2], nl[1],
                  if (row$pool == "axisE2") sprintf("; a 50/50 split leaves ~%d vs %d per half", nl[2] %/% 2L, nl[1] %/% 2L) else "",
                  if (ok) "" else "  -- TOO FEW (run_axisDE needs >= 5 per level)"))
      if (id == "d2_vandeputte") {
        cc <- m$cell_count; cat(sprintf("   cell_count: %d NA; median load %s (cells / g) -- %s\n", sum(is.na(cc)),
          paste(sprintf("%s %.3g", lv, vapply(lv, function(l) stats::median(cc[v == l], na.rm = TRUE), 1)), collapse = ", "),
          "the paper reports lower load in Crohn's disease")) }
      if (!is.null(m$age_months)) cat(sprintf("   age (months) min / median / max %s\n", q3(m$age_months[!is.na(m$age_months)])))
    }
    if (row$pool == "axisRM") cat("   (subjects and visits: rm_suite.R --dry)\n")
    kept[[id]] <- rownames(ct); TRUE
  }, error = function(e) { cat("   LOAD FAILED:", conditionMessage(e), "\n"); FALSE })
  if (!isTRUE(r)) fail <- c(fail, id)
  cat("\n")
}
g <- grep("^e2_gems1_", names(kept), value = TRUE)
if (length(g) > 1L) { cat("GEMS1 cross-site pairs (run_DE.sh crossstudy): taxa shared after each site's prevalence filter\n")
  for (i in seq_along(g)[-length(g)]) for (j in (i + 1L):length(g))
    cat(sprintf("   %s / %s: %d / %d taxa, %d shared\n", g[i], g[j], length(kept[[g[i]]]), length(kept[[g[j]]]), length(intersect(kept[[g[i]]], kept[[g[j]]])))) }
cat("\n", if (length(fail)) paste("conf2 data check FAILED:", paste(unique(fail), collapse = ", ")) else "conf2 data check: every template loads and every contrast has samples", "\n", sep = "")
quit(save = "no", status = as.integer(length(fail) > 0))
