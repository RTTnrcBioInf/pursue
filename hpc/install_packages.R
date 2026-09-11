#!/usr/bin/env Rscript
# Install everything the benchmark can use, in order of importance and tolerating any failure.
#
#   PURSUE itself is installed FIRST and from this checkout, so it lands even on a node with
#   no outbound network (it needs only limma, sandwich, parallel). Everything after it is
#   best-effort: a package that will not install is reported, and the benchmark records that
#   method or simulator as `not_installed` rather than crashing a cell.
#
# Always writes hpc/installed_packages.csv, including when steps fail -- that table plus
# hpc/smoke_*.csv is what says which wrappers are verified on this cluster.
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = max(1L, parallel::detectCores() - 1L), timeout = 600)
# R_LIBS_USER may be a colon-separated list; installs go to the first entry, the rest stay readable.
libs <- Filter(nzchar, strsplit(Sys.getenv("R_LIBS_USER", unset = ""), .Platform$path.sep)[[1]])
if (length(libs)) { dir.create(libs[1], recursive = TRUE, showWarnings = FALSE); .libPaths(c(libs, .libPaths())) }
has <- function(p) requireNamespace(p, quietly = TRUE)
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), ".."), mustWork = FALSE)
# Every failure also goes to hpc/install_log.txt -- console output scrolls away on a cluster,
# and "package X did not install" is useless without the reason.
logf <- file.path(root, "hpc", "install_log.txt")
dir.create(dirname(logf), showWarnings = FALSE)
cat("install_packages.R", format(Sys.time()), "\n", R.version.string, "\n\n", file = logf)
say <- function(...) { msg <- paste0(...); message(msg); cat(msg, "\n", file = logf, append = TRUE) }
try_ <- function(what, expr) tryCatch({ expr; TRUE },
  error   = function(e) { say("  !! ", what, ": ", conditionMessage(e)); FALSE },
  warning = function(w) { say("  ~  ", what, ": ", conditionMessage(w)); has(what) })

cat(">> R", R.version.string, "\n>> library:", .libPaths()[1], "\n>> repo root:", root, "\n\n>> 1. PURSUE (from this checkout)\n")
if (!file.exists(file.path(root, "DESCRIPTION")))
  say("  !! no DESCRIPTION at ", root, " -- is this the repository root? PURSUE cannot install.")
# Run R CMD INSTALL directly and keep every line: "had non-zero exit status" on its own is
# useless, and PURSUE is the one package whose failure stops the whole benchmark.
pi_out <- suppressWarnings(system2(file.path(R.home("bin"), "R"), c("CMD", "INSTALL", shQuote(root)),
                                   stdout = TRUE, stderr = TRUE))
if (!is.null(attr(pi_out, "status")) && attr(pi_out, "status") != 0) {
  say("  !! PURSUE: R CMD INSTALL exited ", attr(pi_out, "status"), ". Full output:")
  cat(paste0("     ", pi_out, collapse = "\n"), "\n", file = logf, append = TRUE)
  cat(paste(utils::tail(pi_out, 25), collapse = "\n"), "\n")
  try_("PURSUE via remotes", { if (!has("remotes")) install.packages("remotes"); remotes::install_local(root, upgrade = "never", force = TRUE) })
} else cat("   PURSUE installed\n")
if (!has("PURSUE")) say("  !! PURSUE did not install. It is the method under test -- nothing can be measured without it.")

cat("\n>> 2. bootstrap\n")
if (!has("BiocManager")) try_("BiocManager", install.packages("BiocManager"))
if (!has("remotes")) try_("remotes", install.packages("remotes"))

# Deriv is a hard dependency of LOCOM2 that CRAN does not pull in automatically here; it must
# be listed before LOCOM2 (smoke test 2026-09-11: "dependency 'Deriv' is not available").
cran <- c("optparse", "jsonlite", "data.table", "sandwich", "lme4", "lmerTest", "pscl", "ranger", "arrow",
          "Deriv", "GUniFrac", "MIDASim", "MicrobiomeStat", "corncob", "LDM", "LOCOM2", "radEmu",
          "vegan", "ape", "matrixStats", "logging", "multcomp")
bioc <- c("limma", "phyloseq", "ANCOMBC", "ALDEx2", "SPsimSeq", "ADAPT", "maaslin3", "SummarizedExperiment", "TreeSummarizedExperiment")
# ADAPT and maaslin3 entered Bioconductor after 3.18, so on an R 4.3 cluster they can only come
# from GitHub. If both keep failing, the real fix is a newer R module (see hpc/README.md).
gh   <- c(SparseDOSSA2 = "biobakery/SparseDOSSA2", LOCOM = "yijuanhu/LOCOM", fastANCOM = "ZRChao/fastANCOM",
          fastEmu = "statdivlab/fastEmu", maaslin3 = "biobakery/maaslin3", ADAPT = "mkbwang/ADAPT")
src <- c(setNames(rep("CRAN", length(cran)), cran), setNames(rep("Bioconductor", length(bioc)), bioc),
         setNames(paste0("GitHub:", gh), names(gh)), PURSUE = "this checkout")

# CRAN and Bioconductor in ONE repository list. LDM and LOCOM2 are CRAN packages that depend
# on BiocParallel; with only CRAN in `repos` they fail with "dependency 'BiocParallel' is not
# available" even though BiocManager is installed (smoke test 2026-09-11).
if (has("BiocManager")) {
  options(repos = BiocManager::repositories())
  cat("\n>> repositories: ", paste(names(getOption("repos")), collapse = ", "), "\n", sep = "")
}
cat("\n>> 3. CRAN\n");         for (p in cran) if (!has(p)) try_(p, install.packages(p))
cat("\n>> 4. Bioconductor\n")
if (has("BiocManager")) { cat("   Bioconductor", as.character(BiocManager::version()), "\n")
  for (p in bioc) if (!has(p)) try_(p, BiocManager::install(p, ask = FALSE, update = FALSE))
} else message("  !! BiocManager unavailable: skipping ", length(bioc), " Bioconductor packages")
cat("\n>> 5. GitHub\n")
if (has("remotes")) { for (p in names(gh)) if (!has(p)) try_(p, remotes::install_github(gh[[p]], upgrade = "never"))
} else message("  !! remotes unavailable: skipping ", length(gh), " GitHub packages")

cat("\n>> 5b. known incompatibilities\n")
# --- known incompatibility: ANCOMBC <-> CVXR --------------------------------------------
# ANCOMBC (2.12.0 on Bioc 3.22) has importFrom(CVXR, solve) in its NAMESPACE, but current CVXR
# no longer exports `solve` -- it was renamed `psolve` to stop colliding with base::solve. The
# install dies at "object 'solve' is not exported by 'namespace:CVXR'" during lazy loading.
# CVXR is used by nothing else here, so stepping it back is safe. Versions are tried newest
# first and the export is checked directly rather than trusted from a version number.
cvxr_ok <- function() has("CVXR") && "solve" %in% getNamespaceExports("CVXR")
if (!has("ANCOMBC")) {
  if (!cvxr_ok()) {
    cur <- if (has("CVXR")) as.character(utils::packageVersion("CVXR")) else "none"
    say("  ~  ANCOMBC needs CVXR::solve, absent from CVXR ", cur, " -- trying older CVXR")
    if (!has("remotes")) try_("remotes", install.packages("remotes"))
    for (v in c("1.0.15", "1.0.14", "1.0.12", "1.0.11")) {
      if (cvxr_ok()) break
      cat("   trying CVXR ", v, "\n", sep = "")
      try_(paste0("CVXR ", v), remotes::install_version("CVXR", version = v, upgrade = "never", quiet = TRUE))
      try(unloadNamespace("CVXR"), silent = TRUE)   # stale namespace would mask the new install
    }
  }
  if (cvxr_ok()) {
    cat(">> CVXR ", as.character(utils::packageVersion("CVXR")), " exports solve; retrying ANCOMBC\n", sep = "")
    try_("ANCOMBC", BiocManager::install("ANCOMBC", ask = FALSE, update = FALSE))
  } else say("  !! no CVXR version exporting `solve` could be installed; ANCOMBC stays unavailable")
  if (!has("ANCOMBC"))
    try_("ANCOMBC from GitHub", remotes::install_github("FrederickHuangLin/ANCOMBC", upgrade = "never"))
}

# PURSUE is attempted FIRST so that it lands even with no network, but R CMD INSTALL refuses a
# package whose Imports are absent -- in a brand-new environment limma and sandwich do not exist
# yet, which is exactly why it failed on 2026-09-11. Retry now that everything else is in.
if (!has("PURSUE")) {
  cat("\n>> 6. PURSUE (retry, now that its dependencies are installed)\n")
  miss_dep <- Filter(function(d) !has(d), c("limma", "sandwich"))
  if (length(miss_dep)) say("  !! PURSUE needs ", paste(miss_dep, collapse = ", "), ", which did not install")
  pi2 <- suppressWarnings(system2(file.path(R.home("bin"), "R"), c("CMD", "INSTALL", shQuote(root)),
                                  stdout = TRUE, stderr = TRUE))
  if (!is.null(attr(pi2, "status")) && attr(pi2, "status") != 0) {
    say("  !! PURSUE retry: R CMD INSTALL exited ", attr(pi2, "status"), ". Full output:")
    cat(paste0("     ", pi2, collapse = "\n"), "\n", file = logf, append = TRUE)
    cat(paste(utils::tail(pi2, 25), collapse = "\n"), "\n")
  } else cat("   PURSUE installed on retry\n")
}

all <- unique(c("PURSUE", cran, bioc, names(gh)))
st <- data.frame(package = all, source = unname(src[all]), installed = vapply(all, has, logical(1)),
                 version = vapply(all, function(p) if (has(p)) as.character(utils::packageVersion(p)) else NA_character_, character(1)),
                 stringsAsFactors = FALSE)
st <- st[order(st$installed, st$package), ]
cat("\n"); print(st, row.names = FALSE)
dir.create(file.path(root, "hpc"), showWarnings = FALSE)
write.csv(st, file.path(root, "hpc", "installed_packages.csv"), row.names = FALSE)
# Re-run each remaining failure in a subprocess with ALL output captured. "installation of
# package 'X' failed" on its own is not actionable and has cost a full round trip more than
# once; the real reason is in R CMD INSTALL's output, which install.packages() only warns about.
miss0 <- st$package[!st$installed & st$package != "PURSUE"]
if (length(miss0)) {
  cat("\n>> 7. diagnosing ", length(miss0), " failure(s) -- full build output -> hpc/install_log.txt\n", sep = "")
  pre <- 'if (requireNamespace("BiocManager", quietly=TRUE)) options(repos = BiocManager::repositories());'
  for (pk in miss0) {
    kind <- unname(src[pk])
    cmd <- if (identical(kind, "CRAN")) sprintf('install.packages("%s")', pk)
           else if (identical(kind, "Bioconductor")) sprintf('BiocManager::install("%s", ask=FALSE, update=FALSE)', pk)
           else sprintf('remotes::install_github("%s", upgrade="never")', sub("^GitHub:", "", kind))
    out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
                                    c("-e", shQuote(paste(pre, cmd))), stdout = TRUE, stderr = TRUE))
    cat("\n=== ", pk, " (", kind, ")\n", sep = "", file = logf, append = TRUE)
    cat(paste0("   ", out, collapse = "\n"), "\n", file = logf, append = TRUE)
    err <- grep("^(ERROR|Error|error:|\\*\\* .*ERROR)", out, value = TRUE)
    cat("   ", pk, ": ", if (length(err)) trimws(paste(utils::head(err, 2), collapse = " | ")) else "see hpc/install_log.txt", "\n", sep = "")
  }
  # a diagnosis run can also succeed; refresh the table so the CSV is not stale
  st$installed <- vapply(st$package, has, logical(1))
  st$version <- vapply(st$package, function(q) if (has(q)) as.character(utils::packageVersion(q)) else NA_character_, character(1))
  write.csv(st, file.path(root, "hpc", "installed_packages.csv"), row.names = FALSE)
}

miss <- st$package[!st$installed]
cat(sprintf("\n%d of %d installed. Missing: %s\n", sum(st$installed), nrow(st), if (length(miss)) paste(miss, collapse = ", ") else "none"))
if (length(miss)) cat("A missing package is not fatal -- its method/simulator is recorded as not_installed.\n",
                      "If many failed at once the node most likely has no outbound network; install from a login node or ask for a local mirror.\n", sep = "")
cat("wrote hpc/installed_packages.csv and hpc/install_log.txt\n")
