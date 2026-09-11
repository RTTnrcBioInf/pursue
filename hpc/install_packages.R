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
try_ <- function(what, expr) tryCatch({ expr; TRUE }, error = function(e) { message("  !! ", what, ": ", conditionMessage(e)); FALSE },
                                                              warning = function(w) { message("  ~  ", what, ": ", conditionMessage(w)); has(what) })
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), ".."), mustWork = FALSE)

cat(">> R", R.version.string, "\n>> library:", .libPaths()[1], "\n\n>> 1. PURSUE (from this checkout)\n")
if (!try_("PURSUE", install.packages(root, repos = NULL, type = "source")))
  try_("PURSUE via remotes", { if (!has("remotes")) install.packages("remotes"); remotes::install_local(root, upgrade = "never", force = TRUE) })

cat("\n>> 2. bootstrap\n")
if (!has("BiocManager")) try_("BiocManager", install.packages("BiocManager"))
if (!has("remotes")) try_("remotes", install.packages("remotes"))

cran <- c("optparse", "jsonlite", "data.table", "sandwich", "lme4", "pscl", "ranger", "arrow",
          "GUniFrac", "MIDASim", "MicrobiomeStat", "corncob", "LDM", "LOCOM2", "radEmu", "vegan", "ape", "matrixStats")
bioc <- c("limma", "phyloseq", "ANCOMBC", "ALDEx2", "SPsimSeq", "ADAPT", "maaslin3", "SummarizedExperiment", "TreeSummarizedExperiment")
gh   <- c(SparseDOSSA2 = "biobakery/SparseDOSSA2", LOCOM = "yijuanhu/LOCOM", fastANCOM = "ZRChao/fastANCOM",
          fastEmu = "statdivlab/fastEmu", maaslin3 = "biobakery/maaslin3")
src <- c(setNames(rep("CRAN", length(cran)), cran), setNames(rep("Bioconductor", length(bioc)), bioc),
         setNames(paste0("GitHub:", gh), names(gh)), PURSUE = "this checkout")

cat("\n>> 3. CRAN\n");         for (p in cran) if (!has(p)) try_(p, install.packages(p))
cat("\n>> 4. Bioconductor\n")
if (has("BiocManager")) { cat("   Bioconductor", as.character(BiocManager::version()), "\n")
  for (p in bioc) if (!has(p)) try_(p, BiocManager::install(p, ask = FALSE, update = FALSE))
} else message("  !! BiocManager unavailable: skipping ", length(bioc), " Bioconductor packages")
cat("\n>> 5. GitHub\n")
if (has("remotes")) { for (p in names(gh)) if (!has(p)) try_(p, remotes::install_github(gh[[p]], upgrade = "never"))
} else message("  !! remotes unavailable: skipping ", length(gh), " GitHub packages")

all <- unique(c("PURSUE", cran, bioc, names(gh)))
st <- data.frame(package = all, source = unname(src[all]), installed = vapply(all, has, logical(1)),
                 version = vapply(all, function(p) if (has(p)) as.character(utils::packageVersion(p)) else NA_character_, character(1)),
                 stringsAsFactors = FALSE)
st <- st[order(st$installed, st$package), ]
cat("\n"); print(st, row.names = FALSE)
dir.create(file.path(root, "hpc"), showWarnings = FALSE)
write.csv(st, file.path(root, "hpc", "installed_packages.csv"), row.names = FALSE)
miss <- st$package[!st$installed]
cat(sprintf("\n%d of %d installed. Missing: %s\n", sum(st$installed), nrow(st), if (length(miss)) paste(miss, collapse = ", ") else "none"))
if (length(miss)) cat("A missing package is not fatal -- its method/simulator is recorded as not_installed.\n",
                      "If many failed at once the node most likely has no outbound network; install from a login node or ask for a local mirror.\n", sep = "")
cat("wrote hpc/installed_packages.csv\n")
