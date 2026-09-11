#!/usr/bin/env Rscript
# Install everything the benchmark can use. Idempotent; prints a final availability table.
# Run inside the environment created by setup_env.sh (or any R >= 4.3 with internet).
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = max(1L, parallel::detectCores() - 1L), timeout = 600)
lib <- Sys.getenv("R_LIBS_USER", unset = ""); if (nzchar(lib)) { dir.create(lib, recursive = TRUE, showWarnings = FALSE); .libPaths(c(lib, .libPaths())) }
has <- function(p) requireNamespace(p, quietly = TRUE)
if (!has("BiocManager")) install.packages("BiocManager")
if (!has("remotes")) install.packages("remotes")
cran <- c("optparse", "jsonlite", "yaml", "data.table", "sandwich", "lme4", "pscl", "ranger", "arrow",
          "GUniFrac", "MIDASim", "MicrobiomeStat", "corncob", "LDM", "LOCOM2", "radEmu", "vegan", "ape", "matrixStats")
bioc <- c("limma", "phyloseq", "ANCOMBC", "ALDEx2", "SPsimSeq", "ADAPT", "maaslin3", "SummarizedExperiment", "TreeSummarizedExperiment")
gh   <- c(SparseDOSSA2 = "biobakery/SparseDOSSA2", LOCOM = "yijuanhu/LOCOM", fastANCOM = "ZRChao/fastANCOM", fastEmu = "statdivlab/fastEmu",
          maaslin3 = "biobakery/maaslin3")     # maaslin3 from GitHub if the Bioconductor release is too old
for (p in cran) if (!has(p)) tryCatch(install.packages(p), error = function(e) message("CRAN failed: ", p))
for (p in bioc) if (!has(p)) tryCatch(BiocManager::install(p, ask = FALSE, update = FALSE), error = function(e) message("Bioc failed: ", p))
for (p in names(gh)) if (!has(p)) tryCatch(remotes::install_github(gh[[p]], upgrade = "never"), error = function(e) message("GitHub failed: ", p))
# PURSUE itself, from this repository
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), ".."), mustWork = FALSE)
tryCatch(remotes::install_local(root, upgrade = "never", force = TRUE), error = function(e) message("PURSUE install failed: ", conditionMessage(e)))
all <- c(cran, bioc, names(gh), "PURSUE")
st <- data.frame(package = all, installed = vapply(all, has, logical(1)),
                 version = vapply(all, function(p) if (has(p)) as.character(packageVersion(p)) else NA_character_, character(1)))
print(st, row.names = FALSE)
write.csv(st, file.path(root, "hpc", "installed_packages.csv"), row.names = FALSE)
