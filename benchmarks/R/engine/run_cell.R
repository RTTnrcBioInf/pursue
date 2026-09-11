#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# run_cell.R -- run every requested method on ONE benchmark cell and write the
# results contract. This is what a SLURM array task calls.
#
#   Rscript R/engine/run_cell.R --axis A --simulator house --template hmp_stool \
#           --regime R00 --replicate 1 --methods pursue,limma_logtss,linda --out results/
#
# Axis B: --simulator implant with --regime naming an implantation spec from
#         benchmarks/implant_specs.tsv.
# Axis C: --simulator nullperm --regime <scheme>   (full | stratified | depth_quintile | cluster)
# Seeds: seed = master_seed * 1e6 + hash(axis, simulator, template, regime, replicate) mod 1e6,
#        master_seed from --master-seed (default 1).
# -----------------------------------------------------------------------------
suppressMessages({ library(optparse); library(jsonlite) })
op <- OptionParser(option_list = list(
  make_option("--axis", type = "character"), make_option("--simulator", type = "character"),
  make_option("--template", type = "character"), make_option("--regime", type = "character"),
  make_option("--replicate", type = "integer", default = 1L), make_option("--methods", type = "character", default = "all"),
  make_option("--out", type = "character", default = "results"), make_option("--master-seed", type = "integer", default = 1L),
  make_option("--cache", type = "character", default = "cache"), make_option("--min-prevalence", type = "double", default = 0.10),
  make_option("--timeout", type = "integer", default = 3600L), make_option("--bench-root", type = "character", default = NULL),
  make_option("--n-cores", type = "integer", default = 1L)))
opt <- parse_args(op)
if (!is.null(opt$`bench-root`)) Sys.setenv(PURSUE_BENCH_ROOT = opt$`bench-root`)
root <- Sys.getenv("PURSUE_BENCH_ROOT", unset = ".")
for (f in c("R/engine/templates.R", "R/engine/regimes.R", "R/engine/metrics.R", "R/engine/io.R", "R/engine/realism.R",
            "R/simulators/sim_house.R", "R/simulators/implant.R", "R/simulators/nullperm.R", "R/simulators/dispatch.R",
            "R/methods/elementary.R", "R/methods/external.R", "R/methods/registry.R")) source(file.path(root, f))

cell_seed <- function(master, ...) { key <- paste(..., sep = "|"); h <- sum(utf8ToInt(key) * (seq_along(utf8ToInt(key)) %% 97 + 1)) %% 1e6
  as.integer(master * 1e6 + h) %% .Machine$integer.max }

cell <- list(axis = opt$axis, simulator = opt$simulator, template = opt$template, regime_id = opt$regime, replicate = opt$replicate)
cell$seed <- cell_seed(opt$`master-seed`, opt$axis, opt$simulator, opt$template, opt$regime, opt$replicate)
methods <- if (opt$methods == "all") method_registry()$id else strsplit(opt$methods, ",")[[1]]
dir.create(opt$cache, recursive = TRUE, showWarnings = FALSE)

tpl <- load_template(opt$template)
sim <- switch(opt$axis,
  A = { reg <- read.delim(file.path(root, "regimes.tsv"), stringsAsFactors = FALSE, comment.char = "#"); reg <- reg[reg$regime_id == opt$regime, ]
        if (nrow(reg) != 1L) stop("regime not found: ", opt$regime)
        if (!simulator_available(opt$simulator)) list(unsupported = "simulator_not_installed") else simulate_cell_data(opt$simulator, tpl, reg, cell$seed, opt$cache) },
  B = { sp <- read.delim(file.path(root, "implant_specs.tsv"), stringsAsFactors = FALSE, comment.char = "#"); sp <- sp[sp$spec_id == opt$regime, ]
        if (nrow(sp) != 1L) stop("implant spec not found: ", opt$regime); implant(tpl, as.list(sp), cell$seed) },
  C = null_permute(tpl, scheme = opt$regime, seed = cell$seed),
  stop("axis must be A, B or C for run_cell.R (D and E have their own drivers)"))

if (!is.null(sim$unsupported)) {
  writeLines(toJSON(c(cell, list(status = "regime_unsupported", unsupported = sim$unsupported)), auto_unbox = TRUE),
             file.path(opt$out, paste0(cell_id(cell), ".skipped.json")))
  cat("skipped:", paste(sim$unsupported, collapse = ","), "\n"); quit(status = 0)
}

# common prevalence filter applied before every method
keep <- rowMeans(sim$counts > 0) >= opt$`min-prevalence` & rowSums(sim$counts > 0) >= 3
counts <- sim$counts[keep, , drop = FALSE]; truth <- sim$truth[match(rownames(counts), sim$truth$feature), ]
meta <- sim$meta
cat(sprintf("cell %s: %d samples x %d features (%d after filter), %d methods\n", cell_id(cell), ncol(sim$counts), nrow(sim$counts), nrow(counts), length(methods)))

contract <- list(); metrics <- list()
for (mid in methods) {
  run <- run_method(mid, counts, meta, sim$formula, sim$tested_term, args = if (mid == "pursue") list(n_cores = opt$`n-cores`) else list(), timeout_s = opt$timeout)
  cc <- assemble_contract(cell, mid, run, truth, run$version)
  mm <- cell_metrics(run$result, truth); mm <- cbind(method = mid, method_version = run$version, runtime_s = run$runtime_s, mm)
  contract[[mid]] <- cc; metrics[[mid]] <- mm
  cat(sprintf("  %-18s %6.1fs  %s\n", mid, run$runtime_s, paste(unique(run$result$status), collapse = "/")))
}
contract <- do.call(rbind, contract); metrics <- cbind(axis = cell$axis, simulator = cell$simulator, template = cell$template,
                                                       regime_id = cell$regime_id, replicate = cell$replicate, do.call(rbind, metrics))
extra <- list(n_samples = ncol(sim$counts), n_features_total = nrow(sim$counts), n_features_tested = nrow(counts),
              tested_term = sim$tested_term, formula = paste(deparse(sim$formula), collapse = ""), min_prevalence = opt$`min-prevalence`)
id <- write_cell(cell, contract, metrics, extra, opt$out)
cat("wrote", id, "\n")
