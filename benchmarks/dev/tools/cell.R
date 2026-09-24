# tools/cell.R -- rebuild one dev-suite cell exactly as devsuite.R does, for interactive diagnosis.
#   source("benchmarks/dev/tools/cell.R"); cl <- dev_cell("house:R08", "twinsuk_stool", 3)
#   -> list(counts (filtered), meta, formula, tested_term, truth (aligned to counts rows), sim)
root <- normalizePath(Sys.getenv("PURSUE_BENCH_ROOT", unset = "benchmarks")); Sys.setenv(PURSUE_BENCH_ROOT = root)
source(file.path(root, "dev", "tools", "stress.R"))
for (f in c("R/engine/templates.R", "R/engine/regimes.R", "R/engine/metrics.R", "R/simulators/sim_house.R",
            "R/simulators/implant.R", "R/simulators/dispatch.R", "R/methods/elementary.R")) source(file.path(root, f))
.cands <- new.env(); register_candidate <- function(name, fn, notes = "") assign(name, list(fn = fn), envir = .cands)
for (f in sort(list.files(file.path(root, "dev", "candidates"), pattern = "\\.R$", full.names = TRUE))) source(f)
.regimes <- read.delim(file.path(root, "regimes.tsv"), stringsAsFactors = FALSE, comment.char = "#")
.spec <- function(signal = "mixed", da = 0.10, balance = "balanced") list(n_per_group = 50, da_frac = da, effect = "medium",
  signal_type = signal, balance = balance, conf_phi = 0, exposure = "binary")
.imp <- list(B_ref = .spec(), B_null = .spec(da = 0), B_prev = .spec("prevalence"), B_abund = .spec("abundance"), B_allup = .spec(balance = "100_0"))
dev_cell <- function(sid, tid, rep, seed = 9001L) {
  cs <- function(...) { k <- paste(..., sep = "|"); as.integer((seed * 1e5 + sum(utf8ToInt(k) * seq_along(utf8ToInt(k)))) %% .Machine$integer.max) }
  tpl <- readRDS(file.path(root, "devdata", paste0(tid, ".rds"))); kind <- sub(":.*", "", sid); rg <- sub(".*:", "", sid)
  sim <- if (kind == "implant") implant(tpl, .imp[[rg]], cs(sid, tid, rep)) else if (kind == "bloom")
    apply_bloom(simulate_cell_data("house", tpl, .regimes[.regimes$regime_id == "R00", ], cs(sid, tid, rep)), as.numeric(sub("x", "", rg)), cs(sid, tid, rep)) else
    simulate_cell_data(kind, tpl, .regimes[.regimes$regime_id == rg, ], cs(sid, tid, rep))
  keep <- rowMeans(sim$counts > 0) >= 0.10 & rowSums(sim$counts > 0) >= 3; ct <- sim$counts[keep, , drop = FALSE]
  meta <- sim$meta[colnames(ct), , drop = FALSE]; if (is.null(meta$depth)) meta$depth <- colSums(sim$counts)
  list(counts = ct, meta = meta, formula = sim$formula, tested_term = sim$tested_term,
       truth = sim$truth[match(rownames(ct), sim$truth$feature), ], sim = sim)
}
