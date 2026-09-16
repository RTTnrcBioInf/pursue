#!/usr/bin/env Rscript
# pursue_presence_ablation.R ---------------------------------------------------
# Axes A and B agree on two things about PURSUE 0.2:
#   (a) equal-weight Cauchy combination costs 11-36% of the power the presence arm alone has;
#   (b) far more important, the presence arm itself is beaten ~3x by plain logistic regression
#       on detection (axis B B_ref: 2.07 true positives per cell vs 6.58).
# (b) dwarfs (a), so the question for 0.3 is not how to combine the arms but WHY a stabilized
# zero-inflated negative binomial on structural presence loses to `glm(detected ~ group)`.
#
# Four candidates, each switched off in turn:
#   ridge      tau_ridge shrinks the presence log-odds toward zero (default SD 5)
#   theta      d0_prior shrinks the NB dispersion across features (default 10)
#   depth      the depth offset, which is what makes presence "structural" rather than observed
#   the model  ZINB-LRT vs plain logistic on the same filtered table
#
# The hypothesis worth killing first: PURSUE tests presence NET of abundance, and in these
# designs the abundance change IS what moves detection -- so the NB component absorbs the group
# effect and the zero-inflation component is left with little. If that is right, no amount of
# tuning fixes it and the estimand itself has to change.
#
#   Rscript hpc/pursue_presence_ablation.R [--reps 10] [--template hmp_stool]
# Runtime: a few minutes (PURSUE is ~15 s per cell).
# ------------------------------------------------------------------------------
suppressMessages(library(optparse))
op <- OptionParser(option_list = list(
  make_option("--reps", type = "integer", default = 10L),
  make_option("--template", default = "hmp_stool"),
  make_option("--out", default = "hpc/pursue_ablation.csv"),
  make_option("--bench-root", default = NULL)))
opt <- parse_args(op)
if (!is.null(opt$`bench-root`)) Sys.setenv(PURSUE_BENCH_ROOT = opt$`bench-root`)
root <- Sys.getenv("PURSUE_BENCH_ROOT", unset = "benchmarks")
for (f in c("R/engine/templates.R", "R/engine/regimes.R", "R/engine/metrics.R",
            "R/simulators/sim_house.R", "R/simulators/implant.R", "R/simulators/dispatch.R",
            "R/methods/elementary.R", "R/methods/external.R", "R/methods/registry.R")) source(file.path(root, f))

tpl <- load_template(opt$template)
specs <- list(
  B_prev  = list(n_per_group = 50, da_frac = 0.10, effect = "medium", signal_type = "prevalence",
                 balance = "balanced", conf_phi = 0, exposure = "binary"),
  B_ref   = list(n_per_group = 50, da_frac = 0.10, effect = "medium", signal_type = "mixed",
                 balance = "balanced", conf_phi = 0, exposure = "binary"))

# Each variant is the presence arm computed a different way, scored identically.
variants <- list(
  default      = list(),
  no_ridge     = list(tau_ridge = 1e6),      # ridge penalty effectively off
  no_theta     = list(d0_prior = 0),         # no cross-feature dispersion shrinkage
  no_depth     = list(depth_adjust = FALSE), # abundance arm's depth covariate off
  no_ridge_theta = list(tau_ridge = 1e6, d0_prior = 0))

rows <- list()
for (sp_id in names(specs)) for (rep in seq_len(opt$reps)) {
  sim <- implant(tpl, specs[[sp_id]], seed = 1000L + rep)
  if (!is.null(sim$unsupported)) next
  keep <- rowMeans(sim$counts > 0) >= 0.10 & rowSums(sim$counts > 0) >= 3
  ct <- sim$counts[keep, , drop = FALSE]
  tr <- sim$truth[match(rownames(ct), sim$truth$feature), ]
  y  <- truth_for(tr, "presence")                       # prevalence-type truth
  ya <- truth_for(tr, "abundance")

  score <- function(label, p, truth_vec) {
    ok <- is.finite(p); if (!any(ok) || sum(truth_vec[ok] == 1) == 0) return(NULL)
    q <- p.adjust(p, "BH"); rej <- !is.na(q) & q <= 0.05
    data.frame(spec = sp_id, rep = rep, variant = label,
               n_pos = sum(truth_vec[ok] == 1), n_rej = sum(rej, na.rm = TRUE),
               tp = sum(rej & truth_vec == 1, na.rm = TRUE),
               fp = sum(rej & truth_vec == 0, na.rm = TRUE), stringsAsFactors = FALSE)
  }

  for (v in names(variants)) {
    r <- tryCatch(method_pursue(ct, sim$meta, sim$formula, sim$tested_term, variants[[v]]),
                  error = function(e) NULL)
    if (is.null(r)) { message("pursue/", v, " failed on ", sp_id, " r", rep); next }
    pr <- r[r$arm == "presence", ]; pr <- pr[match(rownames(ct), pr$feature), ]
    rows[[length(rows) + 1L]] <- score(paste0("pursue_presence:", v), pr$p, y)
    ab <- r[r$arm == "abundance", ]; ab <- ab[match(rownames(ct), ab$feature), ]
    rows[[length(rows) + 1L]] <- score(paste0("pursue_abundance:", v), ab$p, ya)
  }
  # the comparator that is beating it
  lg <- method_logistic_presence(ct, sim$meta, sim$formula, sim$tested_term)
  lg <- lg[match(rownames(ct), lg$feature), ]
  rows[[length(rows) + 1L]] <- score("logistic_presence", lg$p, y)
  cat(sprintf("  %s rep %d done\n", sp_id, rep))
}

D <- do.call(rbind, rows)
if (is.null(D)) { cat("no rows produced\n"); quit(status = 1) }
write.csv(D, opt$out, row.names = FALSE)

agg <- aggregate(cbind(n_pos, n_rej, tp, fp) ~ spec + variant, D, mean)
agg$FDR <- round(agg$fp / pmax(agg$n_rej, 1e-9), 3)
agg <- agg[order(agg$spec, -agg$tp), ]
cat("\n=== true positives per cell, presence-type truth (higher is better) ===\n")
for (s in unique(agg$spec)) {
  cat("\n", s, "\n", sep = "")
  a <- agg[agg$spec == s, ]
  print(data.frame(variant = a$variant, n_pos = round(a$n_pos, 1), TP = round(a$tp, 2),
                   FP = round(a$fp, 2), FDR = a$FDR), row.names = FALSE)
}
cat("\nRead it like this:\n")
cat("  a variant near logistic_presence  -> that component was the whole power loss; fix it.\n")
cat("  every variant still far below it  -> the ZINB structural-presence estimand is the\n")
cat("                                       problem, not its tuning, and 0.3 needs a new arm.\n")
cat("\nwrote ", opt$out, "\n", sep = "")
