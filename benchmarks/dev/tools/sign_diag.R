#!/usr/bin/env Rscript
# -----------------------------------------------------------------------------
# sign_diag.R -- does c14 report the right direction for the taxa it calls? (2026-11-03; no method change)
# conf2 read c14's effect estimates as its weak point: on Vandeputte (measured load) sign agreement 0.85 against LinDA's
# 0.96, and on three replication sets 4-7% of findings called in both halves with opposite signs. Hypothesis: the
# reported estimate always comes from the log kernel, while the call can come from another kernel (detection, sqrt,
# abundance given presence); where those disagree -- above all taxa whose presence and abundance move in opposite
# directions -- the reported sign follows a different signal from the one the test found.
# Development data only (tuning templates, dev seed). For every taxon c14 calls (its own q <= 0.05) that truly changes:
# the true direction, the sign of the reported estimate, which kernel made the call and the sign of that kernel's
# statistic. LinDA alongside, on the same cells.
# Settings: house R00 / R11 / R13 / R16 and implants as the dev suite, plus two implants this tool adds: B_opp (every
# changed taxon moves presence one way and abundance among the present the other) and B_mix4 (abundance / prevalence /
# both / opposite in equal parts). For an opposite taxon the "true direction" depends on the estimand, so four are kept:
# presence, abundance given presence, mean abundance, and mean log abundance with a small pseudocount (the convention of
# axis D's measured-load truth), all from the absolute values before resampling.
#   Rscript benchmarks/dev/tools/sign_diag.R --reps 5 --cores 2 --label signdiag1
# -----------------------------------------------------------------------------
suppressMessages(library(optparse))
op <- OptionParser(option_list = list(
  make_option("--settings", default = "house:R00,house:R11,house:R13,house:R16,implant:B_ref,implant:B_prev,implant:B_abund,implant:B_allup,implant:B_opp,implant:B_mix4"),
  make_option("--templates", default = "hmp_tongue,twinsuk_stool"),
  make_option("--reps", type = "integer", default = 5L), make_option("--cores", type = "integer", default = 2L),
  make_option("--label", default = "signdiag"), make_option("--seed", type = "integer", default = 9001L)))
opt <- parse_args(op)
root <- normalizePath(Sys.getenv("PURSUE_BENCH_ROOT", unset = "benchmarks"), mustWork = TRUE)
Sys.setenv(PURSUE_BENCH_ROOT = root); if (!nzchar(Sys.getenv("PURSUE_DATA_ROOT"))) Sys.setenv(PURSUE_DATA_ROOT = file.path(root, "data"))
for (f in c("R/engine/templates.R", "R/engine/regimes.R", "R/engine/metrics.R", "R/simulators/sim_house.R", "R/simulators/implant.R",
            "R/simulators/dispatch.R", "R/methods/elementary.R")) source(file.path(root, f))
.cands <- new.env(); register_candidate <- function(name, fn, notes = "") assign(name, list(fn = fn, notes = notes), envir = .cands)
for (f in sort(list.files(file.path(root, "dev", "candidates"), pattern = "\\.R$", full.names = TRUE))) source(f)
regimes <- read.delim(file.path(root, "regimes.tsv"), stringsAsFactors = FALSE, comment.char = "#")
cache_dir <- file.path(dirname(root), "cache")

# implant() with an "opp" kind (presence and abundance in opposite directions) and direction truths from the absolute values
implant_dir <- function(template, spec, seed) {
  set.seed(seed); ct <- template$counts; m <- nrow(ct); n <- 2L * spec$n_per_group
  if (n > ncol(ct)) return(list(unsupported = "n"))
  X <- ct[, sample.int(ncol(ct), n), drop = FALSE]; x <- as.integer(rank(stats::runif(n)) > n / 2)
  es <- effect_sizes(spec$effect); n_da <- round(spec$da_frac * m)
  typ <- rep("none", m); a_dir <- numeric(m); p_dir <- numeric(m); lfc2 <- numeric(m); Xf <- X * 1.0
  da <- sample.int(m, n_da)
  kind <- switch(spec$signal_type, abundance = rep("abundance", n_da), prevalence = rep("prevalence", n_da), opp = rep("opp", n_da),
                 mixed = sample(c("abundance", "prevalence", "both"), n_da, TRUE), mixed4 = sample(c("abundance", "prevalence", "both", "opp"), n_da, TRUE))
  p_up <- switch(spec$balance, balanced = 0.5, "100_0" = 1.0); sgn <- ifelse(stats::runif(n_da) < p_up, 1, -1)
  mag_p <- switch(spec$effect, small = 0.2, medium = 0.35, large = 0.5)
  typ[da] <- kind
  for (k in seq_len(n_da)) {
    j <- da[k]
    if (kind[k] %in% c("abundance", "both", "opp")) { lfc2[j] <- sgn[k] * es$lfc2; a_dir[j] <- sgn[k]; Xf[j, ] <- Xf[j, ] * ifelse(x == 1L, 2^lfc2[j], 1) }
    if (kind[k] %in% c("prevalence", "both", "opp")) {
      sp <- if (kind[k] == "opp") -sgn[k] else sgn[k]; p_dir[j] <- sp; up <- if (sp > 0) x == 1L else x == 0L
      nz_down <- which(!up & Xf[j, ] > 0); z_up <- which(up & Xf[j, ] == 0); k_move <- min(length(nz_down), length(z_up), round(mag_p * length(nz_down)))
      if (k_move > 0) { med_nz <- stats::median(Xf[j, Xf[j, ] > 0]); Xf[j, sample(nz_down, k_move)] <- 0; Xf[j, sample(z_up, k_move)] <- med_nz }
    }
  }
  A <- sweep(Xf, 2L, colSums(X), "/")                                   # absolute scale: the template's own depth as unit
  dir_of <- function(v) { c0 <- if (any(v > 0)) min(v[v > 0]) / 2 else 1
    c(mean = sign(log(mean(v[x == 1]) + 1e-300) - log(mean(v[x == 0]) + 1e-300)), logmean = sign(mean(log2(v[x == 1] + c0)) - mean(log2(v[x == 0] + c0)))) }
  dd <- t(vapply(seq_len(m), function(j) if (typ[j] == "none") c(mean = 0, logmean = 0) else dir_of(A[j, ]), numeric(2)))
  N <- colSums(X); Xn <- sapply(seq_len(n), function(i) { p <- Xf[, i]; if (sum(p) <= 0) rep(0L, m) else stats::rmultinom(1L, N[i], p / sum(p))[, 1] })
  dimnames(Xn) <- list(rownames(ct), paste0("S", seq_len(n))); storage.mode(Xn) <- "integer"
  truth <- data.frame(feature = rownames(ct), truth_type = typ, truth_abs = as.integer(typ != "none"), truth_lfc2 = lfc2,
                      dir_presence = p_dir, dir_abund = a_dir, dir_mean = dd[, "mean"], dir_logmean = dd[, "logmean"], stringsAsFactors = FALSE)
  meta <- data.frame(row.names = colnames(Xn), group = factor(ifelse(x == 1L, "case", "control"), levels = c("control", "case")), depth = N)
  list(counts = Xn, meta = meta, truth = truth, tested_term = "group", formula = ~ group)
}
spec <- function(signal = "mixed", balance = "balanced") list(n_per_group = 50, da_frac = 0.10, effect = "medium", signal_type = signal, balance = balance)
IMP <- list(B_ref = spec(), B_prev = spec("prevalence"), B_abund = spec("abundance"), B_allup = spec(balance = "100_0"), B_opp = spec("opp"), B_mix4 = spec("mixed4"))
make_cell <- function(sid, tpl, seed) {
  sm <- sub(":.*", "", sid); rg <- sub(".*:", "", sid)
  if (sm == "implant") return(implant_dir(tpl, IMP[[rg]], seed))
  x <- simulate_cell_data(sm, tpl, regimes[regimes$regime_id == rg, ], seed, cache_dir); t <- x$truth
  t$dir_abund <- sign(t$truth_lfc2); t$dir_presence <- sign(if (is.null(t$truth_logor)) 0 else t$truth_logor)
  t$dir_mean <- ifelse(t$dir_abund != 0, t$dir_abund, t$dir_presence); t$dir_logmean <- t$dir_mean; x$truth <- t; x
}
cell_seed <- function(...) { k <- paste(..., sep = "|"); as.integer((opt$seed * 1e5 + sum(utf8ToInt(k) * seq_along(utf8ToInt(k)))) %% .Machine$integer.max) }
fit_c14 <- function(ct, meta, f, tt) .with_opts({
  v <- if (.is_paired(meta, f, tt)) .epp_fit(ct, meta, f, tt, lin = TRUE, kapfree = 20L, rho1 = TRUE) else
    .epg_fit(ct, meta, f, tt, paths = "cont", adet = TRUE, nodet_conf = TRUE, pk = "pos", minsw = 2L)
  v$q <- .storey_gq(v$p_full); v }, pursue.flstat = "cor", pursue.own = "pos")

sets <- strsplit(opt$settings, ",")[[1]]; tpls <- strsplit(opt$templates, ",")[[1]]
TPL <- setNames(lapply(tpls, function(id) readRDS(file.path(root, "devdata", paste0(id, ".rds")))), tpls)
jobs <- expand.grid(s = sets, tp = tpls, rep = seq_len(opt$reps), stringsAsFactors = FALSE)
odir <- file.path(root, "dev", "results", opt$label); cdir <- file.path(odir, "cells"); dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
cat(sprintf("sign diagnostics '%s': %d cells, %d cores\n", opt$label, nrow(jobs), opt$cores))
run_job <- function(k) {
  s <- jobs$s[k]; tid <- jobs$tp[k]; rep <- jobs$rep[k]; ck <- file.path(cdir, sprintf("%s__%s__r%02d.csv", gsub(":", "-", s), tid, rep))
  if (file.exists(ck)) return(utils::read.csv(ck, stringsAsFactors = FALSE))
  sim <- make_cell(s, TPL[[tid]], cell_seed(s, tid, rep)); if (!is.null(sim$unsupported)) return(NULL)
  keep <- rowMeans(sim$counts > 0) >= 0.10 & rowSums(sim$counts > 0) >= 3; ct <- sim$counts[keep, , drop = FALSE]
  meta <- sim$meta[colnames(ct), , drop = FALSE]; if (is.null(meta$depth)) meta$depth <- colSums(sim$counts)
  tr <- sim$truth[match(rownames(ct), sim$truth$feature), ]
  v <- fit_c14(ct, meta, sim$formula, sim$tested_term); i <- match(rownames(ct), v$feature)
  Tm <- if (!is.null(v$T)) v$T[rownames(ct), , drop = FALSE] else NULL
  ch <- v$choice[i]; ksign <- rep(NA_real_, nrow(ct))
  if (!is.null(Tm)) for (r in seq_len(nrow(ct))) { h <- ch[r]; if (is.na(h)) next
    if (h == "lead") { a <- Tm[r, intersect(c("det", "log"), colnames(Tm))]; ksign[r] <- sign(a[which.max(abs(a))]) } else if (h %in% colnames(Tm)) ksign[r] <- sign(Tm[r, h]) }
  ld <- tryCatch(.cands$ref_linda$fn(ct, meta, sim$formula, sim$tested_term), error = function(e) NULL)
  lq <- if (!is.null(ld)) { q <- if (!is.null(ld$q)) ld$q else stats::p.adjust(ld$p, "BH"); q[match(rownames(ct), ld$feature)] } else rep(NA, nrow(ct))
  le <- if (!is.null(ld) && !is.null(ld$estimate)) ld$estimate[match(rownames(ct), ld$feature)] else rep(NA, nrow(ct))
  out <- data.frame(setting = s, template = tid, rep = rep, feature = rownames(ct), truth_type = tr$truth_type, truth_lfc2 = tr$truth_lfc2,
                    dir_presence = tr$dir_presence, dir_abund = tr$dir_abund, dir_mean = tr$dir_mean, dir_logmean = tr$dir_logmean,
                    c14_q = v$q[i], c14_est = v$est[i], c14_choice = ch, c14_ksign = ksign, c14_tlog = if (!is.null(Tm) && "log" %in% colnames(Tm)) Tm[, "log"] else NA,
                    linda_q = lq, linda_est = le, path = v$path, stringsAsFactors = FALSE)
  out <- out[(is.finite(out$c14_q) & out$c14_q <= 0.05) | (is.finite(out$linda_q) & out$linda_q <= 0.05) | out$truth_type != "none", ]
  tmp <- paste0(ck, ".tmp", Sys.getpid()); utils::write.csv(out, tmp, row.names = FALSE); file.rename(tmp, ck); out
}
res <- parallel::mclapply(seq_len(nrow(jobs)), run_job, mc.cores = opt$cores, mc.preschedule = FALSE)
bad <- vapply(res, function(x) inherits(x, "try-error"), logical(1)); if (any(bad)) { cat("errors:\n"); print(unique(unlist(lapply(res[bad], as.character)))) }
D <- do.call(rbind, Filter(is.data.frame, res)); utils::write.csv(D, file.path(odir, "taxa.csv"), row.names = FALSE)

# ---- summary: called true positives, by truth type and by the kernel that made the call ----
sg <- function(a, b) ifelse(is.finite(a) & is.finite(b) & b != 0, sign(a) == sign(b), NA)
tp <- D[D$truth_type != "none" & is.finite(D$c14_q) & D$c14_q <= 0.05, ]
tp$truth_dir <- ifelse(tp$truth_type == "prevalence", tp$dir_presence, ifelse(tp$truth_type == "opp", tp$dir_mean, tp$dir_abund))
tp$est_ok <- sg(tp$c14_est, tp$truth_dir); tp$kern_ok <- sg(tp$c14_ksign, tp$truth_dir); tp$est_kern <- sg(tp$c14_est, tp$c14_ksign)
f3 <- function(x) sprintf("%.3f", mean(x, na.rm = TRUE))
sink(file.path(odir, "summary.txt"), split = TRUE)
cat(sprintf("sign diagnostics '%s' (%s): %d cells; c14 true positives %d\n\n", opt$label, format(Sys.time(), "%F %T"), length(unique(paste(D$setting, D$template, D$rep))), nrow(tp)))
cat("c14 true positives: share whose reported estimate has the true sign (est), whose calling kernel's statistic has it (kernel),\nand where estimate and calling kernel agree (est=kernel). Truth: change type's own direction; opposite taxa against mean abundance.\n\n")
by <- function(g) { a <- do.call(rbind, lapply(split(tp, tp[[g]]), function(x) data.frame(group = x[[g]][1], n = nrow(x), est = f3(x$est_ok), kernel = f3(x$kern_ok), est_kernel = f3(x$est_kern))))
  names(a)[1] <- g; print(a, row.names = FALSE); cat("\n") }
by("truth_type"); by("c14_choice")
cat("By setting:\n"); by("setting")
op <- tp[tp$truth_type == "opp", ]
if (nrow(op)) { cat("Opposite taxa (presence one way, abundance among the present the other): c14's estimate sign against each estimand\n")
  print(data.frame(n = nrow(op), presence = f3(sg(op$c14_est, op$dir_presence)), abund_given_presence = f3(sg(op$c14_est, op$dir_abund)),
                   mean_abundance = f3(sg(op$c14_est, op$dir_mean)), mean_log_small_pseudocount = f3(sg(op$c14_est, op$dir_logmean)),
                   mean_vs_logmean_agree = f3(op$dir_mean == op$dir_logmean)), row.names = FALSE); cat("\n") }
lt <- D[D$truth_type != "none" & is.finite(D$linda_q) & D$linda_q <= 0.05, ]
lt$truth_dir <- ifelse(lt$truth_type == "prevalence", lt$dir_presence, ifelse(lt$truth_type == "opp", lt$dir_mean, lt$dir_abund))
cat("LinDA true positives on the same cells (est = its estimate has the true sign):\n")
print(do.call(rbind, lapply(split(lt, lt$truth_type), function(x) data.frame(truth_type = x$truth_type[1], n = nrow(x), est = f3(sg(x$linda_est, x$truth_dir))))), row.names = FALSE)
lo <- lt[lt$truth_type == "opp", ]
if (nrow(lo)) { cat("\nLinDA on opposite taxa:\n"); print(data.frame(n = nrow(lo), presence = f3(sg(lo$linda_est, lo$dir_presence)), abund_given_presence = f3(sg(lo$linda_est, lo$dir_abund)),
  mean_abundance = f3(sg(lo$linda_est, lo$dir_mean)), mean_log_small_pseudocount = f3(sg(lo$linda_est, lo$dir_logmean))), row.names = FALSE) }
sink()
