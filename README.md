# PURSUE

**PURSUE** — **P**revalence-abundance **U**nified **R**egression with **S**hrinkage **U**nder an **E**mpirical null —
is a two-part differential abundance analysis (DAA) method for 16S OTU/ASV count tables.
Every feature gets two separately reported tests:

* **Abundance**: a linear model on log relative abundance over the cells where the feature is
  detected, with a centred log-depth covariate, empirical-Bayes variance moderation across
  features (limma) and **empirical-null mixture centring** of the effect estimates, so the
  reported log2 fold change is relative to the community's typical (null) feature rather than
  to a hand-picked reference set.
* **Presence**: a per-feature zero-inflated negative binomial with a depth offset, cross-feature
  dispersion shrinkage and a ridge-penalised presence logit, tested by likelihood ratio on
  *structural* presence — net of the zeros that are merely a consequence of low abundance
  and shallow depth.

Both components accept arbitrary model formulas (covariates, continuous exposures, multi-level
factors). A Cauchy-combined "any effect" p-value is reported alongside the two arms.

Version 0.2 is a ground-up rewrite. The 0.1 design (heuristic reference-taxon selection and
Freedman–Lane permutation inference) is retained in the git history but is no longer part of
the package; the reasoning and the experiments behind the rewrite are in
[`design-experiments/RESULTS.md`](design-experiments/RESULTS.md) and
[`benchmarks/protocol.md`](benchmarks/protocol.md).

## Installation

```r
# limma comes from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("limma")
remotes::install_github("RTTnrcBioInf/pursue")
```

## Minimal run

```r
library(PURSUE)
otu  <- read.csv("examples/otu_table.csv", row.names = 1, check.names = FALSE)   # samples x features
meta <- read.csv("examples/metadata.csv",  row.names = 1)                        # samples x variables

fit <- pursue(otu, meta, formula = ~ brushing_event + confounder, tested_term = "brushing_event")
fit
#> PURSUE two-part differential abundance
#>   100 samples, 500 features (469 tested); tested term: brushing_event
#>   abundance arm: 469 features fitted, prior df 26.6
#>   centre (log2): -0.423; pi0: 0.64
#>   calls at q<=0.05: 8 abundance, 0 presence, 0 both

res <- fit$results
head(res[order(res$any_p), c("feature", "abund_lfc2", "abund_q", "pres_logor", "pres_q", "any_q")])
```

* `otu` — count table, **rows = samples, columns = features** (a matrix or data.frame).
* `meta` — sample metadata with rows aligned to `otu` (matched by row name when present).
* `formula` — every variable in the model, tested and nuisance alike, e.g. `~ group + age + batch`.
* `tested_term` — the term of `formula` to test. Everything else is adjusted for. Continuous
  terms and factors with more than two levels give multi-df tests.

## Results table

One row per feature of the input, in input order. Features below `min_prevalence` are kept with
`NA` statistics and a status explaining why.

| Column | Description |
|---|---|
| `feature` | feature name |
| `prevalence` | fraction of samples in which the feature is detected |
| `abund_lfc2`, `abund_se2` | centred log2 fold change in relative abundance among detected cells, and its SE (includes the uncertainty of the centre) |
| `abund_ci2_lo`, `abund_ci2_hi` | 95 % interval for `abund_lfc2` |
| `abund_t`, `abund_p`, `abund_q` | moderated t statistic, p-value and BH q-value of the abundance arm |
| `pres_logor`, `pres_se`, `pres_ci_lo`, `pres_ci_hi` | change in log odds of structural presence (binary tested term) and its SE / interval |
| `pres_lrt`, `pres_p`, `pres_q` | likelihood-ratio statistic, p-value and BH q-value of the presence arm |
| `pres_theta` | shrunken negative-binomial dispersion used for the feature |
| `any_p`, `any_q` | Cauchy combination of the two arms: "differs in abundance or in presence" |
| `abund_sig`, `pres_sig` | calls at `q_alpha` (default 0.05) |
| `abund_status`, `pres_status` | `ok`, or why the arm was not fitted for this feature |

Two numbers in `fit$centre` deserve a look before interpreting fold changes:

* `delta` — the empirical-null centre (log2). It is the community-wide shift that has been
  subtracted from every abundance effect. A large `|delta|` means the exposure moved the
  bulk of the community, and the fold changes are *relative to that bulk*.
* `pi0` — the estimated fraction of null features. When `pi0` is small (say < 0.5) the null
  cluster is poorly identified and the centre should be treated with caution; the
  compositional identifiability problem has not gone away, it is only being made explicit.

The abundance estimand equals an absolute-abundance fold change only when the null features
form the modal cluster of effects. That assumption is the same one every compositional DAA
method makes in one form or another; PURSUE reports the quantities needed to judge it.

## Parameters you may want to tune

* `min_prevalence` — features detected in fewer than this fraction of samples are not tested
  (default 0.10). `min_nonzero` — minimum detected cells for the abundance arm (default 8).
* `depth_adjust` — centred log-depth as an abundance-arm covariate (default `TRUE`). This is
  what protects the abundance arm from depth confounding; the presence arm always uses depth
  as an offset.
* `robust_se` — HC3 sandwich standard errors instead of moderated variances (default `FALSE`).
  Use for strongly unbalanced designs with heteroscedastic groups; `fit$diagnostics$heteroscedasticity_hint`
  flags when this may matter.
* `center` — empirical-null centring (default `TRUE`). Turning it off reports raw
  log relative abundance effects.
* `d0_prior`, `tau_ridge` — prior weight for dispersion shrinkage and ridge SD in the
  presence arm (defaults 10 and 5). Neither needs changing for ordinary 16S tables.
* `combine` — report the Cauchy-combined `any_p` / `any_q` (default `TRUE`).
* `n_cores` — forked parallelism over features (default 1). A 100 × 500 table takes about
  25 s on one core; runtime is linear in the number of features.

The full list is in `?pursue`.

## Scope

PURSUE 0.2 targets 16S OTU/ASV (and genus-aggregated) count tables. Gene-family / pathway
tables from shotgun data and long-read amplicon tables can be run through it, but the presence
arm is calibrated for the sparsity and dispersion structure of amplicon tables and the
benchmark below does not yet cover those feature types; see `benchmarks/protocol.md` §1 for
what would have to change.

## Benchmark

`benchmarks/` holds the five-axis evaluation protocol (multi-simulator synthetic truth, signal
implantation into real templates, real-data null calibration, biological truth, cross-cohort
replicability) and an engine that runs 17 DAA methods, PURSUE included, under it.
`hpc/README.md` explains how to set it up and run it on a SLURM cluster.

## Repository layout

```
R/, man/, tests/, DESCRIPTION     the PURSUE package
examples/                         small example table + metadata
benchmarks/                       protocol, engine (R/), template / regime / implant registries
hpc/                              environment, data download, smoke test, SLURM scripts
design-experiments/               the in-container experiments that settled the 0.2 design
benchmark scripts/                the original (pre-0.2) per-method benchmark wrappers, kept for reference
```
