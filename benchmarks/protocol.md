# Benchmark protocol — v1.0 (as implemented)

Status: **implemented, awaiting freeze**. v0.9 was the design draft; v1.0 records what the
engine in `benchmarks/R/` actually does, with every departure from v0.9 marked **[impl]**.
The freeze is recorded in §13 the day the first evaluation-pool run is submitted; after that,
changes are dated amendments appended in §13 and nothing above §13 is edited. PURSUE 0.2 was
designed and validated on the in-container generator and the tuning templates only
(`design-experiments/`); it has not been run on the evaluation pool.

This document follows ADEMP (Aims, Data-generating mechanisms, Estimands, Methods,
Performance measures) and Weber et al. 2019. The design decisions from
`plan-v0.1` §5 and the findings of `results-cycle1` are incorporated.

---

## 1. Aims

A1. For each differential-abundance method, estimate its **false discovery rate at nominal
    BH q = 0.05 and 0.10**, its **sensitivity**, and its **p-value calibration under the
    null**, across realistic 16S OTU/ASV data-generating mechanisms.

A2. Estimate **how much each method's performance depends on which mechanism generated the
    data** — the method × simulator interaction — and report it as a first-class result.

A3. Estimate, for two-part methods, the **cross-arm shadow**: how often an abundance-only
    effect is reported as a prevalence effect and vice versa.

A4. Estimate **effect-size recovery** (bias, RMSE, CI coverage) where truth is on a known
    scale.

Non-aims: runtime optimisation beyond reporting; multi-omic integration; longitudinal
trajectory modelling beyond a repeated-measures design.

## 2. Estimands and the definition of truth

Every simulated feature carries three truth fields:

- `truth_abs` — whether its **absolute abundance** changes with the exposure (the
  generator's own perturbation flag).
- `truth_rel` — whether its **relative abundance** changes by more than a tolerance
  (|log₂ fold change| > 0.25 in expectation), computed from the generator's realised
  compositions. Null features acquire `truth_rel = 1` when a bloom or a large
  one-directional shift moves the composition.
- `truth_type` ∈ {none, abundance, prevalence, both}: which component the perturbation
  targeted.

Each method is scored against the truth matching its **stated estimand**, and also against
the other, so that estimand mismatch is a visible cost rather than a hidden one. Methods are
classified by estimand in §6.

## 3. Axes

| Axis | Truth | What it uniquely tests | Implementation |
|---|---|---|---|
| A. Multi-simulator synthetic | exact, per §2 | power across regimes; method × simulator interaction | §4 |
| B. Signal implantation in real data | exact (implanted) | prevalence vs abundance attribution; confounding; effect-size recovery | §5 |
| C. Real-data null calibration | exact null | p-value uniformity; empirical FPR under several permutation schemes | §7 |
| D. Biological / experimental truth | partial | sanity against known biology; spike-in false positives; low-complexity communities | §8 |
| E. Cross-cohort replicability | none | direction/significance agreement across splits and studies | §9 |

## 4. Axis A — multi-simulator synthetic

### 4.1 Simulators (levels of the factor `simulator`)

| id | package / version pin | truth scale | notes |
|---|---|---|---|
| `msq` | GUniFrac::SimulateMSeq (≥ 1.8) | absolute | incumbent; stability-preserving construction |
| `sd2` | sparseDOSSA 2 (Bioconductor, pin) | absolute (ZILN) | fits cached per template; slow |
| `mid` | MIDASim (CRAN, pin), parametric mode | relative shift + cascade | fast; LOCOM2's simulator |
| `sps` | SPsimSeq (Bioconductor, pin) | absolute (resampled) | least parametric; **[impl]** implants no effects of its own -- it selects features that genuinely differ between two groups it is shown and reuses their magnitudes, so it needs a template with real structure. `sps_group` in `templates.tsv` names that variable: a biological grouping where one exists, otherwise sequencing centre (real technical differences). On a template with neither, the cell is recorded `unsupported`. Shown a random grouping it silently returns zero differential features, which would score as every method having no true positives; both that and the analogous sparseDOSSA2 name mismatch now abort the cell instead |
| `house` | in-house ZI lognormal–Poisson, `benchmarks/R/simulators/sim_house.R` | absolute | **[impl]** the only generator that expresses every regime factor of §4.3 (depth confounding, covariate confounding, repeated measures, bloom, heteroscedastic groups); parameters fitted per template |

**[impl]** When an external simulator cannot express a regime factor, the cell is recorded as
skipped with the offending factors listed (`<cell>.skipped.json`) rather than silently run at
the reference level; the summary tables therefore have unequal cell counts across simulators
and the interaction model of §4.5 is fitted on the cells that exist. The `msq` wrapper was
executed in development; `mid`, `sd2` and `sps` were written against package documentation
and are verified by `hpc/smoke_test.R` on the cluster before any array job (§13).

A fifth generator (deep generative) is admitted only if it supports conditional generation
with a known perturbation; otherwise excluded. Realism of every (simulator, template) pair is
measured (§4.4) and reported; no simulator is excluded on realism grounds a priori.

### 4.2 Templates

Public 16S datasets, each processed to a genus- or OTU-level count table with ≥ 100 samples
and ≥ 200 features after a 5%-prevalence filter. Two disjoint pools:

- **Tuning pool** (method development may use): `hmp_tongue` (HMP V35 oral, 365 samples),
  `twinsuk_stool` (1024 stool samples). **[impl]** 2 templates, not 3.
- **Evaluation pool**: `hmp_stool` (388, high-complexity gut), `hmp_gingiva` (359, oral),
  `hmp_skin_ear` (810, skin), `hmp_vagina` (442, low-complexity, *Lactobacillus*-dominated),
  `risk_stool` (166 paediatric-CD stool, deep). **[impl]** 5 templates, not ≥ 6: the
  environmental (soil/freshwater) and infant-gut slots are unfilled because no candidate
  with a curl-able count table was found in the sources below; they are the first
  amendments to make (`hpc/README.md` §2 describes the drop-in mechanism, and a template
  added after the freeze is reported separately from the frozen pool).

Sources: MicrobeDS (phyloseq objects of Qiita 1928 / 1939 / 2014), Zenodo 7382814 (CRC
genus tables), Zenodo 6911027 (MicrobiomeBenchmarkData). The registry is
`benchmarks/templates.tsv` (source, subset expression, prevalence and depth filters);
`hpc/download_templates.sh` writes `benchmarks/data/SHA256SUMS`, which is part of the freeze.

### 4.3 Regime design

A **reference regime** plus one-factor-at-a-time sweeps (star design), to keep the grid
interpretable and affordable. Every regime is run on every simulator × evaluation template.

Reference regime: n = 50 per group, m = 500 features, DA fraction 10%, balanced direction,
medium effect (|log₂ FC| = 1 abundance; presence odds ratio 4), mixed signal (⅓ abundance-
only, ⅓ prevalence-only, ⅓ both), no depth confounding, no covariate confounding, binary
exposure, independent samples, no bloom, homoscedastic.

Sweeps (one factor varied, all others at reference):

| factor | levels |
|---|---|
| n per group | 20, 50, 100, 200 |
| m | 200, 500, 1000 |
| DA fraction | 0 (global null), 0.05, 0.10, 0.25, 0.40 |
| direction balance | balanced, 80/20, 100/0 |
| effect size | small (0.5 / OR 2), medium, large (2 / OR 8), graded U(0.25, 2) |
| signal type | abundance-only, prevalence-only, mixed |
| depth confounding | 1×, 2×, 4×, 9× |
| covariate confounding ϕ | 0, 0.4, 0.7 (confounder affects 10% of non-DA features) |
| exposure | binary, continuous (effects linear in a standardised covariate) |
| design | independent; repeated measures (5 visits × n/5 subjects, ICC 0.3) |
| bloom | none; one feature ×20 in the exposed group |
| heteroscedastic-unbalanced | no; 30/70 groups with 3× variance in the small group |

This is 25 non-reference regimes + reference = 26 per simulator × template (`R00`–`R25` in
`benchmarks/regimes.tsv`; the v0.9 count of 34 was an arithmetic slip). Replicates: 20 for
non-null regimes, 50 for the global null. Full Axis A: 5 simulators × 5 evaluation templates
× 550 = 13 750 cells.

### 4.4 Realism gate

For each (simulator, template): generate 200 unperturbed samples, pool with 200 real samples
from the template, train a random forest (features: relative abundances, richness, depth,
Shannon) with 5-fold CV, record AUROC. Also record: per-feature KS distance between real and
simulated relative-abundance marginals (median), and Frobenius distance between real and
simulated Spearman correlation matrices on the 100 most prevalent features. Report all three
per (simulator, template). AUROC > 0.9 flags the pair as low-realism; flagged pairs are kept
in the analysis but down-weighted (weight 0.5) in the summary ranking and shown separately.
**[impl]** `benchmarks/R/engine/realism.R`: classifier is ranger, else randomForest, else a
PC-logistic fallback; the gate is computed and stored per (simulator, template) in
`results/realism/`; the down-weighting is applied at analysis time, not in the engine.

### 4.5 Analysis model

For each metric Y (FDR, TPR, pAUC) per cell:

    Y ~ method × simulator + method × regime + (1 | template)  [+ realism weight]

Fitted as a linear mixed model on logit(Y) (with continuity correction). Report: marginal
means per method; the method × simulator interaction variance component and its 95% CI;
per-method "fragility" = SD of the method's simulator-specific effects. Rank tables per axis
carry bootstrap CIs on ranks (1000 resamples over templates and replicates).
**[impl]** `benchmarks/R/analysis/aggregate.R` fits the model with lme4 on logit(FDR) and
logit(TPR) and reports the interaction likelihood-ratio test and per-method fragility; rank
bootstrap is done in the analysis notebook, not the engine.

## 5. Axis B — signal implantation

SIMBA (Wirbel et al. 2024) or a re-implementation with identical semantics:

1. Choose a real template; draw two random groups of n each (n ∈ {20, 50, 100, 200}).
2. Implant into a fraction f ∈ {0.05, 0.10} of features:
   - abundance: multiply counts in the exposed group by 2^δ, δ ∈ {0.5, 1, 2}, half up/half down;
   - prevalence: move a fraction s ∈ {0.2, 0.4} of non-zero entries from one group's samples
     to zero and, symmetrically, from zero to a typical non-zero value in the other;
   - both.
3. Confounding: implant a second signal into a disjoint feature set correlated with a
   simulated binary confounder at ϕ ∈ {0.4, 0.7}.
4. Continuous exposure: scale counts by 2^(δ·z) for standardised z.
5. Re-draw depth by multinomial resampling so that implanted counts remain integers and the
   depth distribution matches the template.

Truth is exact and on the implanted scale; `truth_type` recorded. 20 replicates per cell,
all evaluation templates. This is the primary axis for A3 and A4.
**[impl]** Re-implementation in `benchmarks/R/simulators/implant.R`; the 14 specs are in
`benchmarks/implant_specs.tsv` (reference; n 20/100; abundance-only, prevalence-only;
small/large/graded effect; f = 0.05; all-up; ϕ = 0.4/0.7; continuous exposure; global null
with 50 replicates). Prevalence implantation is not defined for the continuous exposure.
Full Axis B: 5 templates × 310 = 1 550 cells.

## 6. Methods

Every wrapper: pinned version, published defaults, same prevalence filter (≥ 10% of
samples) applied *before* the wrapper so filtering is not a method difference, same
metadata, wall-time and memory recorded, one fixed toy table as a unit test to catch API
drift. **[impl]** The registry is `benchmarks/R/methods/registry.R` (17 methods, uniform
wrapper signature); the toy-table test is `hpc/smoke_test.R`; memory is R heap after the
call (`gc()`), not peak RSS. A method whose package is absent is recorded `not_installed`
per cell rather than failing the cell; a wrapper error is `failed`; exceeding the per-method
wall-time cap (`--timeout`, default 1 h) is `timeout`.

| method | estimand | arms | version pin |
|---|---|---|---|
| Wilcoxon on TSS | relative | one | base R |
| t-test / LM on log(TSS + 1e-6·min) | relative | one | base R |
| limma on log-TSS | relative | one | limma ≥ 3.58 |
| logistic on presence | prevalence | one | base R |
| fastANCOM | absolute (median) | one | pin |
| ALDEx2 | relative (CLR) | one | pin |
| ANCOM-BC2 | absolute (sampling fraction) | one (+ structural zeros) | pin |
| LinDA | absolute (median/mode) | one | pin |
| LDM | relative | one | pin |
| LOCOM, LOCOM2 | relative (compositional logistic) | one | pin |
| corncob | relative (beta-binomial) | one | pin |
| ZicoSeq | relative (reference frame) | one | pin |
| MaAsLin 3 | absolute (median) | two, reported separately | pin |
| ADAPT | absolute (median reference) | one | pin |
| radEmu / fastEmu | absolute (typical taxon) | one | pin |
| PURSUE 0.2 | absolute (empirical-null centre) | two, reported separately + Cauchy union | this repository |

**[impl]** PURSUE 0.1 is not in the roster (decision of 2026-09-09: slow, and superseded by
design); its results from the original benchmark stand as the historical comparison. The
wrappers for LinDA, ANCOM-BC2, corncob, LDM, LOCOM, MaAsLin 3, ALDEx2 and ZicoSeq follow the
patterns of the original `benchmark scripts/`; LOCOM2, fastANCOM, ADAPT and fastEmu were
written against documentation and are verified by the smoke test.

Two-part methods are scored per arm (prevalence calls against `truth_type ∈ {prevalence,
both}`, abundance calls against `{abundance, both}`) and on their combined output against
any effect.

## 7. Axis C — null calibration on real data

**[impl]** 10 real datasets (the 5 evaluation templates, the 2 Axis E and the 3 Axis D
datasets; 12 when the tuning pool is included), not ≥ 15. For each: (i) full label shuffle,
(ii) shuffle within a stratifying covariate, (iii) shuffle within depth quintiles, (iv)
shuffle within clusters; where the template has no such variable a random one is drawn so
the scheme still exercises the method's handling of the term. 50 shuffles each, at most 200
samples per shuffle. Record per method: ECDF of p at 0.05, 0.10; KS statistic vs uniform;
number of BH q ≤ 0.05 discoveries (expected ≈ 0). Full Axis C: 10 × 4 × 50 = 2 000 cells.

## 8. Axis D — biological and experimental truth

MicrobiomeBenchmarkData datasets:
- `HMP_2012_16S_gingival_V35` (and V13): enrichment of annotated aerobic taxa in
  supragingival vs anaerobic in subgingival plaque — metric: enrichment odds ratio and
  fraction of calls in the expected direction.
- `Ravel_2011_16S_BV`: *Lactobacillus* decrease / BV-associated increase — same metric; this
  is the low-complexity CLR stress test.
- `Stammler_2016_16S_spikein`: three spike-in taxa at constant load - metric: number of
  spike-in taxa called at q <= 0.05 (any call is a false positive), across 50 random binary
  groupings of the samples. **[impl]** This dataset has **17 samples**, not the 394 stated in
  v0.9 (394 is Ravel's count, carried over in error; confirmed against the
  MicrobiomeBenchmarkData vignette, 2026-09-11). Groupings are therefore 8 vs 9, at which
  size few methods call anything, so this is a **weak supporting null check rather than a
  primary Axis D signal** and must not be ranked on. `make_expected.R` hard-fails unless it
  matches exactly three spike-in ids, because an empty spike-in set would score as zero false
  positives for every method.
- QMP datasets with flow-cytometry totals (Vandeputte 2017 and successors, if
  redistributable): compare each method's log fold changes with load-scaled log fold
  changes; report correlation and sign agreement, separately for absolute- and
  relative-estimand methods. **[impl]** Not implemented: no redistributable count + load
  table was found; open amendment.

**[impl]** Implemented for `mbd_gingival_v35` (V35 only), `mbd_ravel_bv` and
`mbd_stammler_spikein` in `benchmarks/R/engine/run_axisDE.R`; the expected-taxon tables are
built from the MicrobiomeBenchmarkData taxonomy files by `benchmarks/expected/make_expected.R`.

## 9. Axis E — replicability

Following Pelto et al. 2025: 5 random 50/50 splits per dataset; run all methods on both
halves; replication% (same direction, q ≤ 0.05 in both), conflict% (opposite direction,
q ≤ 0.05 in both), NHits. **[impl]** 2 datasets, not ≥ 20: `crc_genus` (Baxter + Zackular +
Zeller, control vs CRC, 708 samples — the three studies also give the cross-study pairs) and
`risk_ileum` (RISK terminal-ileum biopsies, CD 324 vs non-IBD 190). Any 16S dataset with a
binary condition and ≥ 100 samples added to `templates.tsv` with `pool = axisE` joins the
axis.

## 10. Performance measures (definitions)

- **FDR_q**: FP / max(R, 1) among calls with q ≤ 0.05 (and 0.10), per dataset; averaged over
  replicates; reported with Monte-Carlo SE.
- **TPR_q**: TP / P at the same thresholds.
- **pAUC**: partial AUROC for FPR ∈ [0, 0.10], normalised to [0, 1].
- **Calibration**: |ECDF(0.05) − 0.05| and KS statistic on null features' p-values.
- **Shadow_prev**: rejection rate of the prevalence arm on `truth_type = abundance` features.
  **Shadow_abund**: rejection rate of the abundance arm on `truth_type = prevalence`.
- **Effect recovery** (Axes A, B; methods reporting an estimate): bias, RMSE and 95% CI
  coverage of the log₂ fold change on the method's stated scale.
- **Runtime**: wall seconds and peak RSS per dataset.
- **Stability**: mean Jaccard of the q ≤ 0.05 call set across 10 bootstrap resamples of the
  samples (Axis B, reference cells only). **[impl]** Not implemented in the engine; a
  post-hoc script over the `B_ref` cells.

**[impl]** `benchmarks/R/engine/metrics.R` computes FDR/TPR at 0.05 and 0.10, pAUC₁₀,
FPR at nominal 0.05, KS on null p-values, shadow (per arm), effect bias / RMSE / 95 %
coverage, runtime, per method × arm × truth scale (`abs`, `rel`).

No composite score. Per-axis rank tables with bootstrap CIs; the interaction analysis of §4.5.

## 11. Seeds, sealing, and the tuning/evaluation split

- Each cell's seed is `hash(axis, simulator, template, regime_id, replicate)` under the
  master seed (`cell_seed()` in `run_cell.R`). **[impl]** The master seed is **1**
  (`PURSUE_MASTER_SEED`, Robert's decision of 2026-09-09), stated openly rather than sealed;
  the protection against tuning-to-the-benchmark is therefore the rule below alone, plus
  the fact that PURSUE 0.2's design was fixed (this commit) before any evaluation-pool run.
- Tuning pool templates are unrestricted, under any seed.
- Method development may consult Axis A/B results on the tuning pool only. Any run on the
  evaluation pool before the final run is logged as an amendment.

## 12. Results contract

One long table per (axis, regime, replicate) — Parquet when `arrow` is installed, gzipped
CSV otherwise **[impl]** — plus a per-cell metrics CSV and a manifest. Columns:

| column | type | meaning |
|---|---|---|
| axis | chr | A–E |
| simulator | chr | `house`, `msq`, `sd2`, `mid`, `sps`, `implant`, `nullperm` |
| template | chr | template id |
| regime_id | chr | as in `benchmarks/regimes.tsv` |
| replicate | int | |
| seed | int | realised seed |
| method | chr | as in `benchmarks/R/methods/registry.R` |
| method_version | chr | |
| arm | chr | `abundance`, `prevalence`, `combined`, `single` |
| feature | chr | |
| truth_abs, truth_rel | int | 0/1 |
| truth_type | chr | none / abundance / prevalence / both |
| truth_effect | num | true log₂ FC (or log OR for prevalence) on the generator's scale, NA if none |
| p, q | num | |
| estimate, se, ci_lo, ci_hi | num | on the method's stated scale, NA if not reported |
| status | chr | ok / filtered / failed / timeout |
| runtime_s, mem_mb | num | per dataset (repeated across features); mem_mb is R heap after the call |

Manifest: git commit of the benchmark code, R version, package lock hash, cluster node, start
and end times, and the protocol version.

## 13. Amendments

- 2026-09-10 — v1.0: implementation record (all **[impl]** notes above). Open items to
  amend when filled: environmental and infant-gut evaluation templates (§4.2); QMP
  load-truth datasets (§8); additional Axis E datasets (§9); stability metric (§10).
- 2026-09-11 — first cluster smoke test (R 4.3.3, Bioconductor 3.18), recorded in
  `hpc/smoke_*.csv`. Consequences: (a) the ANCOM-BC2 wrapper was rewritten — ANCOMBC 2.4
  dropped the `taxa_are_rows` / `meta_data` calling convention the original benchmark used —
  and its sensitivity flag now gates only the call (q), leaving p raw so that calibration and
  pAUC measure the model (deviation declared under §6); (b) the Stammler sample count is
  corrected in §8 and that check demoted to supporting evidence; (c) `hpc/install_packages.R`
  reordered so PURSUE installs first and every later failure is non-fatal. All 12 templates
  load. Simulators verified on the cluster: `house`, `msq`, `sps`; `mid` and `sd2` remain
  unverified pending installation.
- (freeze date, first evaluation-pool submission, smoke-test CSVs committed) — to be entered
  by Robert.
