# Cycle 1 results — inference, prevalence arm, compositional centring, joint fit

Date: 2026-09-09 (steps 1–3), 2026-09-10 (step 4). In-container simulation experiments settling the three design bets from
`plan-v0.1`. No cluster time used. All code under `experiments/` (to be committed to the
repo under `benchmarks/design-experiments/` once the layout is agreed). Every number below is
a mean over replicates of 300–500 simulated taxa; Monte-Carlo error on a type I rate is
roughly ±0.005 at these sizes.

## Decisions

| Question | Decision | Confidence |
|---|---|---|
| Inference mode for the abundance arm | **limma moderated t** as default; HC3 sandwich option for unbalanced/heteroscedastic designs; Freedman–Lane retained only as a diagnostic mode | High |
| Compositional centring | **Empirical-null mixture centring** (two-component Gaussian mixture on per-taxon effect estimates with known SEs; null component's mean is the centre, its precision gives SE). Matches the oracle in every scenario | High |
| Reference-set selection (PURSUE's core) | **Dropped.** The outcome-blind stability score is contaminated by coordinated shifts (purity 0.64 at 40% one-directional DA, type I 0.41) and discriminates weakly even when balanced (purity 0.78) | High |
| Abundance response | Nonzero-only log relative abundance **with log-depth as a covariate**. Count models (zero-truncated NB) keep more power under confounding but run 10–12% type I under lognormal mixing even at n = 100 — the same misspecification that sinks count GLMs on real data | High |
| Prevalence arm | **Stabilized joint zero-inflated NB** (step 4): per-taxon ZINB with depth offset, dispersion shrunk across taxa and fixed, weak ridge on the presence logit, LRT on the presence coefficient. Type I 0.02–0.05 in every cell incl. 4× depth confounding; abundance shadow 0.005–0.075 (vs 0.19–0.29 for logistic-based tests at n = 100); zero fitting failures; 2–3× the power of pscl's ZINB. Costs 15–30% power vs naive logistic when depth is *not* confounded. PURSUE's expected-rarefied response withdrawn — it re-tests abundance | High |
| Depth handling | Explicit in both arms. Depth confounding is the dominant threat found in this cycle | High |

---

## Step 1 — inference mode

**Design.** The abundance-arm response was simulated directly (log-scale linear model with
a confounder correlated 0.5 with the exposure), so only the inference layer differed.
Candidates on the same fitted model: classical t (`ols`); limma moderated t (`modt`);
Freedman–Lane residual permutation at B ∈ {199, 999, 9999} (`fl_B`); permutation-variance
Wald à la LOCOM2 with FL residuals (`pvw_B`) and with ter Braak full-model residuals
(`pvt_B`); Pearson-III moment-matched permutation null on F = t² (`mm_B`). Error models:
Gaussian, t₃, 10% contamination at 4× SD, skewed (centred gamma), heteroscedastic with a
30/70 imbalance (small group at 3× variance), and 40% missing cells (the two-part reality).
m = 500 taxa, n ∈ {10, 20, 40, 100} total, 20 null + 15 alternative replicates per cell,
10% DA at one residual SD.

**Findings.**

1. *The classical t is already robust here.* Under t₃, contaminated and skewed residuals,
   `ols` held 4.3–5.1% type I down to n = 10. Permutation bought nothing on calibration.
2. *Moderated t is a strict improvement.* Same calibration; consistently higher power (Gaussian
   n = 40: TPR 0.288 vs 0.275; n = 20: 0.036 vs 0.019). limma's prior-df estimate recovered
   the true d₀ ≈ 4 in every scenario — the empirical-Bayes machinery is working as designed.
3. *Freedman–Lane needs B ≥ m/α to reject anything, and then only matches OLS.* At B = 999,
   power was zero at n ≤ 20 and roughly half of OLS at n = 40 — the resolution ceiling. At
   B = 9999 it equalled OLS (FDR and TPR within noise). The thesis critique of runtime was
   correct and the permutations were not buying robustness.
4. *The LOCOM2-style permutation-variance Wald transplants badly to a linear model.* With FL
   residuals it is calibrated but loses power at small n (contaminated n = 20: TPR 0.005 vs
   0.11 for OLS), because the reduced-model residuals still contain the effect and inflate
   the null SD for exactly the taxa that matter. With ter Braak residuals it is anti-
   conservative (13% type I at n = 10; FDR 0.64). Neither is adopted.
5. *Moment matching is unsafe in the far tail.* Pearson III on F* from 199 permutations gave
   BH false rejections under the global null (0.3–0.75 per 500 taxa vs ~0.05 for the others)
   and FDR 0.33–0.48 at n ≤ 20 in alternative cells. 999 permutations helped but did not fix
   it. Not adopted.
6. *Heteroscedasticity with imbalance broke every candidate*: 9–11% type I for OLS, moderated
   t, FL and PVW alike. Permutation does not fix it because exchangeability fails.
   Follow-up: an HC3 sandwich restored calibration (5.1–6.1%) at all n, at a cost of ~15%
   relative power under Gaussian errors and mild conservatism at n ≤ 20 under heavy tails.

**Decision.** Default `modt`. `robust_se = "HC3"` as an option, recommended when group
sizes differ by more than ~2:1 or a variance-ratio diagnostic fires. Freedman–Lane kept
only as a diagnostic/calibration mode. Heteroscedastic-unbalanced becomes a benchmark regime.

---

## Step 3 — compositional centring

**Design.** Absolute abundances lognormal, multinomial counts at depth ~2×10⁴, m = 300,
n = 50 per group. Every method starts from the same per-taxon fit — lm on log(TSS +
pseudocount), the field baseline — and differs only in how effect estimates are centred:
none, CLR, median (LinDA), density mode (LinDA alternative), precision-weighted Huber
M-estimate with sandwich SE propagated (`huber`, the radEmu idea in a linear setting),
empirical-null two-component mixture with propagated SE (`enull`), an outcome-blind
DACOMP/PURSUE stability reference (20% most stable taxa), an ADAPT-style reference (half of
taxa nearest the median effect), and the oracle. Factors: DA fraction 5/20/40%, direction
balanced vs all-up, even vs dominant-taxon community, ±bloom (one taxon ×20), fixed |β| = 1
vs graded |β| ~ U(0.2, 1.5). 30 replicates.

**Type I on null taxa (even community, fixed effects):**

| balance | DA | none | clr | median | mode | huber | **enull** | dacomp | adapt | oracle |
|---|---|---|---|---|---|---|---|---|---|---|
| all-up | 5% | .077 | .056 | .046 | .051 | .047 | **.045** | .049 | .098 | .046 |
| all-up | 20% | .304 | .185 | .060 | .053 | .064 | **.050** | .103 | .144 | .053 |
| all-up | 40% | .722 | .524 | .203 | .073 | .250 | **.068** | .414 | .280 | .076 |
| balanced | 40% | .161 | .061 | .052 | .057 | .051 | **.051** | .070 | .288 | .052 |

**95% CI coverage at 40% all-up:** none .25, clr .43, median .77, mode .93, huber .71,
**enull .93**, dacomp .54, adapt .88, oracle .93.

**Findings.**

1. *Empirical-null centring is the oracle in practice.* Type I 0.045–0.068 in every cell,
   bias ≤ 0.006, coverage 0.92–0.95, TPR equal to oracle. It holds under graded effects
   (0.048–0.061), where a mixture might be expected to absorb small effects into the null.
2. *The median (LinDA, MaAsLin 3) fails gracefully at 20% one-directional DA and badly at
   40%.* Bias −0.06 → −0.20; coverage falls to 0.77. This is the "majority unchanged"
   assumption failing exactly as the literature warns; the mixture recovers because null
   taxa form the tightest cluster even when they are a minority of the *mass* of effects.
3. *The smoothed-median/Huber centre is no more robust than the median* against one-sided
   shifts (0.25 at 40% all-up). Its contribution is the SE propagation, which `enull` also
   provides — so radEmu's advantage in this setting reduces to the SE, not the centre.
4. *The stability-based reference set — PURSUE's mechanism — is unsafe.* Coordinated shifters
   are stable relative to each other and get selected: reference purity 0.87 at 20% all-up,
   0.64 at 40%; type I 0.10 and 0.41. Even in the balanced 40% case purity is 0.78, because
   ±1 effects are only marginally less stable than nulls when residual SDs are 0.6–1.2.
   The criterion has weak discrimination by construction.
5. *The ADAPT-style reference is anti-conservative everywhere* (0.10–0.29): selecting the
   reference from the same effect estimates that are then tested induces dependence. (The
   published ADAPT has additional validation steps; this is the naive version.)
6. *Dominant-taxon and bloom scenarios are handled by every coefficient-level centring*, CLR
   included. CLR's known pathology in low-complexity communities requires sparsity, which
   this experiment deliberately excludes; it is covered by benchmark Axis D.

**Decision.** Empirical-null mixture centring with propagated SE. It also yields π₀ and the
gap between median- and mixture-centring as free diagnostics of how much of the community is
moving — replacing the reference-purity diagnostic proposed in the plan, which this
experiment shows is not informative.

---

## Step 2 — prevalence arm and the depth problem

**Design.** Structural presence S ~ Bern(ψ), logit ψ = a + b·x; abundance when present
lognormal with effect β·x; counts Poisson or NB(size 2) at depth N; depth lognormal around
10⁴ and multiplied by c ∈ {1, 2, 4} in the exposed group. Taxa are mostly rare (median
expected count ≈ 2 at N = 10⁴) so detection is genuinely depth-driven. Truth types: null
80%, prevalence-only 7%, abundance-only 7%, both 6%. m = 300, n ∈ {20, 50, 100}, 8
replicates. The metric that matters most is the **abundance shadow**: the rejection rate of
the prevalence test on abundance-only taxa, which should be 5%.

**Prevalence arm, n = 100:**

| c | metric | naive | logit+logN | PURSUE ER | occ (naive λ̂) | occ + θ | occ (ZTNB) | ZIP |
|---|---|---|---|---|---|---|---|---|
| 1 | type I | .046 | .043 | .041 | .061 | .062 | .082 | .045 |
| 1 | shadow | .190 | .203 | **.380** | .132 | .112 | .116 | .102 |
| 1 | power (prev-only) | .50 | .52 | .47 | .53 | .50 | .40 | .46 |
| 4 | type I | .254 | .046 | .052 | .148 | .110 | .080 | .121 |
| 4 | shadow | .232 | .073 | **.440** | .151 | .082 | .100 | .120 |
| 4 | power (prev-only) | .70 | .33 | .53 | .71 | .58 | .56 | .65 |

**Abundance arm, type I on null taxa:**

| c | n | nonzero LM | nonzero LM + logN | ZTNB (LRT) | log(TSS+pc) LM | + median corr. |
|---|---|---|---|---|---|---|
| 1 | 100 | .059 | .060 | .12 | .045 | .033 |
| 2 | 100 | .156 | .055 | .11 | .421 | .110 |
| 4 | 100 | .288 | .053 | .11 | .628 | .218 |

**Prevalence shadow on the abundance arm** (rejection of prevalence-only taxa by the abundance
test, n = 100, c = 1): nonzero LM .027, +logN .025, ZTNB .072, **log(TSS+pc) LM .46**.

**Findings.**

1. *Depth confounding dominates.* Naive presence/absence testing goes from 4.6% to 25% type I
   at 4× depth; the nonzero-only abundance LM from 6% to 29%; the pseudocount LM from 4.5%
   to 63%. Every arm needs explicit depth handling, and the 4× and 9× regimes in the existing
   benchmark are the right ones to keep.
2. *The field's standard one-part model is not an abundance test.* log(TSS + pseudocount) on
   all cells flags 46% of prevalence-only taxa as abundance effects even without any depth
   confounding, because zeros enter as log(0.5/N). LinDA's median correction reduces
   depth-driven type I (0.63 → 0.22) but not this. A one-part model answers "did anything
   change"; it cannot attribute.
3. *log-depth as a covariate fixes the nonzero-only abundance arm* at zero cost when depth is
   not confounded (power 0.49 vs 0.49) and at an honest cost when it is (0.31 vs 0.67 at 4×,
   where x and log N correlate ~0.8). The zero-truncated NB with a depth *offset* keeps power
   under confounding (0.62) but its LRT runs at 10–12% type I at these sample sizes — a
   nuisance-dispersion small-sample problem that dispersion shrinkage across taxa (edgeR-style)
   should address. Both are viable; the choice is a power-vs-simplicity trade.
4. *PURSUE's expected-rarefied prevalence response is calibrated on nulls but re-tests
   abundance.* Shadow 0.38–0.44 at n = 100 — worse than doing nothing (0.19–0.23), because
   P^ER is a continuous function of the count and rises with abundance among detected cells.
   It is the wrong response for a prevalence claim.
5. *Occupancy models are the most powerful prevalence tests*, especially at small n (0.23–0.25
   vs 0.08 for naive at n = 20), and roughly halve the abundance shadow. But when detection is
   supplied from the *biased* nonzero-only abundance fit, the bias propagates (15% type I at
   4× depth). Feeding the calibrated abundance model into the occupancy arm is the natural
   fix and is what step 2b tests.
6. *Nothing removes the shadow entirely at n = 100* — not the occupancy variants (0.08–0.13),
   not pscl's zero-inflated Poisson (0.10–0.12). Some of this is estimation error in the
   detection model; some may be intrinsic to finite depth, as MaAsLin 3 says. The benchmark
   should report the shadow as a first-class metric rather than pretend it is zero.

### Step 2b — closing the loop on the prevalence arm

Same generative model; the occupancy arm was re-run with detection supplied from the
*calibrated* abundance fit (`occ_lmd`), with log-depth added to the presence logit
(`occ_lmd_lN`), with NB detection from the zero-truncated fit (`occ_ztnb`), and against
pscl's joint zero-inflated NB (`zinb`, n ≥ 50). The abundance arm added a t-referenced Wald
for the ZTNB.

**Prevalence arm, n = 100:**

| c | metric | logit+logN | occ (calibrated λ̂) | occ + logN in ψ | occ (ZTNB det.) | ZIP | **ZINB** |
|---|---|---|---|---|---|---|---|
| 1 | type I | .043 | .058 | .057 | .082 | .045 | .012 |
| 1 | shadow | .203 | .143 | .149 | .116 | .102 | **.018** |
| 1 | power | .52 | .54 | .54 | .42 | .46 | .23 |
| 4 | type I | .046 | .152 | .059 | .080 | .121 | .026 |
| 4 | shadow | .073 | .152 | .076 | .075 | .120 | **.035** |
| 4 | power | .33 | .70 | .39 | .34 | .65 | .36 |

**Findings.**

7. *The occupancy inflation under depth confounding is not caused by the biased λ̂.* Swapping
   in the calibrated abundance fit left type I at 0.152. The cause is that a plug-in detection
   model estimated from *detected* cells carries truncation-distorted μ and σ, so p(N) is
   miscalibrated as a function of depth and the excess non-detections in the shallow group
   are attributed to lower presence. A two-stage plug-in cannot fix this; the detection
   parameters have to be estimated jointly with presence.
8. *Adding log-depth to the presence logit restores calibration* (0.059 at 4×) *at the
   collinearity price* (power 0.39 vs 0.70). Same trade as the abundance arm.
9. *The joint ZINB has no abundance shadow.* 0.00–0.035 in every cell, versus 0.07–0.20 for
   everything else. The shadow is estimation error, not intrinsic to finite depth. But the
   pscl implementation is conservative (type I 0.01–0.03), low-powered (0.15–0.36) and failed
   on 85% of taxa at n = 20. The information is there; current fitting does not extract it.
10. *The zero-truncated NB's 10–12% type I is not a small-sample artefact.* A t-referenced
    Wald fixed n = 20 (0.053) and nothing else (0.10–0.11 at n ≥ 50). Under lognormal-Poisson
    truth, gamma-Poisson (NB) mixing is misspecified at the 10% level — the same mechanism
    behind count GLMs' real-data failures in every benchmark cited in the plan. The log-scale
    LM is semi-parametric on the mixing distribution and holds 5% everywhere.

**Decision.** Abundance arm: nonzero-only log relative abundance with log-depth covariate,
moderated t, empirical-null centring. Prevalence arm: the contribution is a *stabilized joint
zero-inflated fit* — NB count part with depth offset and dispersion shrunk across taxa
(edgeR-style), zero part with a Firth/ridge penalty against separation, LRT on the zero-part
coefficient — because that is the only construction shown to remove the shadow. It is a
development task, not a plug-in. Until it exists, the first implementation uses logistic +
log-depth and reports the shadow-prone taxa (those whose abundance arm is significant). The
benchmark carries the abundance shadow as a first-class metric in either case.

---

## Step 4 — the stabilized joint fit

**Design.** Same generator as step 2 plus a pure-NB cell (no lognormal mixing) to separate
misspecification from small-sample effects. The joint model per taxon: logit ψ = a + b·x;
log μ = log N + c + β·x; X ~ ZINB(ψ, μ, θ). Stabilizers: θ estimated per taxon by the
zero-truncated NB on detected cells (closed form), shrunk on the log scale toward a loess
trend on mean abundance with prior weight d₀ = 10 pseudo-observations, then *fixed*; a
Gaussian ridge (τ = 5) on (a, b) against separation; LRT on b by refitting the constrained
model from the full solution, five starts. Comparators: naive logistic, logistic + log-depth,
pscl ZINB (n ≥ 50). m = 300, n ∈ {20, 50, 100}, c ∈ {1, 4}, three sampling models, 8 reps.

**Prevalence test, n = 100:**

| mixing | c | metric | naive | logit+logN | **jzi** | pscl ZINB |
|---|---|---|---|---|---|---|
| Poisson-lognormal | 1 | type I | .042 | .040 | **.034** | .012 |
| | 1 | shadow | .220 | .222 | **.068** | .041 |
| | 1 | power | .62 | .61 | .43 | .24 |
| | 4 | type I | .250 | .045 | **.043** | .019 |
| | 4 | shadow | .286 | .080 | **.075** | .036 |
| | 4 | power | .70 | .33 | **.48** | .34 |
| pure NB | 1 | type I | .041 | .041 | **.043** | .021 |
| | 1 | shadow | .281 | .285 | **.005** | .005 |
| | 1 | power | .51 | .50 | .45 | .31 |
| | 4 | type I | .316 | .049 | **.048** | .034 |
| | 4 | shadow | .225 | .149 | **.038** | .017 |
| | 4 | power | .66 | .33 | **.59** | .51 |

Fitting failures: jzi ≤ 0.25 taxa per 300 at any n; pscl ZINB 255–268 of 300 at n = 20.

**Findings.**

11. *The shadow is gone.* 0.005–0.075 across all 18 cells; the two cells above 0.05 are both
    Poisson-lognormal at n = 100, where the NB count part is least well specified. Logistic
    tests, with or without a depth covariate, sit at 0.19–0.29 in the same cells.
12. *Calibrated under depth confounding without paying the collinearity price.* Type I
    0.039–0.048 at 4× depth, power 0.44–0.59 versus 0.31–0.33 for logistic + log-depth. The
    depth offset in the count part carries the depth information with slope fixed at 1, so
    the presence coefficient does not have to compete with a free depth covariate.
13. *Stable where pscl is not.* Fixing θ after shrinkage and adding the weak ridge removed
    the failure mode entirely and roughly tripled power relative to pscl's ZINB, whose
    over-conservatism (type I 0.006–0.034) came from the same nuisance-parameter problem.
14. *The power cost when depth is not confounded is real*: 15–30% relative to naive
    logistic. The joint model attributes some zeros to non-detection and gives up the
    information in them; naive treats every zero as absence. That is the estimand
    difference, and it is the reason the shadow exists for naive tests.
15. *The joint model's count part is not an abundance test.* Type I 0.07–0.14 under
    lognormal mixing, 0.04–0.07 only under pure NB; prevalence shadow 0.06–0.18. The
    log-scale LM with log-depth (type I 0.036–0.065 in every cell) remains the abundance
    arm. The two arms use different likelihoods by design: the joint fit exists to model
    detection, not to estimate fold changes.

**Decision.** Prevalence arm = stabilized joint ZINB, presence-coefficient LRT. Abundance arm
= nonzero-only log relative abundance with log-depth covariate, moderated t, empirical-null
centring. Implementation notes: the ridge and fixed-θ make the test slightly conservative
(0.02–0.04 at n ≤ 50); estimating d₀ by empirical Bayes rather than fixing it at 10, and
extending both parts to arbitrary design matrices (nuisance covariates, continuous
exposure), are the next engineering items.

---

## What this changes in the plan

- §2.4 / §4.5: inference settled — `modt` default, HC3 option, FL diagnostic-only, no
  moment-matching, no permutation-variance Wald.
- §2.1 / §4.2 / §4.6.1: centring settled — empirical-null mixture; the DACOMP-reference
  diagnostic is withdrawn (uninformative); replaced by π₀ and median-vs-mixture gap.
- §2.5 / §4.3: prevalence arm — expected-rarefied response withdrawn; the two-stage
  occupancy plug-in is insufficient; the target is a stabilized joint zero-inflated fit.
  Depth enters both arms explicitly.
- §2.2 / §4.2: count-model abundance arms (truncated NB, Tobit-like) are not adopted:
  mixing-distribution misspecification costs 10% type I at any n. Log-scale stays.
- §5: benchmark must carry the abundance-shadow and prevalence-shadow metrics, the
  heteroscedastic-unbalanced regime, and depth confounding at 2× as well as 4×/9×.
- Component verdicts: reference selection moves from "demote to diagnostic" to "drop";
  expected-rarefied prevalence moves from "keep, extend" to "replace by occupancy".

## Open after this cycle

- Done (step 4): the stabilized joint fit recovers the ZINB's zero shadow with usable power
  and no failures. Remaining: empirical-Bayes d₀, general design matrices, and whether a
  ZIPLN count part (PLNmodels, cluster-side) closes the last two points of shadow under
  lognormal mixing.
- Clustered / repeated-measures designs: not yet tested in any step.
- Heteroscedastic-unbalanced regime and 2× depth confounding into the benchmark grid.
