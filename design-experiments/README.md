# Design experiments (cycle 1)

In-container simulation studies that settled the first three design questions for the
PURSUE successor. Each directory holds the simulation script, an analysis script, and the
result tables it produced. `RESULTS.md` is the write-up.

| dir | question | verdict |
|---|---|---|
| `step1_inference/`  | which inference mode for the abundance arm | limma moderated t; HC3 option; permutation diagnostic-only |
| `step2_prevalence/` | how to test prevalence without re-testing abundance or depth | two-stage plug-ins insufficient; joint fit needed (step 4); expected-rarefied response withdrawn |
| `step3_centering/`  | where the compositional reference goes | empirical-null mixture centring with propagated SE; stability reference dropped |
| `step4_joint/`      | can a joint zero-inflated fit be made stable | yes: theta shrunk and fixed + ridge on presence logit; shadow removed, zero failures |

Reproduce: `Rscript <script> <out.csv>` (add `quick` for a smoke run). Requires R >= 4.3 with
limma, pscl, sandwich, parallel. Seeds are fixed inside each script.
