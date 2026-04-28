# PURSUE

**PURSUE** is a reference-normalized two-part regression framework for microbial differential abundance analysis (DAA) that models prevalence and abundance separately and computes empirical p-values using a Freedman–Lane residual permutation scheme.

## Installation

```r
remotes::install_github("RTTnrcBioInf/pursue")
```

## Minimal run

```r
out <- run_pursue(
  otu = otu_table,
  meta = metadata,
  formula = ~ group,
  tested_term = "group"
)
```

Where:

`otu` is a feature OTU table (Rows = Samples; Columns = Features)

`meta` is sample Metadata (Rows = Samples; Columns = Variables)

`formula` is where you provide the name of every variable in your statistical analysis, including the tested variables and nuisance variables. 

`tested_term` is where you specify which of the variables present in `formula` are going to be tested. Variables that are in `formula`, but not in `tested_term` will not be tested, but rather adjusted for as nuisance variation. If there is only one variable in both `formula` and `tested_term`, only that variable will be tested, and additional nuisance variation adjustment will not be performed.


To view the final results table with q-values:
```r
res <- out$results_table
```

## Additional parameters you're likely to tune

- `min_prevalence_filter` — Features below the specified proportion are filtered from the OTU table. Used to remove uninformative taxa. Default: **off**.

- `min_total_count_filter` — Features below the specified total count threshold are filtered from the OTU table. Used to remove uninformative taxa. Default: **off**.

- `n_perm` — Number of permutations used for empirical p-value estimation. Low values produce coarse (discrete) p-values with ties. Increase `n_perm` until p-values stabilize. Primary driver of runtime. Default: **1999**.

- `depth_adjust` — Appends `log(depth)` to the model formula to adjust for sequencing depth differences. Default: **on**.

- `parallelize_permutations`, `permutation_n_cores` — Enable and configure parallel execution of permutations to reduce runtime. Default: **off**.

- `keep_diagnostics` — If `TRUE`, returns detailed step-wise diagnostics; if `FALSE`, returns only the results table.

-  `tested_vars` — Character vector specifying the variable or variables to be permuted in the Freedman–Lane scheme. If `NULL`, defaults to `tested_term`.

- `cluster_id` — Optional clustering variable defining correlated sample groups for clustered permutation. Use when samples are not exchangeable at the individual-sample level, such as repeated measures or paired designs. Default: `NULL`.

- `strata` — Optional stratification variable restricting permutations to occur within strata. Used to preserve higher-level structure during permutation. Default: `NULL`.

- `cluster_exposure_mode` — Controls whether permutation is performed at the sample level or cluster level. `"auto"` detects cluster-level exposure structure automatically, `"sample"` forces sample-level permutation, and `"cluster"` forces cluster-level permutation. Default: `"auto"`.
