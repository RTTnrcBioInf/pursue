# PURSUE

**PURSUE** (**P**revalence-abundance **U**nified **R**eference-normalized **S**ignificance testing **U**sing **E**mpirical permutation) 
is a reference-normalized two-part regression framework for microbial differential abundance analysis (DAA) that 
models prevalence and abundance separately and computes empirical p-values using a Freedman–Lane residual permutation scheme.

## Installation

```r
remotes::install_github("RTTnrcBioInf/pursue")
```

## Minimal run

```r
out <- run_pursue(
  otu = otu_table,
  meta = metadata,
  formula = ~ group + confounder,
  tested_term = "group"
)
```

Where:

`otu` is a feature OTU table (Rows = Samples; Columns = Features)

`meta` is sample Metadata (Rows = Samples; Columns = Variables)

`formula` is where you provide the name of every variable in your statistical analysis, including the tested variables and nuisance variables. 

`tested_term` is where you specify which of the variables present in `formula` are going to be tested. Variables that are in `formula`, but not in `tested_term` will not be tested, but rather adjusted for as nuisance variation. If there is only one variable in both `formula` and `tested_term`, only that variable will be tested, and additional nuisance variation adjustment will not be performed.

You may use the `otu_table.csv` and `metadata` in /examples.

To view the final results table with q-values:
```r
res <- out$results_table
```

## Main results table


```r
head(res)
```

Example output:

| taxon_index | taxon_name | taxon_index_model | prev_empirical_p | abund_empirical_p | union_p_value | union_q_value | rejected_bh | rejected_bh_raw | reason_NA |
|---:|---|---:|---:|---:|---:|---:|---|---|---|
| 1 | OTU_001 | 1 | 0.004 | 0.031 | 0.007 | 0.021 | TRUE | TRUE | NA |
| 2 | OTU_002 | 2 | 0.412 | 0.088 | 0.154 | 0.231 | FALSE | FALSE | NA |
| 3 | OTU_003 | NA | NA | NA | NA | NA | NA | NA | reference |

The front columns summarize the main taxon-level results:

| Column | Description |
|---|---|
| `taxon_index` | Taxon index in the original input feature table. |
| `taxon_name` | Taxon or feature name. |
| `taxon_index_model` | Taxon index among modeled non-reference taxa. Excluded taxa have `NA`. |
| `prev_empirical_p` | Empirical p-value from the prevalence arm. |
| `abund_empirical_p` | Empirical p-value from the abundance arm. |
| `union_p_value` | Combined taxon-level p-value from the prevalence and abundance arms. |
| `union_q_value` | Final BH-adjusted q-value used for the main decision. |
| `rejected_bh` | Final rejection call at the selected `q_alpha`. |
| `rejected_bh_if` | Rejection call after independent filtering, when independent filtering is enabled. |
| `rejected_bh_raw` | Rejection call from raw BH correction before independent filtering. |
| `reason_NA` | Reason why a taxon was not modeled, such as `reference` or `filtered`. Modeled taxa have `NA`. |

Further columns include model status, Wald statistics, and other diagnostics.

## Additional parameters you're likely to tune

- `min_prevalence_filter` — Features below the specified proportion are filtered from the OTU table. Used to remove uninformative taxa. Default: **off**.

- `min_total_count_filter` — Features below the specified total count threshold are filtered from the OTU table. Used to remove uninformative taxa. Default: **off**.

- `n_perm` — Number of permutations used for empirical p-value estimation. Low values produce coarse (discrete) p-values with ties. Increase `n_perm` until p-values stabilize. Primary driver of runtime. Default: **1999**.

- `depth_adjust` — Appends `log(depth)` to the model formula to adjust for sequencing depth differences. Default: **on**. **NOTE** that this is used by reference selection and the abundance arm, not by the prevalence arm. The prevalence arm removes any column named `log_depth` from the formula.

- `parallelize_permutations`, `permutation_n_cores` — Enable and configure parallel execution of permutations to reduce runtime. Default: **off**.

- `keep_diagnostics` — If `TRUE`, returns detailed step-wise diagnostics; if `FALSE`, returns only the results table.

-  `tested_vars` — Character vector specifying the variable or variables to be permuted in the Freedman–Lane scheme. If `NULL`, defaults to `tested_term`.

- `cluster_id` — Optional clustering variable defining correlated sample groups for clustered permutation. Use when samples are not exchangeable at the individual-sample level, such as repeated measures or paired designs. Default: `NULL`.

- `strata` — Optional stratification variable restricting permutations to occur within strata. Used to preserve higher-level structure during permutation. Default: `NULL`.

- `cluster_exposure_mode` — Controls whether permutation is performed at the sample level or cluster level. `"sample"` forces sample-level permutation, and `"cluster"` forces cluster-level permutation. Default: `"sample"`.

- `reference_group_cleanup` — Controls whether reference group-shift cleanup is enabled. `TRUE` by default . Should be disabled when the tested term is not binary.

- `abund_resid_winsorize` — Controls whether residual winsorization in the abundance arm is enabled. `FALSE` by default.

The full list of parameters can be viewed through `?run_pursue()`.
