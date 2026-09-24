# benchmarks/dev — the inner loop of PURSUE R&D

Rules live in the project charter (`claude/rd-charter.md`); the running log in `claude/rd-notebook.md`.

- `devsuite.R` scores candidate DA procedures on the **tuning-pool templates only**
  (`benchmarks/devdata/hmp_tongue.rds`, `twinsuk_stool.rds`), with seeds disjoint from the benchmark's.
- `candidates/*.R` each call `register_candidate(name, fn, notes)`. `fn(counts, meta, formula, tested_term)`
  gets a prevalence-filtered features x samples table and returns `data.frame(feature, p)`.
- `results/<label>/` holds `cells.csv`, `summary.csv` and `scoreboard.csv` for each run, plus per-cell
  checkpoints in `cells/`: re-running with the same label resumes, losing only cells in flight.
- `score.R` is the charter's scoring (gate, then power), shared by `devsuite.R` and `compare.R`.
- `compare.R --labels a,b` ranks candidates from several runs together. Cells are identical across runs,
  so results merge without re-running; it also reads an in-progress run's checkpoints.

House and implant settings run in Claude's sandbox. On the server, run the msq/mid half:

    Rscript benchmarks/dev/devsuite.R --sims msq,mid --reps 5 --cores 16 --label <label>

(16, not 32: the dev suite has no memory guard, and `twinsuk_stool` is large.)

Diagnosis tools:

- `diag_fp.R --candidates a,b --settings house:R00,implant:B_ref --reps 3` rebuilds dev cells from the
  same seeds and describes each candidate's false positives against the other null features
  (realised relative fold change, detection shift, prevalence, abundance).
- `tools/cell.R`: `source()` it, then `dev_cell("house:R08", "twinsuk_stool", 3)` returns that exact
  cell (filtered counts, meta, formula, truth) with every candidate registered in `.cands`.

Candidate files, by iteration (details and results in the notebook):

| file | iteration | candidates |
|---|---|---|
| `00_baselines.R` | it0 | pursue02, pursue02_presence, logistic, logistic_depth, lm_logtss |
| `05_enull.R` | it3 | shared helpers: `enull2` (estimated null scale), `enull_robust` (π₀ ≥ 0.5) |
| `10_sharedlink.R` | it1, it3 | sharedlink, sharedlink_raw, sharedlink_rc |
| `20_nbglm.R` | it2 | nbglm, nbglm_raw, nbglm_hc3 |
| `30_enull2.R` | it3 | sharedlink_en, nbglm_en, nbglm_hc3_en |
| `40_firth.R` | it3 | firth_depth, sharedlink_firth, sharedlink_firth_raw |
| `50_pmix.R` | it4 | pln, pln_raw, pln_sw, pln_rs, pln_rs_raw, pmix, pmix_raw |
| `60_erd.R` | it5 | erd_lm, erd_wx, erd_qb, erd_qb_raw, erd_qp, erl_lm, erl_lm_raw |
