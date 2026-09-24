# benchmarks/dev — the inner loop of PURSUE R&D

Rules live in the project charter (`claude/rd-charter.md`); the running log in `claude/rd-notebook.md`.

- `devsuite.R` scores candidate DA procedures on the **tuning-pool templates only**
  (`benchmarks/devdata/hmp_tongue.rds`, `twinsuk_stool.rds`), with seeds disjoint from the benchmark's.
- `candidates/*.R` each call `register_candidate(name, fn, notes)`. `fn(counts, meta, formula, tested_term)`
  gets a prevalence-filtered features x samples table and returns `data.frame(feature, p)`.
- `results/<label>/` holds `cells.csv`, `summary.csv` and `scoreboard.csv` for each run.

House and implant settings run in Claude's sandbox. On the server, run the msq/mid half:

    Rscript benchmarks/dev/devsuite.R --sims msq,mid --reps 5 --cores 32 --label <label>
