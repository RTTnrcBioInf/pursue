# Running the benchmark on the HPC

Everything below runs from the repository root on the cluster. Steps 1–4 are one-time
setup; 5–8 are the benchmark. Nothing needs a login to any data service.

```bash
git clone https://github.com/RTTnrcBioInf/pursue.git && cd pursue
```

Everything below runs from **that clone's root** — the directory holding `DESCRIPTION`.
A partial copy of `benchmarks/` and `hpc/` is not enough: PURSUE is installed from this
directory, so without `DESCRIPTION` at the top it cannot install and there is no method under
test.

## The loop

Code goes cluster-ward and results come back through git, so neither side hand-copies files:

| where | what |
|---|---|
| laptop | `git add -A && git commit && git push` |
| cluster | `git pull` → run the step below |
| cluster | `git add -A && git commit && git push` — the verification record travels back |
| laptop | `git pull` |

`.gitignore` is set up for exactly this: `hpc/installed_packages.csv`, `hpc/install_log.txt`,
`hpc/prep_caches.csv`, `hpc/smoke_*.csv`, `hpc/smoke_errors.log` and `results/summary/` are
tracked on purpose, while the per-cell results, the downloaded templates, the simulator caches
and the SLURM logs are not (they are large and reproducible from the seed).

## 1. Environment

```bash
bash hpc/setup_env.sh                 # builds the pursue-bench env, then installs R packages
```

The environment is **always rebuilt from scratch** — an existing `pursue-bench` prefix is
deleted first. Updating one in place is what failed on 2026-09-11: the old env still held
bioconda `bioconductor-*` builds from an R 4.3 attempt, raising `r-base` to 4.5 forced their
rebuilds, and their post-link scripts ran `R CMD INSTALL` against a half-updated R
(`ERROR: loading failed for 'R', 'R.c~'`). conda's own logger then crashed while formatting
that error (`ValueError: unsupported format character 'T'`), hiding the real cause. If the
script cannot delete the prefix it stops and says so — remove it by hand and re-run:

```bash
rm -rf ~/.conda/envs/pursue-bench
```

**The R version is not a detail.** On R 4.3 / Bioconductor 3.18 the benchmark silently loses
three comparators — **LOCOM2** (its dependency `Deriv` now requires R >= 4.5), **MaAsLin 3**
and **ADAPT** (both Bioconductor >= 3.20). LOCOM2 is the direct competitor PURSUE is measured
against and MaAsLin 3 is the closest rival to the two-part framing, so running without them
leaves the headline comparison unanswerable. `environment.yml` pins **R 4.5 / Bioc 3.21**, and
the script warns loudly if it ends up on anything older.

Solver preference is **micromamba > mamba > conda**, and with none present micromamba is
bootstrapped into `~/.local/bin` (a single static binary, no admin rights). That order is
deliberate: plain conda crashed here linking a bioconda package with `ValueError: unsupported
format character 'T'` — a bug in conda's own prefix-replacement path that mamba and micromamba
do not share. "prefix already exists" is handled too: the script updates the existing
environment rather than failing.

`environment.yml` is deliberately minimal — R, a compiler toolchain, and the system libraries
that are painful to build. Every R package comes from `hpc/install_packages.R` via
CRAN/Bioconductor, which resolve against Bioc 3.21. Pulling the `bioconductor-*` builds
through conda instead made the solve large and fragile for no benefit.

`install_packages.R` installs **PURSUE first**, from this checkout with base R, so it lands
even on a node with no outbound network; then every other package, each best-effort. It always
writes `hpc/installed_packages.csv` (package, source, installed, version) and
`hpc/install_log.txt` with the reason for every failure, even when steps fail. A package that
will not install is reported, not fatal: its method or simulator is recorded as
`not_installed` in the results rather than crashing a cell.

If a whole class of packages fails at once, the node has no outbound network — install from a
login node, or ask for a local CRAN/Bioconductor mirror.

## 2. Data

```bash
bash hpc/download_templates.sh          # ~45 MB into benchmarks/data/, writes SHA256SUMS
Rscript benchmarks/expected/make_expected.R   # Axis D annotation tables from the taxonomy files
```

Sources (all verified public, fetched with curl):

| template ids | file | source |
|---|---|---|
| `hmp_stool`, `hmp_tongue`, `hmp_gingiva`, `hmp_skin_ear`, `hmp_vagina` | `MicrobeDS/HMPv35.rda` | MicrobeDS (Battaglia), Qiita 1928 |
| `risk_stool`, `risk_ileum` | `MicrobeDS/RISK_CCFA.rda` | MicrobeDS, Qiita 1939 |
| `twinsuk_stool` | `MicrobeDS/TwinsUK.rda` | MicrobeDS, Qiita 2014 |
| `crc_genus` | `zenodo_crc/genus.csv` + `metadata.csv` | Zenodo 7382814 (Baxter, Zackular, Zeller) |
| `mbd_gingival_v35`, `mbd_ravel_bv`, `mbd_stammler_spikein` | `mbd/*.tsv` | Zenodo 6911027 (MicrobiomeBenchmarkData) |

To add a template of your own (an environmental 16S table, for instance): put a
features × samples TSV named `<name>_count_matrix.tsv` and `<name>_sample_metadata.tsv`
under `benchmarks/data/own/`, add a row to `benchmarks/templates.tsv` with
`format = tsv_features_x_samples`, and it is picked up everywhere.

## 3. Smoke test — do not skip

```bash
Rscript hpc/smoke_test.R
```

Loads every template, runs every simulator on a tiny regime, runs every method on a tiny
dataset, and writes `hpc/smoke_{templates,simulators,methods}.csv` plus
`hpc/smoke_errors.log` with the full messages (the CSVs carry a truncated one-line status,
the package each method needs, and the number of finite p-values returned). External wrappers
written against package documentation rather than executed are where API drift shows up —
the 2026-09-11 run caught exactly that in ANCOM-BC2. Anything reporting `error:` needs its
wrapper fixed before it goes into an array job; commit the three CSVs and
`installed_packages.csv` — together they are the record of what was verified on this cluster.

## 4. Warm the simulator caches — do not skip either

```bash
Rscript hpc/prep_caches.R                     # evaluation pool, mid + sd2
Rscript hpc/prep_caches.R --pool all --timeout 21600
```

MIDASim needs a setup object per (template, feature count) and sparseDOSSA2 a fit per
template; the sparseDOSSA2 fit is minutes to hours on a full template. These are cached under
`cache/`, but the cache has to exist **before** the array jobs start — otherwise all 200
concurrent tasks miss at the same instant and refit the same template in parallel, which is
the most expensive mistake available here. Writes `hpc/prep_caches.csv`.

Run it on one node (a compute node with a long wall time, or an interactive session). It is
idempotent: re-running skips whatever is already cached.

## 5. Task lists

```bash
# pilot (protocol §5.12): 2 simulators x evaluation templates x 7 regimes x 5 reps
Rscript hpc/make_tasklist.R --pool evaluation --simulators house,mid \
        --regimes R00,R01,R03,R09,R13,R16,R19 --max-rep 5 --out hpc/tasks

# full protocol
Rscript hpc/make_tasklist.R --pool evaluation --simulators house,msq,mid,sd2,sps
```

Prints the cell counts. One line per cell; the array index is the line number.

## 6. Submit

```bash
# all at once, in dependency order (realism/caching first, then A; B, C, D/E independent)
bash hpc/submit_all.sh
# or individually
sbatch --array=1-$(wc -l < hpc/tasks/realism.txt)      hpc/slurm/realism.sbatch
sbatch --array=1-$(wc -l < hpc/tasks/axisA.txt)%200    hpc/slurm/axisA.sbatch
```

Edit `hpc/slurm/env.sh` once for the cluster (conda vs module). The sbatch headers assume
2 CPUs / 8 GB / 6 h per Axis A cell; ZicoSeq, LOCOM and ALDEx2-GLM are the slow methods —
if cells time out, raise `--time` or drop those methods for the pilot with
`METHODS=pursue,limma_logtss,lm_logtss,wilcoxon_tss,logistic_presence,linda,maaslin3,ancombc2`.

Seeds: every cell's seed is derived from `PURSUE_MASTER_SEED` (default **1**) and the cell
coordinates, so a re-run reproduces the same data. Set the variable before submitting if
you ever need a second independent draw.

## 7. What comes back

Per cell, under `results/axis*/`:

- `<cell>.features.csv.gz` (or `.parquet` if `arrow` is installed) — the results contract,
  one row per method × arm × feature (protocol §12).
- `<cell>.metrics.csv` — FDR/TPR at q = 0.05 and 0.10, pAUC, calibration, cross-arm
  shadow, effect-size bias/RMSE/coverage, runtime, per method × arm × truth scale.
- `<cell>.manifest.json` — git commit, R version, node, timestamps.
- `<cell>.skipped.json` — when a simulator cannot express a regime (recorded, not an error).

`results/realism/*.csv` holds the realism gate per (simulator, template).

## 8. Aggregate

```bash
Rscript benchmarks/R/analysis/aggregate.R results
```

Writes `results/summary/`: the long metrics table, per-axis summary tables with
Monte-Carlo SEs, and — once at least two simulators have run — the method × simulator
interaction analysis (`interaction.txt`). Hand `results/summary/` back and the analysis
continues from there; the per-cell files are only needed for drill-down.

## Resuming and re-running

Array tasks are idempotent per cell (same seed → same data). To re-run failures, list the
missing cell ids and resubmit only those lines:

```bash
comm -23 <(seq 1 $(wc -l < hpc/tasks/axisA.txt)) \
         <(ls results/axisA/*.manifest.json | sed 's/.*r\([0-9]*\).manifest.json/\1/' | sort -n | uniq) > redo.txt
sbatch --array=$(paste -sd, redo.txt) hpc/slurm/axisA.sbatch
```
