# Running the benchmark on the HPC

Everything below runs from the repository root on the cluster. Steps 1–4 are one-time
setup; 5–7 are the benchmark. Nothing needs a login to any data service.

```bash
git clone https://github.com/RTTnrcBioInf/pursue.git && cd pursue
```

## 1. Environment

Pick whichever the cluster supports.

```bash
# A) conda / mamba (self-contained; recommended)
bash hpc/setup_env.sh conda            # creates env "pursue-bench", then installs R packages
conda activate pursue-bench

# B) cluster R module (>= 4.3) + a user library
bash hpc/setup_env.sh module R/4.3.3   # second argument = the module name on your cluster
export R_LIBS_USER=$HOME/R/pursue-bench-lib
```

Both routes end with `hpc/install_packages.R`. It installs **PURSUE first**, from this
checkout with base R, so it lands even on a node with no outbound network (it needs only
limma, sandwich, parallel); then every CRAN, Bioconductor and GitHub package the benchmark
can use, each best-effort. It always writes `hpc/installed_packages.csv` — package, where it
comes from, installed, version — even when steps fail. A package that will not install is
reported, not fatal: its method or simulator is recorded as `not_installed` in the results
rather than crashing a cell. Re-run the script after fixing anything.

Route B builds `$HOME/R/pursue-bench-lib` and puts it **in front of** your existing library
rather than replacing it, so packages you already have on the cluster stay visible. Export
the same `R_LIBS_USER` before submitting (`hpc/slurm/env.sh` does it for you).

If a whole class of packages fails at once, the node has no outbound network — install from a
login node, or ask for a local CRAN/Bioconductor mirror. Note that MaAsLin 3 is only in
Bioconductor 3.20+; on an older R it has to come from GitHub (`biobakery/maaslin3`), which
the script already attempts.

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

## 4. Task lists

```bash
# pilot (protocol 5.12): 2 simulators x evaluation templates x 7 regimes x 5 reps
Rscript hpc/make_tasklist.R --pool evaluation --simulators house,mid \
        --regimes R00,R01,R03,R09,R13,R16,R19 --max-rep 5 --out hpc/tasks

# full protocol
Rscript hpc/make_tasklist.R --pool evaluation --simulators house,msq,mid,sd2,sps
```

Prints the cell counts. One line per cell; the array index is the line number.

## 5. Submit

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

## 6. What comes back

Per cell, under `results/axis*/`:

- `<cell>.features.csv.gz` (or `.parquet` if `arrow` is installed) — the results contract,
  one row per method × arm × feature (protocol §12).
- `<cell>.metrics.csv` — FDR/TPR at q = 0.05 and 0.10, pAUC, calibration, cross-arm
  shadow, effect-size bias/RMSE/coverage, runtime, per method × arm × truth scale.
- `<cell>.manifest.json` — git commit, R version, node, timestamps.
- `<cell>.skipped.json` — when a simulator cannot express a regime (recorded, not an error).

`results/realism/*.csv` holds the realism gate per (simulator, template).

## 7. Aggregate

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
