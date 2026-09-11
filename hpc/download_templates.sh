#!/usr/bin/env bash
# Download every template dataset into $PURSUE_DATA_ROOT (default: benchmarks/data).
# All sources are public and fetched with curl; nothing needs a login.
#   MicrobeDS (Battaglia): phyloseq objects for HMP V35, RISK-CCFA, TwinsUK (GitHub raw)
#   Zenodo 7382814: three-study CRC genus table (Baxter, Zackular, Zeller) + metadata
#   Zenodo 6911027: MicrobiomeBenchmarkData count matrices + metadata (Axis D)
set -euo pipefail
ROOT="${PURSUE_DATA_ROOT:-$(cd "$(dirname "$0")/.." && pwd)/benchmarks/data}"
mkdir -p "$ROOT/MicrobeDS" "$ROOT/zenodo_crc" "$ROOT/mbd"
get() { local url="$1" out="$2"; if [ -s "$out" ]; then echo "  have $out"; else echo "  GET  $out"; curl -fsSL --retry 5 --retry-delay 5 -o "$out" "$url"; fi; }

echo "== MicrobeDS (phyloseq .rda) =="
for ds in HMPv35 RISK_CCFA TwinsUK; do
  get "https://raw.githubusercontent.com/twbattaglia/MicrobeDS/master/data/$ds.rda" "$ROOT/MicrobeDS/$ds.rda"
done

echo "== Zenodo 7382814: CRC 16S genus table =="
for f in genus.csv metadata.csv; do
  get "https://zenodo.org/records/7382814/files/$f?download=1" "$ROOT/zenodo_crc/$f"
done

echo "== Zenodo 6911027: MicrobiomeBenchmarkData =="
for ds in HMP_2012_16S_gingival_V35 Ravel_2011_16S_BV Stammler_2016_16S_spikein; do
  for part in count_matrix sample_metadata taxonomy_table; do
    get "https://zenodo.org/records/6911027/files/${ds}_${part}.tsv?download=1" "$ROOT/mbd/${ds}_${part}.tsv"
  done
done

echo "== checksums =="
( cd "$ROOT" && find . -type f -name '*.rda' -o -type f -name '*.csv' -o -type f -name '*.tsv' | sort | xargs sha256sum ) > "$ROOT/SHA256SUMS"
echo "done -> $ROOT (see SHA256SUMS)"
