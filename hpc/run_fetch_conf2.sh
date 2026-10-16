#!/usr/bin/env bash
# run_fetch_conf2.sh -- download the candidate real datasets for conf2's axes D / E (chosen 2026-10-15, before the
# confirmation; notebook "New real data for conf2's D / E"). Download, checksum and list only: nothing is parsed or
# analysed, so the data stay pristine until the confirmation picks its set. Files go to benchmarks/data/conf2/<source>/
# (git-ignored); the manifest (file, bytes, md5, URL, status; archive listings) is pushed as
# benchmarks/dev/results/conf2_fetch/manifest.txt. Re-runnable: files already present with the right checksum are kept.
# GitHub files are pinned to the commits checked on 2026-10-15. Network and ~300 MB of disk only (no memory cap needed).
#   git pull && mkdir -p logs && (nohup bash hpc/run_fetch_conf2.sh > logs/run_fetch_conf2.log 2>&1 &)
# Not scriptable: Wagner et al. 2025 (Front Microbiol 16:1481197) supplementary Tables 1, 2, 4, 8, 9 -- the article page
# builds its download links in the browser. Save them by hand into benchmarks/data/conf2/wagner2025/.
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"
D=benchmarks/data/conf2; R=benchmarks/dev/results/conf2_fetch; mkdir -p "$D/wagner2025" "$R" logs
M="$R/manifest.txt"; : > "$M.tmp"
note() { echo "$*" | tee -a "$M.tmp"; }
get() {  # get <source> <file name> <url> [md5]
  local src="$1" f="$2" url="$3" want="${4:-}" out="$D/$1/$2" have st
  mkdir -p "$(dirname "$out")"
  if [ -s "$out" ] && { [ -z "$want" ] || [ "$(md5sum "$out" | cut -d' ' -f1)" = "$want" ]; }; then st=kept
  elif curl -fsSL --retry 3 --retry-delay 5 -m 1800 -o "$out.part" "$url" && mv "$out.part" "$out"; then st=fetched
  else rm -f "$out.part"; note "FAIL  $src/$f  $url"; return 0; fi
  have="$(md5sum "$out" | cut -d' ' -f1)"
  [ -n "$want" ] && [ "$have" != "$want" ] && st="MD5-MISMATCH(expected $want)"
  note "$(printf '%-8s %-60s %10s  %s  %s' "$st" "$src/$f" "$(stat -c %s "$out")" "$have" "$url")"
}
note "conf2 candidate data -- $(date '+%F %T')"; note ""

# --- axis E: MicrobiomeHD (Duvallet et al. 2017), Zenodo 10.5281/zenodo.1146764, CC BY-NC 4.0; md5 from the record ---
Z=https://zenodo.org/records/1146764/files
get microbiomehd cdi_schubert_results.tar.gz     "$Z/cdi_schubert_results.tar.gz?download=1"     4655b22631158c23cf1d0fb702823c56
get microbiomehd edd_singh_results.tar.gz        "$Z/edd_singh_results.tar.gz?download=1"        9e369e1e62cb87a53f56059194e74b99
get microbiomehd par_scheperjans_results.tar.gz  "$Z/par_scheperjans_results.tar.gz?download=1"  b0789854ef54157f4f41255340c880de
get microbiomehd t1d_alkanani_results.tar.gz     "$Z/t1d_alkanani_results.tar.gz?download=1"     556eb450ee0d4c4a435f9167308cf0d0
get microbiomehd ob_zupancic_results.tar.gz      "$Z/ob_zupancic_results.tar.gz?download=1"      e3de99636c5f811765f3adade2c01ce6
get microbiomehd ob_gordon_2008_v2_results.tar.gz "$Z/ob_gordon_2008_v2_results.tar.gz?download=1" 23ebfb3e6c78d2e9794dbe0659813314

# --- axis E: GEMS1 diarrhoea, 4 sites (microbiomeDB/microbiomeData @ ced9da6) ---
G=https://raw.githubusercontent.com/microbiomeDB/microbiomeData/ced9da6f012842596ffa411b473e2781f296f384
get gems1 GEMS1.rda                     "$G/data/GEMS1.rda"
get gems1 GEMS_16S_rRNA_V1-V2_assay.txt "$G/data-raw/GEMS1/assayData/GEMS_16S_rRNA_%28V1-V2%29_assay.txt"
get gems1 GEMS_Participant.txt          "$G/data-raw/GEMS1/sampleMetadata/GEMS_Participant.txt"
get gems1 GEMS_Sample.txt               "$G/data-raw/GEMS1/sampleMetadata/GEMS_Sample.txt"
get gems1 GEMS_OntologyMetadata.txt     "$G/data-raw/GEMS1/ontology/GEMS_OntologyMetadata.txt"

# --- axis E: autism multi-cohort, uniformly processed (Morton et al. 2023; mortonjt/asd_multiomics_analyses @ 2df2e18) ---
A=https://raw.githubusercontent.com/mortonjt/asd_multiomics_analyses/2df2e187d8614a2304d816a6ba2763fc88a009c7/sfari/data/sra/Combined
get asd_morton age_sex_matched-v2.biom     "$A/age_sex_matched-v2.biom"
get asd_morton age_sex_match_metadata.txt  "$A/age_sex_match_metadata.txt"
get asd_morton sibling_matched-v2.biom     "$A/sibling_matched-v2.biom"
get asd_morton sibling_match_metadata.txt  "$A/sibling_match_metadata.txt"
get asd_morton feature-table.biom          "$A/deblur/vsearch/clustered_table/feature-table.biom"

# --- axis E: Nearing et al. 2022, all 38 datasets' tables (figshare 14531724, CC BY 4.0); the Arctic freshwater and
#     office sets are inside ---
get nearing2022 16S_rRNA_Microbiome_Datasets.tar.gz "https://ndownloader.figshare.com/files/27857901"   # a gzipped tar, not a zip

# --- axis D: Vandeputte et al. 2017 (Nature 551:507): deblurred OTU table + flow-cytometry counts of the disease cohort
#     (barakbri/CompositionalAnalysis_CodeBase @ 252e220, taken from the paper's supplement), and the paper's
#     supplementary tables and source data (disease labels) ---
V=https://raw.githubusercontent.com/barakbri/CompositionalAnalysis_CodeBase/252e220510be4f94ee1b8310762f937fa1cbf62c/Crohn
get vandeputte2017 gut_otu_table.RData      "$V/gut_otu_table.RData"
get vandeputte2017 FlowCytometryCounts.csv  "$V/FlowCytometryCounts.csv"
get vandeputte2017 README_bundle.md         "$V/README.md"
N="https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fnature24460/MediaObjects/41586_2017_BFnature24460"
get vandeputte2017 MOESM2_supplementary_tables.xlsx "${N}_MOESM2_ESM.xlsx"
for k in 7 8 9 10 11 12 13 14 15 16 17 18; do get vandeputte2017 "MOESM${k}_source_data.csv" "${N}_MOESM${k}_ESM.csv"; done

# --- axis D: Wagner et al. 2025 -- by hand (see header); listed here if present ---
for f in "$D"/wagner2025/*; do [ -f "$f" ] && note "$(printf '%-8s %-60s %10s  %s  %s' manual "wagner2025/$(basename "$f")" "$(stat -c %s "$f")" "$(md5sum "$f" | cut -d' ' -f1)" "(saved by hand)")"; done

note ""; note "archive contents (file names only):"
for t in "$D"/microbiomehd/*.tar.gz; do [ -f "$t" ] && { note "  $(basename "$t"):"; tar -tzf "$t" 2>/dev/null | sed 's/^/    /' | tee -a "$M.tmp" >/dev/null; }; done
z="$D/nearing2022/16S_rRNA_Microbiome_Datasets.tar.gz"
[ -f "$z" ] && { note "  $(basename "$z") (datasets under Hackathon/Studies):"; tar -tzf "$z" 2>/dev/null | grep "^Hackathon/Studies/" | awk -F/ 'NF > 3 {print $3}' | sort -u | sed 's/^/    /' | tee -a "$M.tmp" >/dev/null; }
note ""; note "total on disk: $(du -sh "$D" | cut -f1)"
mv "$M.tmp" "$M"; cp logs/run_fetch_conf2.log "$R/run.log" 2>/dev/null || true
git add "$R" && git commit -q -m "conf2 candidate real data (D / E): download manifest only; data stay on the server" \
  && git pull -q --rebase --autostash && git push -q
echo ">> done $(date '+%F %T')"
