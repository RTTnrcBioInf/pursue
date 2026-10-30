#!/usr/bin/env bash
# run_fetch_conf3.sh -- download candidate data for the NEXT confirmation (rd-charter, data discipline, 2026-10-30: reserved,
# pristine; format and design checks only until that confirmation is pre-registered). Download, checksum and list only:
# nothing is parsed or analysed. Files go to benchmarks/data/conf3/<source>/ (git-ignored); the manifest (status, file,
# bytes, md5, URL) is pushed as benchmarks/dev/results/conf3_fetch/manifest.txt. Re-runnable: files already present are kept.
#   - MicrobiomeDB curated infant 16S series (microbiomeDB/microbiomeData @ ced9da6, the commit GEMS1 came from; same
#     file layout, so the GEMS1 loader nearly covers them): repeated measures with exposures that vary within a child.
#   - IBDMDB / HMP2 16S (ibdmdb.org "Products", 2018 release): merged taxonomic profiles and the sample metadata.
# Already on the server and also reserved (no download): HMP V35 tonsils / hard palate / saliva / anterior nares (MicrobeDS
# HMPv35 + devdata/hmp16s_subjects.csv), Morton's ASD sibling pairs (conf2/asd_morton), MicrobiomeHD ob_gordon_2008_v2,
# ob_zupancic, edd_singh (conf2/microbiomehd), Nearing 2022's unused studies (conf2/nearing2022).
# Network and ~120 MB of disk only.
#   git pull && mkdir -p logs && (nohup bash hpc/run_fetch_conf3.sh > logs/run_fetch_conf3.log 2>&1 &)
set -o pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"
D=benchmarks/data/conf3; R=benchmarks/dev/results/conf3_fetch; mkdir -p "$D" "$R" logs
M="$R/manifest.txt"; : > "$M.tmp"
note() { echo "$*" | tee -a "$M.tmp"; }
get() {  # get <source> <file name> <url>
  local src="$1" f="$2" url="$3" out="$D/$1/$2" st
  mkdir -p "$(dirname "$out")"
  if [ -s "$out" ]; then st=kept
  elif curl -fsSL --retry 3 --retry-delay 5 -m 1800 -o "$out.part" "$url" && mv "$out.part" "$out"; then st=fetched
  else rm -f "$out.part"; note "FAIL     $src/$f  $url"; return 0; fi
  note "$(printf '%-8s %-72s %10s  %s  %s' "$st" "$src/$f" "$(stat -c %s "$out")" "$(md5sum "$out" | cut -d' ' -f1)" "$url")"
}
note "conf3 candidate data -- $(date '+%F %T')"; note ""

# --- MicrobiomeDB infant series: assay (genus counts), participant, sample and repeated-measure tables, ontology ---
G=https://raw.githubusercontent.com/microbiomeDB/microbiomeData/ced9da6f012842596ffa411b473e2781f296f384/data-raw
mdb() {  # mdb <folder> <file prefix> <16S region>
  local fo="$1" p="$2" rg="$3"
  get "microbiomedb/$fo" "${p}_16S_rRNA_${rg}_assay.txt" "$G/$fo/assayData/${p}_16S_rRNA_%28${rg}%29_assay.txt"
  for s in Participant Sample participant_repeated_measure; do get "microbiomedb/$fo" "${p}_$s.txt" "$G/$fo/sampleMetadata/${p}_$s.txt"; done
  get "microbiomedb/$fo" "${p}_OntologyMetadata.txt" "$G/$fo/ontology/${p}_OntologyMetadata.txt"
}
mdb DiabImmune     DiabImmune            V4
mdb Bangladesh     Bangladesh_healthy_5yr V4
mdb DailyBaby      DailyBaby             V4
mdb ECAM           EMP_10249_ECAM        V4
mdb EcoCF          EcoCF                 V4
mdb MALED_2yr      MALED_healthy         V4
mdb MALED_diarrhea MALED_diarrhea        V4

# --- IBDMDB / HMP2 16S (merged tables; per-sample files exist if these fail) ---
I=https://g-227ca.190ebd.75bc.data.globus.org/ibdmdb
get ibdmdb taxonomic_profiles.tsv.gz  "$I/products/HMP2/16S/2018-01-07/taxonomic_profiles.tsv.gz"
get ibdmdb taxonomic_profiles.biom.gz "$I/products/HMP2/16S/2018-01-07/taxonomic_profiles.biom.gz"
get ibdmdb hmp2_metadata_2018-08-20.csv "$I/metadata/hmp2_metadata_2018-08-20.csv"

note ""; note "total on disk: $(du -sh "$D" | cut -f1)"
mv "$M.tmp" "$M"; cp logs/run_fetch_conf3.log "$R/run.log" 2>/dev/null || true
git add "$R" && git commit -q -m "conf3 candidate data (reserved): download manifest only; data stay on the server" \
  && git pull -q --rebase --autostash && git push -q
echo ">> done $(date '+%F %T')"
