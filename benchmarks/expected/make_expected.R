#!/usr/bin/env Rscript
# Build the Axis D "expected" tables from the MicrobiomeBenchmarkData taxonomy files
# (downloaded by hpc/download_templates.sh). Run once on the HPC:
#   Rscript benchmarks/expected/make_expected.R
# Produces:
#   expected/gingival_aerobes.tsv   feature, expected_direction (+1 = enriched in supragingival)
#   expected/bv_taxa.tsv            feature, expected_direction (+1 = enriched in BV)
#   expected/stammler_spikein_ids.txt  comma-separated feature ids of the three spike-in taxa
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), ".."))
d <- file.path(Sys.getenv("PURSUE_DATA_ROOT", file.path(root, "data")), "mbd")
rd <- function(f) read.delim(file.path(d, f), stringsAsFactors = FALSE, quote = "", check.names = FALSE)

# 1. gingival: MicrobiomeBenchmarkData annotates taxa as aerobic / anaerobic / facultative
tx <- rd("HMP_2012_16S_gingival_V35_taxonomy_table.tsv")
ann_col <- grep("annotation|oxygen|aerob", colnames(tx), ignore.case = TRUE, value = TRUE)[1]
if (is.na(ann_col)) stop("no annotation column found in gingival taxonomy table; columns: ", paste(colnames(tx), collapse = ", "))
feat_col <- if ("taxon" %in% colnames(tx)) "taxon" else colnames(tx)[1]
ann <- tolower(tx[[ann_col]]); dir <- ifelse(grepl("^aerob", ann), 1L, ifelse(grepl("anaerob", ann), -1L, NA_integer_))
g <- data.frame(feature = tx[[feat_col]], expected_direction = dir)[!is.na(dir), ]
write.table(g, file.path(root, "expected", "gingival_aerobes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("gingival:", nrow(g), "annotated taxa (", sum(g$expected_direction == 1), "aerobic,", sum(g$expected_direction == -1), "anaerobic )\n")

# 2. BV: Lactobacillus decreases; the BV-associated genera increase (Ravel 2011)
tx <- rd("Ravel_2011_16S_BV_taxonomy_table.tsv"); feat_col <- if ("taxon" %in% colnames(tx)) "taxon" else colnames(tx)[1]
gen_col <- grep("genus", colnames(tx), ignore.case = TRUE, value = TRUE)[1]
genus <- if (!is.na(gen_col)) tx[[gen_col]] else tx[[feat_col]]
bv_up <- c("Gardnerella", "Atopobium", "Prevotella", "Sneathia", "Megasphaera", "Mobiluncus", "Dialister", "Aerococcus", "Peptoniphilus", "Eggerthella", "Parvimonas", "Gemella", "BVAB1", "BVAB2", "BVAB3", "Leptotrichia", "Mycoplasma", "Ureaplasma")
dir <- ifelse(grepl("Lactobacillus", genus), -1L, ifelse(sapply(genus, function(z) any(startsWith(z, bv_up))), 1L, NA_integer_))
b <- data.frame(feature = tx[[feat_col]], expected_direction = dir)[!is.na(dir), ]
write.table(b, file.path(root, "expected", "bv_taxa.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("BV:", nrow(b), "annotated taxa\n")

# 3. Stammler spike-ins: Salinibacter ruber, Rhizobium radiobacter, Alicyclobacillus acidiphilus
tx <- rd("Stammler_2016_16S_spikein_taxonomy_table.tsv"); feat_col <- if ("taxon" %in% colnames(tx)) "taxon" else colnames(tx)[1]
hit <- apply(tx, 1, function(r) any(grepl("Salinibacter|Rhizobium|Agrobacterium|Alicyclobacillus", r, ignore.case = TRUE)))
ids <- tx[[feat_col]][hit]
writeLines(paste(ids, collapse = ","), file.path(root, "expected", "stammler_spikein_ids.txt"))
cat("spike-ins:", length(ids), "feature ids:", paste(ids, collapse = ", "), "\n")
