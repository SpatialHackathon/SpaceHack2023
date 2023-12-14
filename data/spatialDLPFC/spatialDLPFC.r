#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Niklas Mueller-Boetticher; contributed code

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(
    c("-o", "--out_dir"),
    type = "character", default = NULL,
    help = "Output directory to write files to."
  )
)

# TODO adjust description
description <- "Load data for ..."

opt_parser <- OptionParser(
  usage = description,
  option_list = option_list
)
opt <- parse_args(opt_parser)

out_dir <- opt$out_dir

# The folder structure should look like the following
# out_dir
# |_______sample_1  (sample name can be chosen freely)
#         |_____coordinates.tsv
#         |_____features.tsv
#         |_____observations.tsv
#         |_____counts.mtx  (use Matrix::writeMM)
#         |_____labels.tsv  (optional)
#         |_____H_E.(tiff/png/...)  (optional)
#         |_____H_E.json  (optional, required if H_E is provided)
# |_______sample_2
#         | ...
# |_______samples.tsv
# |_______experiment.json
# if additional output files are required write it also to out_dir


## Your code goes here
technology <- "Visium"

suppressPackageStartupMessages(library(spatialLIBD))
suppressPackageStartupMessages(library(magrittr))

write_tsv <- function(df, path) {
  write.table(df, path, sep = "\t", col.names = NA, quote = FALSE)
}

write_SpatialExperiment_to_folder <- function(spe, path, obs_col, assay_name = "counts") {
  dir.create(path, showWarnings = FALSE, recursive = TRUE)

  colData(spe)[obs_col] %>%
    as.data.frame() %>%
    write_tsv(file.path(path, "observations.tsv"))

  rowData(spe) %>%
    as.data.frame() %>%
    write_tsv(file.path(path, "features.tsv"))

  coords <- spatialCoords(spe)
  mode(coords) <- "integer"
  as.data.frame(coords) %>%
    dplyr::rename(x = "pxl_col_in_fullres", y = "pxl_row_in_fullres") %>%
    write_tsv(file.path(path, "coordinates.tsv"))

  assay(spe, assay_name) %>%
    t() %>%
    Matrix::writeMM(file.path(path, "counts.mtx"))
}

spe <- fetch_data("spatialDLPFC_Visium")

keep_cols <- c("sample_id", "subject", "position", "sex", "age", "row", "col")
colData(spe) <- colData(spe)[, keep_cols]

keep_rows <- c("gene_name", "gene_version", "source", "gene_type")
rowData(spe) <- rowData(spe)[, keep_rows]

for (sample in unique(colData(spe)$sample_id)) {
  spe_sample <- spe[, spe$sample_id == sample]
  write_SpatialExperiment_to_folder(
    spe_sample,
    file.path(out_dir, sample),
    obs_col = c("row", "col")
  )
}

samples_df <- colData(spe) %>%
  as.data.frame() %>%
  dplyr::mutate(replicate = NA) %>%
  dplyr::select(patient = subject, sample = sample_id, position, replicate, sex, age) %>%
  dplyr::distinct() %>%
  dplyr::mutate(directory = sample) %>%
  `row.names<-`(NULL)


## Metadata files
row.names(samples_df) <- NULL
write.table(samples_df, file = file.path(out_dir, "samples.tsv"), sep = "\t", col.names = NA, quote = FALSE)

json <- file(file.path(out_dir, "experiment.json"))
writeLines(c(paste0('{"technology": "', technology, '"}')), json)
close(json)
