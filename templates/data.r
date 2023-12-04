#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: ENTER YOUR NAME AND CONTRIBUTION HERE

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
# TODO
# features_df = ...  # data.frame with rownames (gene-id/name) and n columns (?)
# observations_df = ...  # data.frame with rownames (cell-id/barcode) and n columns (?)
# coordinates_df = ...  # data.frame with rownames (cell-id/barcode) and 2/3 columns (x, y, z?)
# counts = ...  # matrix with #observations rows x #features columns
# labels_df = None  # optional, data.frame with rownames (cell-id/barcode) and 1 column (label)
# img = None  # optional
# technology = ...  # "Visium", "ST", "imaging-based"
# samples_df = ...  # data.frame with information on samples. columns: (patient, sample, position, replicate, directory, n_clusters), columns can be NA

# Make sure to use consistent indexes for the data.frames
# i.e. the index (not necessarily the order) of observations and coordinates should match
# But the order of observations and features must match counts (observations x features)

# Example how a sample could be written
write_sample <- function(
    path, sample, coordinates_df, observations_df, features_df, counts, labels_df = NA, img = NA) {
  write_tsv <- function(df, path) {
    write.table(df, path, sep = "\t", col.names = NA, quote = FALSE)
  }

  sample_path <- file.path(path, sample)
  dir.create(sample_path, showWarnings = FALSE, recursive = TRUE)

  write_tsv(coordinates_df, file.path(sample_path, "coordinates.tsv"))
  write_tsv(features_df, file.path(sample_path, "features.tsv"))
  write_tsv(observations_df, file.path(sample_path, "observations.tsv"))

  Matrix::writeMM(counts, file.path(sample_path, "counts.mtx"))

  if (!is.na(labels_df)) {
    colnames(labels_df) <- c("label")
    write_tsv(labels_df, file.path(sample_path, "labels.tsv"))
  }

  if (!is.na(img)) {
    # TODO write to image_file
    # H_E.json must contain the scale
  }
}


## Metadata files
samples_df <- samples_df[c("patient", "sample", "position", "replicate", "directory", "n_clusters")]
row.names(samples_df) <- NULL
write.table(samples_df, file = file.path(out_dir, "samples.tsv"), sep = "\t", col.names = NA, quote = FALSE)

json <- file(file.path(out_dir, "experiment.json"))
writeLines(c(paste0('{"technology": "', technology, '"}')), json)
close(json)
