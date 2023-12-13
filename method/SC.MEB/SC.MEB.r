#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Søren Helweg Dam; implemented method

suppressPackageStartupMessages({
    library(optparse)
    library(SingleCellExperiment)
    library(SC.MEB)
})

option_list <- list(
  make_option(
    c("-c", "--coordinates"),
    type = "character", default = NULL,
    help = "Path to coordinates (as tsv)."
  ),
  make_option(
    c("-m", "--matrix"),
    type = "character", default = NA,
    help = "Path to (transformed) counts (as mtx)."
  ),
  make_option(
    c("-f", "--features"),
    type = "character", default = NULL,
    help = "Path to features (as tsv)."
  ),
  make_option(
    c("-o", "--observations"),
    type = "character", default = NULL,
    help = "Path to observations (as tsv)."
  ),
  make_option(
    c("-n", "--neighbors"),
    type = "character", default = NA,
    help = "Path to neighbor definitions. Square matrix (not necessarily symmetric) where each row contains the neighbors of this observation (as mtx)."
  ),
  make_option(
    c("-d", "--out_dir"),
    type = "character", default = NULL,
    help = "Output directory."
  ),
  make_option(
    c("--dim_red"),
    type = "character", default = NA,
    help = "Reduced dimensionality representation (e.g. PCA)."
  ),
  make_option(
    c("--image"),
    type = "character", default = NA,
    help = "Path to H&E staining."
  ),
  make_option(
    c("--n_clusters"),
    type = "integer", default = NULL,
    help = "Number of clusters to return."
  ),
  make_option(
    c("--technology"),
    type = "character", default = NULL,
    help = "The technology of the dataset (Visium, ST, imaging-based)."
  ),
  make_option(
    c("--seed"),
    type = "integer", default = NULL,
    help = "Seed to use for random operations."
  ),
  make_option(
    c("--config"),
    type = "character", default = NA,
    help = "Optional config file (json) used to pass additional parameters."
  )
)

description <- "Spatial clustering with hidden Markov random field using empirical Bayes"

opt_parser <- OptionParser(
  usage = description,
  option_list = option_list
)
opt <- parse_args(opt_parser)

out_dir <- opt$out_dir

# Output files
label_file <- file.path(out_dir, "domains.tsv")
embedding_file <- file.path(out_dir, "embedding.tsv")
# if additional output files are required write it also to out_dir

# Use these filepaths as input ...
coord_file <- opt$coordinates
feature_file <- opt$features
observation_file <- opt$observations

if (!is.na(opt$neighbors)) {
  neighbors_file <- opt$neighbors
}
if (!is.na(opt$matrix)) {
  matrix_file <- opt$matrix
}
if (!is.na(opt$dim_red)) {
  dimred_file <- opt$dim_red
  dimred <- read.delim(dimred_file, stringsAsFactors = FALSE, row.names = 1)
}
if (!is.na(opt$image)) {
  image_file <- opt$image
}
if (!is.na(opt$config)) {
  config_file <- opt$config
}

technology <- opt$technology
n_clusters <- opt$n_clusters

# You can get SpatialExperiment directly
get_SingleCellExperiment <- function(feature_file, observation_file, matrix_file, dimred_file) {
  rowData <- read.delim(feature_file, stringsAsFactors = FALSE, row.names = 1)
  colData <- read.delim(observation_file, stringsAsFactors = FALSE, row.names = 1)

  # Filter features and samples
  if ("selected" %in% colnames(rowData)) {
    rowData <- rowData[as.logical(rowData$selected), ]
  }
  if ("selected" %in% colnames(colData)) {
    colData <- colData[as.logical(colData$selected), ]
  }

  dimRed <- read.delim(dimred_file, stringsAsFactors = FALSE, row.names = 1)
  dimRed <- as.matrix(dimRed[rownames(colData), ])

  sce <- SingleCellExperiment(
    rowData = rowData,
    colData = colData,
    reducedDims = list("dimRed" = dimRed)
  )
  return(sce)
}

sce <- get_SingleCellExperiment(feature_file, observation_file, matrix_file, dimred_file)


# Seed
seed <- opt$seed
set.seed(seed)

## Your code goes here
# The data.frames with observations may contain a column "selected" which you need to use to
# subset and also use to subset coordinates, neighbors, (transformed) count matrix
# I used their neighbor method because I couldn't get the input neighbors to work. Feel free to give it a go.
neighbors <- find_neighbors2(sce, technology)

fit <- SC.MEB(as.matrix(dimred), neighbors, K_set = n_clusters, num_core = 2)
label_df <- data.frame("label" = unlist(fit["x", ]), row.names = rownames(dimred))


## Write output
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

colnames(label_df) <- c("label")
write.table(label_df, file = label_file, sep = "\t", col.names = NA, quote = FALSE)

if (exists("embedding_df")) {
  write.table(embedding_df, file = embedding_file, sep = "\t", col.names = NA, quote = FALSE)
}
