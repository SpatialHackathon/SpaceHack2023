#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created script

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
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
    help = "Reduced dimensionality representation (e.g. PCA).",
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
    help = "The technology of the dataset (Visium, ST, ...)."
  ),
  make_option(
    c("--seed"),
    type = "integer", default = NULL,
    help = "Seed to use for random operations."
  ),
  make_option(
    c("--config"),
    type = "character", default = NA,
    help = "Optional config file used to pass additional parameters."
  )
)

description <- "BayesSpace (https://www.nature.com/articles/s41587-021-00935-2)
Requirements: Visium or ST. PCA file is needed."

opt_parser <- OptionParser(
  usage = description,
  option_list = option_list
)
opt <- parse_args(opt_parser)

out_dir <- opt$out_dir

# Output files
label_file <- file.path(out_dir, "domains.tsv")
embedding_file <- file.path(out_dir, "embedding.tsv")

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
}
if (!is.na(opt$image)) {
  image_file <- opt$image
}
if (!is.na(opt$config)) {
  config_file <- opt$config
  config <- fromJSON(config_file)
}

technology <- opt$technology
n_clusters <- opt$n_clusters

# Seed
seed <- opt$seed
set.seed(seed)

## Your code goes here

suppressPackageStartupMessages(library(BayesSpace))

get_SingleCellExperiment <- function(
    feature_file, 
    observation_file, 
    matrix_file, 
    dimred_file,
    assay_name = "counts",
    reducedDim_name = "reducedDim") {
  rowData <- read.delim(feature_file, stringsAsFactors = FALSE, row.names = 1)
  colData <- read.delim(observation_file, stringsAsFactors = FALSE, row.names = 1)

  coordinates <- read.delim(coord_file, sep = "\t", row.names = 1)
  coordinates <- as.matrix(coordinates[rownames(colData), ])
  coordinates[,c(1:2)] <- as.numeric(coordinates[,c(1:2)])

    sce <- SingleCellExperiment::SingleCellExperiment(
    rowData = rowData, colData = colData, metadata = list("spatialCoords" = coordinates))

  if (!is.na(matrix_file)) {
    assay(sce, assay_name, withDimnames = FALSE) <- as(Matrix::t(Matrix::readMM(matrix_file)), "CsparseMatrix")
    assay(sce, "logcounts", withDimnames = FALSE) <- log1p(as(Matrix::t(Matrix::readMM(matrix_file)), "CsparseMatrix"))
  }

  # Filter features and samples
  if ("selected" %in% colnames(rowData(sce))) {
    sce <- sce[as.logical(rowData(sce)$selected), ]
  }
  if ("selected" %in% colnames(colData(sce))) {
    sce <- sce[, as.logical(colData(sce)$selected)]
  }

  return(sce)
}

sce <- get_SingleCellExperiment(feature_file, observation_file, matrix_file, dimred_file)

# Source: https://www.ezstatconsulting.com/BayesSpace/articles/maynard_DLPFC.html

dec <- scran::modelGeneVar(sce)
top <- scran::getTopHVGs(dec, n = 2000)

sce <- scater::runPCA(sce, subset_row=top)

sce <- spatialPreprocess(
  sce,
  platform = technology,
  skip.PCA = TRUE
)

nPcs = config$nPCs

sce <- spatialCluster(
  sce,
  q = n_clusters,
  d = nPcs,
  platform = technology
)

label_df <- as.data.frame(colData(sce))[c("spatial.cluster")]


## Write output
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

colnames(label_df) <- c("label")
write.table(label_df, file = label_file, sep = "\t", col.names = NA, quote = FALSE)
