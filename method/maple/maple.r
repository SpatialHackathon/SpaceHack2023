#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Giorgia Moranzoni, implemented method. 

suppressPackageStartupMessages({
    library(optparse)
    library(SpatialExperiment)
    library(Seurat)
    library(maple)
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

# TODO adjust description
description <- "MAPLE: Bayesian spatial finite mixture models for identification of cell sub-populations in multi-sample spatial transcriptomics experiments"

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
get_SpatialExperiment <- function(
    feature_file,
    observation_file,
    coord_file,
    matrix_file = NA,
    reducedDim_file = NA,
    assay_name = "counts",
    reducedDim_name = "PCs") {
  rowData <- read.delim(feature_file, stringsAsFactors = FALSE, row.names = 1)
  colData <- read.delim(observation_file, stringsAsFactors = FALSE, row.names = 1)

  coordinates <- read.delim(coord_file, sep = "\t", row.names = 1)
  coordinates <- as.matrix(coordinates[rownames(colData), ])
  coordinates[,c(1:2)] <- as.numeric(coordinates[,c(1:2)])

    
    if (!is.na(reducedDim_file)) {
    dimRed <- read.delim(reducedDim_file, stringsAsFactors = FALSE, row.names = 1)
    #reducedDim(spe, reducedDim_name) <- as.matrix(dimRed[colnames(spe), ])
  }
  spe <- SpatialExperiment::SpatialExperiment(
    rowData = rowData, colData = colData, spatialCoords = coordinates,
      reducedDims = list(pca = as.matrix(dimRed[rownames(colData), ]))
  )

  if (!is.na(matrix_file)) {
    assay(spe, "counts", withDimnames = FALSE) <- as(Matrix::t(Matrix::readMM(matrix_file)), "CsparseMatrix")
    assay(spe, "logcounts", withDimnames = FALSE) <- as(Matrix::t(Matrix::readMM(matrix_file)), "CsparseMatrix")
  }

  # Filter features and samples
  if ("selected" %in% colnames(rowData(spe))) {
    spe <- spe[as.logical(rowData(spe)$selected), ]
  }
  if ("selected" %in% colnames(colData(spe))) {
    spe <- spe[, as.logical(colData(spe)$selected)]
  }

  
  return(spe)
}

# Seed
seed <- opt$seed
set.seed(seed)
# TODO if the method requires the seed elsewhere please pass it on
# You can use the data as SpatialExperiment
spe <- get_SpatialExperiment(feature_file = feature_file,observation_file = observation_file,
                                    coord_file = coord_file, reducedDim_file = dimred_file, matrix_file = matrix_file)

## Your code goes here
# TODO
seurat_obj <- as.Seurat(spe)

# Insert image coordinates
seurat_obj@images$image =  new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "coordinates_",
  coordinates = as.data.frame(colData(spe))
)

# Run maple
maple_results <- fit_maple(
  seurat_obj, 
  K = n_clusters, 
  emb = "PCs", 
  covars = "sample_id" #, n_dim = d 
  )

# The data.frames with observations may contain a column "selected" which you need to use to
# subset and also use to subset coordinates, neighbors, (transformed) count matrix
# label_df = ...  # data.frame with row.names (cell-id/barcode) and 1 column (label)
# embedding_df = NULL  # optional, data.frame with row.names (cell-id/barcode) and n columns

# save data
label_df <- data.frame("label" = maple_results$z, row.names=colnames(seurat_obj))
embedding_df <- as.data.frame(maple_results$Y)

## Write output
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

colnames(label_df) <- c("label")
write.table(label_df, file = label_file, sep = "\t", col.names = NA, quote = FALSE)

if (exists("embedding_df")) {
  write.table(embedding_df, file = embedding_file, sep = "\t", col.names = NA, quote = FALSE)
}
