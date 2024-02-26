#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Søren Helweg Dam; implemented method

suppressPackageStartupMessages({
    library(optparse)
    library(jsonlite)
    library(SpatialExperiment)
    library(BASS)
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

description <- "BASS: Bayesian Analytics for Spatial Segmentation"

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
  matrix <- as(Matrix::readMM(neighbors_file), "CsparseMatrix")
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

get_SpatialExperiment <- function(
    feature_file,
    observation_file,
    coord_file,
    matrix_file = NA,
    reducedDim_file = NA,
    assay_name = "logcounts",
    reducedDim_name = "reducedDim") {
  rowData <- read.delim(feature_file, stringsAsFactors = FALSE, row.names = 1)
  colData <- read.delim(observation_file, stringsAsFactors = FALSE, row.names = 1)

  coordinates <- read.delim(coord_file, sep = "\t", row.names = 1)
  coordinates <- as.matrix(coordinates[rownames(colData), ])
  coordinates[,c(1:2)] <- as.numeric(coordinates[,c(1:2)])
    
  spe <- SpatialExperiment::SpatialExperiment(
    rowData = rowData, colData = colData, spatialCoords = coordinates
  )

  if (!is.na(matrix_file)) {
    assay(spe, assay_name, withDimnames = FALSE) <- as(Matrix::t(Matrix::readMM(matrix_file)), "CsparseMatrix")
    #assay(spe, "logcounts", withDimnames = FALSE) <- log1p(as(Matrix::t(Matrix::readMM(matrix_file)), "CsparseMatrix"))
  }

  # Filter features and samples
  if ("selected" %in% colnames(rowData(spe))) {
    spe <- spe[as.logical(rowData(spe)$selected), ]
  }
  if ("selected" %in% colnames(colData(spe))) {
    spe <- spe[, as.logical(colData(spe)$selected)]
  }

  if (!is.na(reducedDim_file)) {
    dimRed <- read.delim(reducedDim_file, stringsAsFactors = FALSE, row.names = 1)
    reducedDim(spe, reducedDim_name) <- as.matrix(dimRed[colnames(spe), ])
  }
  return(spe)
}


technology <- opt$technology
n_clusters <- opt$n_clusters

# Seed
seed <- opt$seed
set.seed(seed)

## Your code goes here
R <- n_clusters
C <- config$C
init_method <- config$init_method
beta_method <- config$beta_method
geneSelect <- config$geneSelect
scaleFeature <- config$scaleFeature
nPC <- config$nPC
nSE <- config$nSE

# SpatialExperiment
spe <- get_SpatialExperiment(
    feature_file = feature_file,
    observation_file = observation_file,
    coord_file = coord_file,
    matrix_file = matrix_file,
    reducedDim_file = dimred_file
)


# Set up BASS object
#list("experiment" = SummarizedExperiment::assay(spe, "logcounts")),
BASS <- createBASSObject(
    list("experiment" = t(SingleCellExperiment::reducedDim(spe))), 
    list("experiment" = SpatialExperiment::spatialCoords(spe)),
    C = C, R = R,
    beta_method = beta_method, init_method = init_method, 
    burnin = 2000, nsample = 10000)

# Data pre-processing:

#BASS <- BASS.preprocess(
#    BASS, doLogNormalize = FALSE,
#    doBatchCorrect = FALSE,
#    geneSelect = geneSelect, nHVG = nSE, nSE = nSE, doPCA = TRUE, 
#    scaleFeature = as.logical(scaleFeature), nPC = nPC)
BASS@X_run <- t(SingleCellExperiment::reducedDim(spe))
# Run BASS algorithm
BASS <- BASS.run(BASS)

# post-process posterior samples:
# 1.Adjust for label switching with the ECR-1 algorithm
# 2.Summarize the posterior samples to obtain the spatial domain labels
BASS <- BASS.postprocess(BASS)

#clabels <- BASS@results$c # cell type clusters
#zlabels <- BASS@results$z # spatial domain labels
labels <- unlist(BASS@results$z)
label_df <- data.frame("label" = labels, row.names = rownames(BASS@xy[[1]]))


## Write output
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

colnames(label_df) <- c("label")
write.table(label_df, file = label_file, sep = "\t", col.names = NA, quote = FALSE)

if (exists("embedding_df")) {
  write.table(embedding_df, file = embedding_file, sep = "\t", col.names = NA, quote = FALSE)
}
