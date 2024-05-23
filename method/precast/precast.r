#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Søren Helweg Dam; implemented method

suppressPackageStartupMessages({
    library(optparse)
    library(jsonlite)
    library(SingleCellExperiment)
    library(Seurat)
    library(PRECAST)
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

description <- "PRECAST: a probabilistic embedding and clustering with alignment for spatial transcriptomics data integration."

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
    neighbors <- as(Matrix::readMM(neighbors_file), "CsparseMatrix")
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
if(!(technology %in% c("Visium", "ST"))) technology <- "Other_SRT"
n_clusters <- opt$n_clusters

# You can get SpatialExperiment directly
get_SingleCellExperiment <- function(
    feature_file,
    observation_file,
    coord_file,
    matrix_file = NA,
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


# Seed
seed <- opt$seed
set.seed(seed)

# You can use the data as SingleCellExperiment
sce <- get_SingleCellExperiment(feature_file = feature_file, observation_file = observation_file,
                                    coord_file = coord_file, matrix_file = matrix_file)

## Your code goes here
seurat_obj <- as.Seurat(sce)
if (!("row" %in% colnames(seurat_obj@meta.data) & "col" %in% colnames(seurat_obj@meta.data))){
    seurat_obj@meta.data$col <- sce@metadata$spatialCoords[rownames(seurat_obj@meta.data), 1]
    seurat_obj@meta.data$row <- sce@metadata$spatialCoords[rownames(seurat_obj@meta.data), 2]
}

n_var = config$gene_number
method = config$method
maxIter = config$maxIter
postminspots = config$postminspots
postminfeatures = config$postminfeatures

if (nrow(seurat_obj)>n_var){
  PRECASTObj <- CreatePRECASTObject(seuList = list(seurat_obj), 
                                  gene.number=n_var,
                                  selectGenesMethod=method, 
                                  postmin.features = postminfeatures, 
                                  postmin.spots = postminspots)
} else {
  PRECASTObj <- CreatePRECASTObject(seuList = list(seurat_obj), 
                                    customGenelist = rownames(rowData(sce)),
                                  postmin.features = postminfeatures, 
                                  postmin.spots = postminspots)
}

platform <- ifelse(technology %in% c("Visium", "ST"), technology, "Other_SRT")
PRECASTObj <- AddAdjList(PRECASTObj, platform = platform)
PRECASTObj <- AddParSetting(PRECASTObj, Sigma_equal = FALSE, coreNum = 10, maxIter = maxIter, verbose = TRUE)

PRECASTObj <- PRECAST(PRECASTObj, K = n_clusters)
PRECASTObj <- SelectModel(PRECASTObj)

# The data.frames with observations may contain a column "selected" which you need to use to
# subset and also use to subset coordinates, neighbors, (transformed) count matrix
label_df <- data.frame("label" = unlist(PRECASTObj@resList$cluster[1]), row.names=colnames(PRECASTObj@seulist[[1]]))
# data.frame with row.names (cell-id/barcode) and 1 column (label)
embedding_df <- as.data.frame(PRECASTObj@resList$hZ[[1]], row.names=colnames(PRECASTObj@seulist[[1]]))  # optional, data.frame with row.names (cell-id/barcode) and n columns


## Write output
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

colnames(label_df) <- c("label")
write.table(label_df, file = label_file, sep = "\t", col.names = NA, quote = FALSE)

if (exists("embedding_df")) {
  write.table(embedding_df, file = embedding_file, sep = "\t", col.names = NA, quote = FALSE)
}
