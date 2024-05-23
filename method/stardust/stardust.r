#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Søren Helweg Dam, implemented method. 

suppressPackageStartupMessages({
    library(optparse)
    library(jsonlite)
    library(SpatialExperiment)
    library(Seurat)
    library(stardust)
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

# Description
description <- "stardust: SNN + louvain on integrated minkowski distance matrices of PCA and spatial coordinates"

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
  config <- fromJSON(config_file)
}


# You can get SpatialExperiment directly
get_SpatialExperiment <- function(
    feature_file,
    observation_file,
    coord_file,
    matrix_file = NA,
    reducedDim_file = NA,
    assay_name = "counts",
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

# Study configuration
technology <- opt$technology
n_clusters <- opt$n_clusters
method <- config$method
pcaDimensions <- config$npcs

# Seed
seed <- opt$seed
set.seed(seed)

# Resolution optimization
binary_search <- function(
    seurat_obj,
    do_clustering,
    n_clust_target,
    resolution_update = 2,
    resolution_init = 1,
    resolution_boundaries=NULL,
    num_rs = 100,
    tolerance = 1e-3,
    ...) {

  # Initialize boundaries
  lb <- rb <- NULL
  n_clust <- -1

  if (!is.null(resolution_boundaries)){
    lb <- resolution_boundaries[1]
    rb <- resolution_boundaries[2]
  } else {
    res <-  resolution_init
    results <- do_clustering(seurat_obj, res = res, ...)
    # Adjust cluster_ids extraction per method
    n_clust <- length(unique(results@active.ident))
    if (n_clust > n_clust_target) {
      while (n_clust > n_clust_target && res > 1e-5) {
        rb <- res
        res <- res/resolution_update
        results <- do_clustering(seurat_obj, res = res, ...)
        n_clust <- length(unique(results@active.ident))
      }
      lb <- res
    } else if (n_clust < n_clust_target) {
      while (n_clust < n_clust_target) {
        lb <- res 
        res <- res*resolution_update
        results <- do_clustering(seurat_obj, res = res, ...)
        n_clust <- length(unique(results@active.ident))
      }
      rb <- res
    }
    if (n_clust == n_clust_target) {lb = rb = res }
  }

  i <- 0
  while ((rb - lb > tolerance || lb == rb) && i < num_rs) {
    mid <- sqrt(lb * rb)
    message("Resolution: ", mid)
    set.seed(seed)
    results <- do_clustering(seurat_obj, res = mid, ...)
    cluster_ids <- results@active.ident
    # Adjust cluster_ids extraction per method
    n_clust <- length(unique(cluster_ids))
    message("Cluster: ", n_clust)
    if (n_clust == n_clust_target || lb == rb) break
    if (n_clust > n_clust_target) {
      rb <- mid
    } else {
      lb <- mid
    }
    i <- i + 1
  }

  # Warning if target not met
  if (n_clust != n_clust_target) {
    warning(sprintf("Warning: n_clust = %d not found in binary search, return best approximation with res = %f and n_clust = %d. (rb = %f, lb = %f, i = %d)", n_clust_target, mid, n_clust, rb, lb, i))
  }
  return(results)
}
## Your code goes here

# SpatialExperiment
spe <- get_SpatialExperiment(
    feature_file = feature_file,
    observation_file = observation_file,
    coord_file = coord_file,
    matrix_file = matrix_file#,
    #reducedDim_file = dimred_file
)

## Your code goes here
countMatrix <- as.data.frame(SummarizedExperiment::assay(spe, "counts"))
if (technology %in% c("ST", "Visium")){
    spotPositions <- as.data.frame(SummarizedExperiment::colData(spe)[,c("row", "col")])
} else {
    spotPositions <- as.data.frame(SpatialExperiment::spatialCoords(spe)[,c(1,2)])
}

# Clustering
options(future.globals.maxSize = 100 * 1000 * 1024^2)

countMatrix <- countMatrix[,sort(colnames(countMatrix))]
d <- dim(spotPositions)[2]
spotPositions <- spotPositions[,(d-1):d]
spotPositions <- spotPositions[sort(rownames(spotPositions)),]

seuratObject <- Seurat::CreateSeuratObject(countMatrix)
seuratObject <- suppressWarnings(Seurat::SCTransform(seuratObject, assay = "RNA", verbose = FALSE))
seuratObject <- Seurat::RunPCA(seuratObject, assay = "SCT", verbose = FALSE)
if(pcaDimensions <= 2){
  pcaDimensions = 2
}
m <- seuratObject@reductions[["pca"]]@cell.embeddings[,1:pcaDimensions]
distPCA = dist(m,method="minkowski",p=2)
distCoord <- dist(spotPositions,method="minkowski",p=2)
distCoord <- distCoord*(max(distPCA)/max(distCoord))

if (method == "weight"){
  spaceWeight <- config$weight
  distCoord <- distCoord*((max(distPCA)*as.double(spaceWeight))/(max(distCoord)))
} else {
  expr_norm <- (distPCA - min(distPCA)) / (max(distPCA) - min(distPCA))
  distCoord <- (distCoord)*(as.double(as.vector(expr_norm)))
}

finalDistance <- as.matrix(distPCA + distCoord)
neighbors <- suppressWarnings(Seurat::FindNeighbors(finalDistance))
neighbors <- list(neighbors_nn=neighbors$nn,neighbors_snn=neighbors$snn)
seuratObject@graphs <- neighbors

do_clustering <- Seurat::FindClusters

output <- binary_search(
    countMatrix, 
    n_clust_target = n_clusters, 
    do_clustering = do_clustering, 
    # stardust specific
    verbose = FALSE,
    graph.name = "neighbors_snn")

# save data
label_df <- data.frame("label" = output@active.ident, row.names=colnames(output))

## Write output
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

colnames(label_df) <- c("label")
write.table(label_df, file = label_file, sep = "\t", col.names = NA, quote = FALSE)

if (exists("embedding_df")) {
  write.table(embedding_df, file = embedding_file, sep = "\t", col.names = NA, quote = FALSE)
}
