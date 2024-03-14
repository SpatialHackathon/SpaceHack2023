#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Søren Helweg Dam; implemented method

suppressPackageStartupMessages({
    library(optparse)
    library(jsonlite)
    library(SingleCellExperiment)
    library(Matrix)
    library(SpatialExperiment)
    library(Banksy)
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

description <- "BANKSY cluster cells in a feature space"

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

technology <- opt$technology
n_clusters <- opt$n_clusters

# You can get SpatialExperiment directly
get_SpatialExperiment <- function(
    feature_file,
    observation_file,
    coord_file,
    use_features = TRUE,
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
    #assay(spe, "logcounts", withDimnames = FALSE) <- as(Matrix::t(Matrix::readMM(matrix_file)), "CsparseMatrix")
  }

  # Filter features and samples
  if (("selected" %in% colnames(rowData(spe))) && use_features) {
    spe <- spe[as.logical(rowData(spe)$selected), ]
  }
  if (("selected" %in% colnames(colData(spe))) && use_features) {
    spe <- spe[, as.logical(colData(spe)$selected)]
  }

  if (!is.na(reducedDim_file)) {
    dimRed <- read.delim(reducedDim_file, stringsAsFactors = FALSE, row.names = 1)
    reducedDim(spe, reducedDim_name) <- as.matrix(dimRed[colnames(spe), ])
  }
  return(spe)
}


# Load configuration
lambda <- config$lambda
k_geom <- config$k_geom
npcs <- config$npcs
method <- config$method
use_agf <- config$use_agf
assay_name <- "normcounts"
use_features <- config$use_features

# Seed
seed <- opt$seed

# You can use the data as SpatialExperiment
spe <- get_SpatialExperiment(feature_file = feature_file, observation_file = observation_file,
                                    coord_file = coord_file, matrix_file = matrix_file, use_features = use_features)


# Extract proper coordinates
if (technology %in% c("ST", "Visium")){
    coord_names <- c("row", "col")
} else {
    coord_names <- NULL
}

## Your code goes here
set.seed(seed)

if (!use_features) {
    library(Seurat)
    counts <- assay(spe, "counts")

    scale.factor <- median(colSums(counts))
    seu <- CreateSeuratObject(counts = counts)
    seu <- NormalizeData(seu, normalization.method = 'RC', scale.factor = scale.factor, verbose = FALSE)
    seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
    spe <- spe[VariableFeatures(seu), ]
}

# Normalization to mean library size
spe <- scuttle::computeLibraryFactors(spe)
assay(spe, assay_name) <- scuttle::normalizeCounts(spe, log = FALSE)

# Run BANKSY
spe <- Banksy::computeBanksy(spe, assay_name = assay_name, k_geom = k_geom, compute_agf = use_agf, coord_names = coord_names)
spe <- Banksy::runBanksyPCA(spe, lambda = lambda, npcs = npcs, use_agf = use_agf)


# Resolution optimization
binary_search <- function(
    spe,
    do_clustering,
    n_clust_target,
    resolution_boundaries,
    num_rs = 100,
    tolerance = 1e-3,
    ...) {

  # Initialize boundaries
  lb <- rb <- NULL
  n_clust <- -1

  lb <- resolution_boundaries[1]
  rb <- resolution_boundaries[2]

  i <- 0
  while ((rb - lb > tolerance || lb == rb) && i < num_rs) {
    mid <- sqrt(lb * rb)
    message("Resolution: ", mid)
    result <- do_clustering(spe, resolution = mid, ...)
    # Adjust cluster_ids extraction per method
    cluster_ids <- colData(result)[, clusterNames(result)]
    n_clust <- length(unique(cluster_ids))
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
  return(result)
}


# The data.frames with observations may contain a column "selected" which you need to use to
# subset and also use to subset coordinates, neighbors, (transformed) count matrix
#cnames <- colnames(colData(spe))
bank <- binary_search(spe, n_clust_target = n_clusters, resolution_boundaries = c(0.1, 2), 
                        do_clustering = Banksy::clusterBanksy, 
                        # Banksy specific
                        lambda = lambda, 
                        use_pcs = TRUE, 
                        npcs = npcs, 
                        seed = seed, 
                        method = method, 
                        assay_name = assay_name, 
                        use_agf = use_agf)

label_df <- data.frame("label" = colData(bank)[, clusterNames(bank)], row.names=rownames(colData(bank)))  # data.frame with row.names (cell-id/barcode) and 1 column (label)
if (use_agf) embedding_df <- as.data.frame(t(assay(bank, "H1")))  # optional, data.frame with row.names (cell-id/barcode) and n columns


## Write output
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

colnames(label_df) <- c("label")
write.table(label_df, file = label_file, sep = "\t", col.names = NA, quote = FALSE)

if (exists("embedding_df")) {
  write.table(embedding_df, file = embedding_file, sep = "\t", col.names = NA, quote = FALSE)
}
