#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Søren Helweg Dam; implemented method

suppressPackageStartupMessages({
    library(optparse)
    library(jsonlite)
    library(MERINGUE)
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

description <- "Spatially informed clustering with igraph"

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

if (!is.na(opt$dim_red)) {
  dimred_file <- opt$dim_red
  dimred <- read.delim(dimred_file, stringsAsFactors = FALSE, row.names = 1)
}
if (!is.na(opt$image)) {
  image_file <- opt$image
}
if (!is.na(opt$config)) {
  config_file <- opt$config
    config <- fromJSON(config_file)
}

if (!is.na(opt$matrix)) {
  matrix_file <- opt$matrix
    #matrix <- as(Matrix::t(Matrix::readMM(matrix_file)), "CsparseMatrix")
}

technology <- opt$technology
n_clusters <- opt$n_clusters

if (technology %in% c("Visium", "ST")){
    pos_file <- opt$observations
    positions <- read.delim(pos_file, stringsAsFactors = FALSE, row.names = 1)
} else {
    pos_file <- opt$coordinates
    positions <- as.matrix(read.delim(pos_file, sep = "\t", row.names = 1))
    positions[,c(1:2)] <- as.numeric(positions[,c(1:2)])
}

if ("selected" %in% colnames(positions)) {
        positions <- positions[as.logical(positions$selected), c(1, 2)]
  }

# Seed
seed <- opt$seed
set.seed(seed)

## Your code goes here
k <- config$k
alpha <- config$alpha
beta <- config$beta

# Resolution optimization
binary_search <- function(
    dimred,
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
    new_louvain <- igraph::cluster_louvain
    formals(new_louvain)$resolution <- mid
    cluster_ids <- do_clustering(dimred, method = new_louvain, ...)
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
  return(cluster_ids)
}

# Spatial neighbors
W <- getSpatialNeighbors(positions, filterDist = 2)

# Clustering
result <- binary_search(dimred, n_clust_target = n_clusters, resolution_boundaries = c(0.1, 2), 
                        do_clustering = getSpatiallyInformedClusters, 
                        # Meringue specific
                        W = W, 
                        k = k,
                        alpha = alpha,
                        beta = beta)

#labels <- getSpatiallyInformedClusters(dimred, W = neighbors, k = k, alpha = alpha, beta = beta)
label_df <- data.frame("label" = result, row.names = rownames(dimred))


## Write output
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

colnames(label_df) <- c("label")
write.table(label_df, file = label_file, sep = "\t", col.names = NA, quote = FALSE)

if (exists("embedding_df")) {
  write.table(embedding_df, file = embedding_file, sep = "\t", col.names = NA, quote = FALSE)
}
