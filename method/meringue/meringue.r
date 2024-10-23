#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Søren Helweg Dam; implemented method

suppressPackageStartupMessages({
    library(optparse)
    library(jsonlite)
    library(MERINGUE)
})
# Get script path
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg_name <- "--file="
script_path <- dirname(sub(file_arg_name, "", initial_options[grep(file_arg_name, initial_options)]))
# Source binary search function
source(file.path(script_path, "../search_res.r"))

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
    c("--n_pcs"),
    type = "integer", default = NULL,
    help = "Number of PCs to use."
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
if (!is.na(opt$neighbors)) {
  neighbors_file <- opt$neighbors
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
}


technology <- opt$technology
n_clusters <- opt$n_clusters

if (technology %in% c("Visium", "ST")){
    pos_file <- opt$observations
    positions <- read.delim(pos_file, stringsAsFactors = FALSE, row.names = 1,
                           numerals = "no.loss")
} else {
    pos_file <- opt$coordinates
    positions <- as.matrix(read.delim(pos_file, sep = "\t", row.names = 1,
                                     numerals = "no.loss"))
    #positions[,c(1:2)] <- as.numeric(positions[,c(1:2)])
    mode(positions) = "numeric"
}
positions <- positions[, c(1:2)]

if ("selected" %in% colnames(positions)) {
        positions <- positions[as.logical(positions$selected), c(1, 2)]
  }

# Load configuration
k <- config$k
alpha <- config$alpha
beta <- config$beta
n_pcs <- config$n_pcs
n_pcs <- ifelse(is.null(opt$n_pcs), n_pcs, opt$n_pcs)
filterDist <- config$filterDist
min.reads <- config$min.reads
min.lib.size <- config$min.lib.size

# Seed
seed <- opt$seed
set.seed(seed)

# Preprocessing
counts <- Matrix::t(Matrix::readMM(matrix_file))
colnames(counts) <- rownames(positions)
counts <- cleanCounts(counts = counts, 
                      min.reads = min.reads, 
                      min.lib.size = min.lib.size, 
                      plot=FALSE,
                      verbose=TRUE)

positions <- positions[colnames(counts),]

counts <- normalizeCounts(counts, log=FALSE, verbose=TRUE)

# PCA
pcs.info <- prcomp(t(log10(as.matrix(counts)+1)), center=TRUE)
pcs <- pcs.info$x[,1:n_pcs]


# Spatial neighbors
W <- getSpatialNeighbors(positions, filterDist = filterDist)

# Clustering
extract_nclust <- function(result){
    length(unique(result))
}
do_clustering <- function(pcs, resolution, ...){
    new_louvain <- igraph::cluster_louvain
    formals(new_louvain)$resolution <- resolution
    getSpatiallyInformedClusters(pcs, method = new_louvain, ...)
}

result <- binary_search(
    pcs, 
    n_clust_target = n_clusters, 
    extract_nclust = extract_nclust,
    do_clustering = do_clustering, 
    # Meringue specific
    W = W, 
    k = k,
    alpha = alpha,
    beta = beta)


label_df <- data.frame("label" = result, row.names = rownames(pcs))

## Write output
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

colnames(label_df) <- c("label")
write.table(label_df, file = label_file, sep = "\t", col.names = NA, quote = FALSE)

if (exists("embedding_df")) {
  write.table(embedding_df, file = embedding_file, sep = "\t", col.names = NA, quote = FALSE)
}
