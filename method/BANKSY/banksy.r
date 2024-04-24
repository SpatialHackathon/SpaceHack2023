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
    library(Seurat)
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
    c("--n_genes"),
    type = "integer", default = NULL,
    help = "Number of genes to use."
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

# Opt params
technology <- opt$technology
n_clusters <- opt$n_clusters
seed <- opt$seed

# Load configuration
lambda <- config$lambda
k_geom <- config$k_geom
method <- config$method
use_agf <- config$use_agf
n_pcs <- config$n_pcs
n_pcs <- ifelse(is.null(opt$n_pcs), n_pcs, opt$n_pcs)
n_genes <- config$n_genes
n_genes <- ifelse(is.null(opt$n_genes), n_genes, opt$n_genes)
assay_name <- "normcounts"


# You can use the data as SpatialExperiment
spe <- get_SpatialExperiment(
    feature_file = feature_file,
    observation_file = observation_file,
    coord_file = coord_file,
    matrix_file = matrix_file)


# Extract proper coordinates
if (technology %in% c("ST", "Visium")){
    coord_names <- c("row", "col")
} else {
    coord_names <- NULL
}

set.seed(seed)

# Variable features (if given - opt takes priority)
if (nrow(spe) >= n_genes){
    counts <- assay(spe, "counts")
    seu <- CreateSeuratObject(counts = counts)
    seu <- NormalizeData(seu, normalization.method = 'RC', scale.factor = median(colSums(counts)), verbose = FALSE)
    seu <- FindVariableFeatures(seu, nfeatures = n_genes, verbose = FALSE)
    spe <- spe[VariableFeatures(seu), ]
}


# Normalization to mean library size
spe <- scuttle::computeLibraryFactors(spe)
assay(spe, assay_name) <- scuttle::normalizeCounts(spe, log = FALSE)

# Run BANKSY
spe <- Banksy::computeBanksy(spe, assay_name = assay_name, k_geom = k_geom, compute_agf = use_agf, verbose = FALSE, coord_names = coord_names)
spe <- Banksy::runBanksyPCA(spe, lambda = lambda, npcs = n_pcs, use_agf = use_agf)


# Resolution optimization
extract_nclust <- function(result){
    length(unique(colData(result)[, clusterNames(result)]))
}

result <- binary_search(
    spe, 
    n_clust_target = n_clusters, 
    extract_nclust = extract_nclust,
    do_clustering = Banksy::clusterBanksy, 
                        # Banksy specific
                        lambda = lambda, 
                        use_pcs = TRUE, 
                        seed = seed, 
                        method = method,
                        assay_name = assay_name, 
                        use_agf = use_agf)

label_df <- data.frame("label" = colData(result)[, clusterNames(result)], row.names=rownames(colData(result)))  # data.frame with row.names (cell-id/barcode) and 1 column (label)
if (use_agf) embedding_df <- as.data.frame(t(assay(result, "H1")))  # optional, data.frame with row.names (cell-id/barcode) and n columns


## Write output
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

colnames(label_df) <- c("label")
write.table(label_df, file = label_file, sep = "\t", col.names = NA, quote = FALSE)

if (exists("embedding_df")) {
  write.table(embedding_df, file = embedding_file, sep = "\t", col.names = NA, quote = FALSE)
}
