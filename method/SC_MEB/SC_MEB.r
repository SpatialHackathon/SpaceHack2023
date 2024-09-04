#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Søren Helweg Dam; implemented method

suppressPackageStartupMessages({
    library(optparse)
    library(jsonlite)
    library(SingleCellExperiment)
    library(scuttle)
    library(scran)
    library(scater)
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
}
if (!is.na(opt$image)) {
  image_file <- opt$image
}
if (!is.na(opt$config)) {
  config_file <- opt$config
  config <- fromJSON(config_file)
}

technology <- opt$technology
if (!technology %in% c("Visium", "ST")) {
    warning("This method was built for Visium and ST. 
            Setting technology to ST.")
    technology <- "ST"
    # ST: Square spots
    # Visium: Hexagonal spots
}
n_clusters <- opt$n_clusters
seed <- opt$seed
# Load configuration
n_pcs <- ifelse(is.null(opt$n_pcs), config$n_pcs, opt$n_pcs)
n_genes <- ifelse(is.null(opt$n_genes), config$n_genes, opt$n_genes)

# SingleCellExperiment
get_SingleCellExperiment <- function(
    feature_file,
    observation_file,
    coord_file,
    matrix_file = NA,
    assay_name = "counts",
    reducedDim_name = "reducedDim") {
  rowData <- read.delim(feature_file, stringsAsFactors = FALSE, row.names = 1, numerals="no.loss")
  colData <- read.delim(observation_file, stringsAsFactors = FALSE, row.names = 1, numerals="no.loss")

  coordinates <- read.delim(coord_file, sep = "\t", row.names = 1, numerals="no.loss")
  coordinates <- as.matrix(coordinates[rownames(colData), ])
  mode(coordinates) = "numeric"

  sce <- SingleCellExperiment::SingleCellExperiment(
  rowData = rowData, colData = colData, metadata = list("spatialCoords" = coordinates))

  if (!is.na(matrix_file)) {
    assay(sce, assay_name, withDimnames = FALSE) <- as(Matrix::t(Matrix::readMM(matrix_file)), "CsparseMatrix")
    assay(sce, "logcounts", withDimnames = FALSE) <- as(Matrix::t(Matrix::readMM(matrix_file)), "CsparseMatrix")
  }

  # Filter features and samples
  if ("selected" %in% colnames(rowData(sce)) && FALSE) { # Don't subset here
    sce <- sce[as.logical(rowData(sce)$selected), ]
  }
  if ("selected" %in% colnames(colData(sce))) {
    sce <- sce[, as.logical(colData(sce)$selected)]
  }

  return(sce)
}

sce <- get_SingleCellExperiment(
    feature_file = feature_file, 
    observation_file = observation_file,
    coord_file = coord_file, 
    matrix_file = matrix_file)

# Preprocess
set.seed(seed)
sce <- spatialPreprocess(
    sce,
    platform = technology,
    n.PCs = n_pcs,
    n.HVGs = n_genes)

# Fidn neighbors
Adj_sp  <- find_neighbors2(sce, platform = technology)

y <- reducedDim(sce, "PCA")[,seq_len(n_pcs)]

## Your code goes here
fit <- SC.MEB(y, Adj_sp, K_set = n_clusters, num_core = 2)
label_df <- data.frame("label" = unlist(fit["x", ]), row.names = colnames(sce))


## Write output
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

colnames(label_df) <- c("label")
write.table(label_df, file = label_file, sep = "\t", col.names = NA, quote = FALSE)

if (exists("embedding_df")) {
  write.table(embedding_df, file = embedding_file, sep = "\t", col.names = NA, quote = FALSE)
}
