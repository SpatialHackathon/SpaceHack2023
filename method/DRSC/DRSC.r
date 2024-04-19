#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Søren Helweg Dam; implemented method

suppressPackageStartupMessages({
    library(optparse)
    library(SingleCellExperiment)
    library(jsonlite)
    library(Seurat)
    library(DR.SC)
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

description <- "DR.SC: Joint dimension reduction and spatial clustering for single-cell/spatial transcriptomics data"

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

# Creating SingleCellExperiment
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

# Opt params
technology <- opt$technology
if(!(technology %in% c("Visium", "ST"))) technology <- "Other_SRT"
n_clusters <- opt$n_clusters
seed <- opt$seed
set.seed(seed)

# Load configuration
feature_method <- config$feature_method
n_genes <- congig$n_genes


# You can use the data as SingleCellExperiment
sce <- get_SingleCellExperiment(
    feature_file = feature_file, 
    observation_file = observation_file,
    coord_file = coord_file, 
    matrix_file = matrix_file)

## Create Seurat and normalize
seurat_obj <- as.Seurat(sce)
seurat_obj <- NormalizeData(seurat_obj, verbose = F)

# Variable features (if given - opt takes priority)
n_genes <- ifelse(is.null(opt$n_genes), n_genes, opt$n_genes)

if (nrow(seurat_obj) > n_genes){
    if (feature_method == "FindVariableFeatures") {
        seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = n_genes, verbose = FALSE)
    } else if (feature_method == "FindSVGs") {
        seurat_obj <- FindSVGs(seurat_obj, nfeatures = n_genes)
    } else if (feature_method == "All genes") {
        seurat_obj[["originalexp"]]@var.features <- rownames(rowData(sce))
    } else {
        stop("Invalid feature selection method in 'desc' configuration.")
    }
} else{
    seurat_obj[["originalexp"]]@var.features <- rownames(rowData(sce))
}



# Fit the DRSC requirement
if (!("row" %in% colnames(seurat_obj@meta.data) & "col" %in% colnames(seurat_obj@meta.data))){
    seurat_obj@meta.data$col <- sce@metadata$spatialCoords[rownames(seurat_obj@meta.data), 1]
    seurat_obj@meta.data$row <- sce@metadata$spatialCoords[rownames(seurat_obj@meta.data), 2]
}
seu_drsc <- DR.SC::DR.SC(seurat_obj, K = n_clusters, platform = technology)

# The data.frames with observations may contain a column "selected" which you need to use to
# subset and also use to subset coordinates, neighbors, (transformed) count matrix
label_df <- data.frame("label" = seu_drsc$spatial.drsc.cluster, row.names=colnames(seu_drsc))
# data.frame with row.names (cell-id/barcode) and 1 column (label)
embedding_df <- as.data.frame(Embeddings(seu_drsc, reduction = "dr-sc"))  # optional, data.frame with row.names (cell-id/barcode) and n columns


## Write output
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

colnames(label_df) <- c("label")
write.table(label_df, file = label_file, sep = "\t", col.names = NA, quote = FALSE)

if (exists("embedding_df")) {
  write.table(embedding_df, file = embedding_file, sep = "\t", col.names = NA, quote = FALSE)
}
