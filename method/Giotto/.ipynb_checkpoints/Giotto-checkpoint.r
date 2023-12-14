#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Lucie Pfeiferova; functions for Giotto HMRF spatial domain exploring

suppressPackageStartupMessages(library(optparse))

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
## -- make option for betas (eg c(0,8,10) and betas_to_add (eg which one to use)
# TODO adjust description
description <- "Method Giotto HMRF spatial domain"

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


# Seed
seed <- opt$seed
set.seed(seed)
# TODO if the method requires the seed elsewhere please pass it on

# You can use the data as SpatialExperiment
spe <- spe <- get_SpatialExperiment(feature_file = feature_file,observation_file = observation_file,
                                    coord_file = coord_file,matrix_file = matrix_file)

## Your code goes here
# TODO
# The data.frames with observations may contain a column "selected" which you need to use to
# subset and also use to subset coordinates, neighbors, (transformed) count matrix
# label_df = ...  # data.frame with row.names (cell-id/barcode) and 1 column (label)
# embedding_df = NULL  # optional, data.frame with row.names (cell-id/barcode) and n columns
##raw expression counts expected
createGiotto_fn = function(sce,annotation = FALSE, selected_clustering = NULL){
  raw_expr = SummarizedExperiment::assay(sce, "counts")
  #colnames(raw_expr) <- colData(sce)[,"Barcode"]
  #norm_expression <- SummarizedExperiment::assay(sce, "logcounts")
  
  cell_metadata <- SingleCellExperiment::colData(sce)
  cell_metadata$cell_ID <- rownames(SingleCellExperiment::colData(sce))
  colnames(cell_metadata)[c(1,2)] <- c("sdimx","sdimy")
  cell_metadata <- cell_metadata[,4,1,2,3]
  gene_metadata <- as.data.frame(SingleCellExperiment::rowData(sce))
  colnames(gene_metadata)[[1]] <- c("gene_ID")
  #rownames(raw_epxr) <- gene_metadata$gene_ID
  if (annotation) {
    rownames(raw_expr) <- c(SingleCellExperiment::rowData(sce)[,"SYMBOL"])
    #rownames(norm_expression) <- c(SingleCellExperiment::rowData(sce)[,"SYMBOL"])
  }
  gobj = Giotto::createGiottoObject(raw_exprs = raw_expr,
                                    cell_metadata = cell_metadata,
                                    spatial_locs = as.data.frame(SpatialExperiment::spatialCoords(spe)),
                                    gene_metadata = gene_metadata)
  return(gobj)
}
my_giotto_object <- createGiotto_fn(spe)
# my_giotto_object = createGiottoObject(raw_exprs = feature_file,
#                                       spatial_locs = coord_file,norm_expr = matrix_file)

my_giotto_object <- Giotto::normalizeGiotto(my_giotto_object)
# create network (required for binSpect methods)
my_giotto_object = Giotto::createSpatialNetwork(gobject = my_giotto_object,minimum_k = 10)

# identify genes with a spatial coherent expression profile
km_spatialgenes = Giotto::binSpect(my_giotto_object, bin_method = 'rank')

###2. Run HMRF

#if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)

my_spatial_genes = km_spatialgenes[1:100]$genes
HMRF_spatial_genes = Giotto::doHMRF(gobject = my_giotto_object,
                                    betas = c(0, 2, 10),
                            spatial_genes = my_spatial_genes
                            )
#viewHMRFresults2D(gobject = my_giotto_object,
#                    HMRFoutput = HMRF_spatial_genes,
#                   k = 9, betas_to_view = i,
#                  point_size = 2)



## Write output
my_giotto_object = addHMRF(gobject = my_giotto_object,
                           HMRFoutput = HMRF_spatial_genes,
                           k = 10, betas_to_add = c(10),
                           hmrf_name = 'HMRF')



## Write output
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
label_df<- as.data.frame(my_giotto_object@cell_metadata[,c("cell_ID","HMRF_k10_b.10")])
rownames(label_df) <- label_df$cell_ID
label_df <- label_df[,2, drop=FALSE]

colnames(label_df) <- c("label")
write.table(label_df, file = label_file, sep = "\t", col.names = NA, quote = FALSE)


