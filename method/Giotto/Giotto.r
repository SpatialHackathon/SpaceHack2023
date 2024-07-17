#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Lucie Pfeiferova; functions for Giotto HMRF spatial domain exploring
# Author_and_contribution: Søren Helweg Dam; created environment setup script, updated environment yaml, added configs, tidied code, HMFR works when adding PCA using runPCA

suppressPackageStartupMessages({
    library(optparse)
    library(jsonlite)
    library(SpatialExperiment)
    library(Giotto)
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
## -- make option for betas (eg c(0,8,10) and betas_to_add (eg which one to use)
description <- "Giotto: Implementation of HMRF from smfishHmrf" # , Leiden, Louvain, randomwalk, SNNclust, kmeans, and hierarchical were removed

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
  #neighbors <- as(Matrix::readMM(neighbors_file), "CsparseMatrix")
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
    #spe <- scuttle::logNormCounts(spe)
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

# SpatialExperiment
spe <- get_SpatialExperiment(
    feature_file = feature_file,
    observation_file = observation_file,
    coord_file = coord_file,
    matrix_file = matrix_file
)

## Configuration
method <- "HMRF"
k <- config$k
n_pcs <- ifelse(is.null(opt$n_pcs), config$n_pcs, opt$n_pcs)
n_genes <- ifelse(is.null(opt$n_genes), config$n_genes, opt$n_genes)
beta <- config$beta
betas <- config$betas
bin_method <- config$bin_method
## Giotto instructions
python_path <- Sys.which(c("python"))
instrs <- createGiottoInstructions(save_dir = out_dir,
                                   save_plot = FALSE,
                                   show_plot = FALSE,
                                   python_path = python_path)



## raw expression counts expected
createGiotto_fn = function(spe, annotation = FALSE, selected_clustering = NULL, instructions = NULL){
  raw_expr <- SummarizedExperiment::assay(spe, "counts")
  #colnames(raw_expr) <- colData(sce)[,"Barcode"]
  #norm_expression <- SummarizedExperiment::assay(spe, "logcounts")
  
  cell_metadata <- SingleCellExperiment::colData(spe)
  cell_metadata$cell_ID <- rownames(SingleCellExperiment::colData(spe))
  colnames(cell_metadata)[c(1,2)] <- c("sdimx", "sdimy")
  cell_metadata <- as.data.frame(cell_metadata[,c(4,1,2,3)])
  feat_metadata <- tibble::as_tibble(SingleCellExperiment::rowData(spe),rownames = "feat_ID")
  #colnames(feat_metadata)[[1]] <- c("feat_ID")
  #rownames(raw_epxr) <- gene_metadata$gene_ID
  if (annotation) {
    rownames(raw_expr) <- c(SingleCellExperiment::rowData(spe)[, "SYMBOL"])
    #rownames(norm_expression) <- c(SingleCellExperiment::rowData(sce)[,"SYMBOL"])
  }
  gobj = Giotto::createGiottoObject(
      expression = list("raw" = raw_expr),
      cell_metadata = cell_metadata,
      spatial_locs = as.data.frame(SpatialExperiment::spatialCoords(spe)),
      feat_metadata = feat_metadata,
      instructions = instructions
  )
  return(gobj)
}
# Convert to Giotto object
gobj <- createGiotto_fn(spe, instructions = instrs)

# Normalize
gobj <- Giotto::normalizeGiotto(gobj, scalefactor = 6000)

## highly variable features (genes)
gobj <- calculateHVF(gobj)

# PCA
gobj <- runPCA(gobject = gobj, center = TRUE, scale_unit = TRUE, name = "PCA", feats_to_use = NULL)
# Using provided PCA throws:
# (Error in cov(y[lclust[[i]], ]) :
#   supply both 'x' and 'y' or a matrix-like 'x')



### Run clustering

#if(!file.exists(hmrf_folder)) dir.create(hmrf_folder, recursive = T)
message("Running ", method, " clustering")

# create spatial network (tricky to use input network)
gobj <- Giotto::createSpatialNetwork(
    gobject = gobj,
    k = k
    #minimum_k = 10
)

# identify genes with a spatial coherent expression profile (Not necessary - default uses all 'selected' features)
km_spatialgenes <- Giotto::binSpect(gobj, bin_method = bin_method)
my_spatial_genes <- ifelse(n_genes == 0, km_spatialgenes$feats, km_spatialgenes[1:n_genes]$feats)


HMRF_spatial_genes <- Giotto::doHMRF(
    gobject = gobj,
    spat_unit = "cell",
    feat_type = "rna",
    betas = betas,
    # expression_values = "normalized", # not used when dim_reduction_to_use is given
    spatial_genes = my_spatial_genes,
    dim_reduction_to_use = "pca",
    dim_reduction_name = "PCA",
    dimensions_to_use = 1:n_pcs,
    k = n_clusters,
    name = method,
    seed = seed
    )

# Manually inspect results:
#viewHMRFresults2D(gobject = gobj,
#                    HMRFoutput = HMRF_spatial_genes,
#                   k = 9, betas_to_view = i,
#                  point_size = 2)

## Add HMRF
gobj <- addHMRF(
    gobject = gobj,
    HMRFoutput = HMRF_spatial_genes,
    k = k, betas_to_add = beta,
    hmrf_name = method)
    

## Write output
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
label_df <- as.data.frame(Giotto::getCellMetadata(gobj, output = "data.table"))
label_df <- data.frame(label = label_df[length(colnames(label_df))], row.names = label_df[[1]])
colnames(label_df) <- "label"

#print(table(label_df$label))

write.table(label_df, file = label_file, sep = "\t", col.names = NA, quote = FALSE)

# Clean HMRF_output folder
system("rm -rf HMRF_output")
