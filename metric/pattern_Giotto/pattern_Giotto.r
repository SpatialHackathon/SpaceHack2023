#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Lucie Pfeiferova; simulation dataset from Giotto:runPatternSimulation, matrix needing

suppressPackageStartupMessages(library(optparse))

# TODO adjust description
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

# TODO adjust description
description <- "Create simulated dataset based on spatial pattern and gene selection"

opt_parser <- OptionParser(
  usage = description,
  option_list = option_list
)
opt <- parse_args(opt_parser)

# Use these filepaths as input
label_file <- opt$labels

if (!is.na(opt$ground_truth)) {
  groundtruth_file <- opt$ground_truth
}
if (!is.na(opt$embedding)) {
  embedding_file <- opt$embedding
}
if (!is.na(opt$config)) {
  config_file <- opt$config
}
outfile <- file(opt$out_file)

## Your code goes here
# metric = ... #  output of the metric as float
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
cell_ids <- as.data.frame(my_giotto_object@spatial_locs)
pattern_cell_ids <- cell_ids[as.numeric(cell_ids$sdimx) < mean(as.numeric(cell_ids$sdimx))/4*3,]
random_genes <- sample(km_spatialgenes$genes, 100)
pattern_sim <- Giotto::runPatternSimulation(my_giotto_object,pattern_colors = c(`in` = "green", out = "red"),reps = 2,
                                            pattern_cell_ids = pattern_cell_ids$cell_ID,
                                            gene_names = random_genes,
                                            spat_methods_params = list(bin_method = 'rank'),
                                            spat_methods_names = c("binSpect_single" ),
                                            spat_methods = c('binSpect_single'),save_raw = T,save_plot = F,
                                            save_dir=dirname(out_file))

## Write output
## output is new matrix with generated data (raw counts), output in out_file directory

