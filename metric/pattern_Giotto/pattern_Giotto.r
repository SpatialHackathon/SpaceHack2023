#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: ENTER YOUR NAME AND CONTRIBUTION HERE

suppressPackageStartupMessages(library(optparse))

# TODO adjust description
option_list <- list(
  make_option(
    c("-l", "--labels"),
    type = "character", default = NULL,
    help = "Labels from domain clustering."
  ),
  make_option(
    c("-g", "--ground_truth"),
    type = "character", default = NA,
    help = "Groundtruth labels."
  ),
  make_option(
    c("-e", "--embedding"),
    type = "character", default = NA,
    help = "Embedding of points in latent space. Potential usage for metrics without groundtruth."
  ),
  # format should be json
  make_option(
    c("-c", "--config"),
    type = "character", default = NA,
    help = "Optional config file (json) used to pass additional parameters."
  ),
  make_option(
    c("-o", "--out_file"),
    type = "character", default = NULL,
    help = "Output file."
  )
)

# TODO adjust description
description <- "Calculate metric ..."

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
km_spatialgenes = Giotto::binSpect(my_giotto_object, bin_method = 'rank')
my_spatial_genes = km_spatialgenes[1:500]$genes

pattern_cell_ids = my_giotto_object@cell_metadata$cell_ID[1:500]
pattern_sim <- Giotto::runPatternSimulation(my_giotto_object,pattern_colors = c(`in` = "green", out = "red"),reps = 2,
                                            pattern_cell_ids = my_giotto_object@cell_metadata$cell_ID[1:500],
                                            gene_names = my_spatial_genes,
                                            spat_methods_params = list(bin_method = 'rank'),
                                            spat_methods_names = c("binSpect_single" ),
                                            spat_methods = c('binSpect_single'))

## Write output
##outfile <- file(opt$out_file)
##dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)

