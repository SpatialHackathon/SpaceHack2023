#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Mark D. Robinson; coded the domain-specific F1

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
description <- "Calculate domain-specific F1 score (returns JSON with vector: F1 for each true domain)"

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


## Code for calculating metric goes here
## --------------------------------------
library(cluster)

# # for testing - start
label_file <- "~/scratch/SpaceHack2/method_results/LIBD_DLPFC/Br5292_151507/SpaGCN/domains.tsv"
outfile <- "./cluster-specific-silhouette.json"
embedding_file <- "~/scratch/SpaceHack2/method_results/LIBD_DLPFC/Br5292_151507/log1p/hvg/pca_20/dimensionality_reduction.tsv"
# groundtruth_file <- "data/libd_dlpfc/Br5595_151670/labels.tsv"
# # for testing - stop

domains <- read.delim(label_file, sep="\t", row.names = 1)
embedding <- read.delim(embedding_file, sep="\t", row.names = 1)

rn <- intersect(rownames(domains), rownames(embedding))

# subset to common set
embedding <- embedding[rn,,drop = FALSE]
domains <- domains[rn,,drop = FALSE]

# calculate silhouette score on the embeddings per cell
sil <- silhouette(list(clustering=domains$label),
                  dist=dist(embedding))

agg_median <- aggregate(sil[,"sil_width",drop=FALSE], 
                        list(cluster=sil[,"cluster"]), FUN = median)
agg_mean <- aggregate(sil[,"sil_width",drop=FALSE], 
                      list(cluster=sil[,"cluster"]), FUN = mean)

(df <- data.frame(cluster = agg_mean$cluster,
                  mean_sil_width = agg_mean$sil_width,
                  median_sil_width = agg_median$sil_width))

## Write output
outfile <- file(opt$out_file)
write_json(df,outfile)

# from template
# dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
# 
# writeLines(format(metric, digits = 6, scientific = TRUE), outfile)
# close(outfile)
