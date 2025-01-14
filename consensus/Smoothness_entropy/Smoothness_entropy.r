#!/usr/bin/env Rscript

# Author_and_contribution: Jieran Sun & Mark Robinson; Create the script

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(
    c("-i", "--input_file"),
    type = "character", default = NULL,
    help = "Input containing the aggregated labels."
  ),
  make_option(
    c("-c", "--coordinates"),
    type = "character", default = NULL,
    help = "file path to the spatial coordinates of the spots/cells"
  ),
  make_option(
    c("-o", "--output_file"),
    type = "character", default = NULL,
    help = "desired output file"
  ),
  make_option(
    c("-n", "--neighbors"),
    type = "integer", default = NULL,
    help = "Number of neighbors to calculate the smoothness"
  ),
  make_option(
    c("-s", "--seed"),
    type = "integer", default = NULL,
    help = "seed for neighboring algorithm"
  )
)

description <- "Calculate overall smoothness of the clustering"

opt_parser <- OptionParser(
  usage = description,
  option_list = option_list
)
opt <- parse_args(opt_parser)

# Use these filepaths as input
input_file <- opt$input_file
output_file <- opt$output_file
coord_file <- opt$coordinates
neighbors <- ifelse(is.null(opt$neighbors), 6, opt$neighbors)
seed <- ifelse(is.null(opt$seed), 2025, opt$seed)



set.seed(seed)

suppressPackageStartupMessages({
  library(dbscan)
})

##### Define function
calc_entropy <- function(u) {
  p <- u[u>0]
  p <- p/sum(p)
  -sum(p*log(p))
}

spot_entropy <- function(spatial_coords, label, k) {
  suppressPackageStartupMessages(require(dbscan))
  knns <- dbscan::kNN(spatial_coords, k=k)
  label <- as.factor(label)
  neighb_labels <- apply(knns$id, 2, function(u) label[u])
  apply(neighb_labels, 
        1, function(u) calc_entropy(table(factor(u, levels=levels(label)))))
}

##### Load files
label_df <- read.delim(input_file, stringsAsFactors = FALSE, row.names = 1, numerals="no.loss")
coord_df <- read.delim(coord_file, stringsAsFactors = FALSE, row.names = 1, numerals="no.loss")
coord_df <- coord_df[row.names(label_df), ]

# Calculate point-wise entropy with its neighbor
spot_entropy_df <- apply(label_df, 2,
                         function(u) spot_entropy(coord_df, u, neighbors))

# calculate colmeans and save it to a dataframe
sm_df <- data.frame(smoothness = colMeans(spot_entropy_df))

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

# Save the results
write.table(sm_df, file = output_file, sep = "\t", col.names = NA, quote = FALSE)
