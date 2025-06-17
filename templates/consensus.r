#!/usr/bin/env Rscript

# Author_and_contribution: Jieran Sun & Mark Robinson; implmented method
# Author_and_contribution: Peiying Cai; created template
# Author_and_contribution: ENTER YOUR NAME AND CONTRIBUTION HERE

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(
    c("-i", "--input_file"),
    type = "character", default = NULL,
    help = "Input containing the aggregated labels."
  ),
  make_option(
    c("--seed"),
    type = "integer", default = NULL,
    help = "seed for input"
  ),
  make_option(
    c("-b", "--base_clusterings"),
    type = "character", default = NULL,
    help = "Path to base-clustering ranking file"
  ), 
  make_option(
    c("--n_clusters"),
    type = "integer", default = NULL,
    help = "Desired number of clusters in the consensus output"
  ),
  make_option(
    c("--n_bcs"),
    type = "integer", default = NULL,
    help = "Desired number of base clustering results feed into the algorithm"
  ),
  make_option(
    c("-o", "--output_file"),
    type = "character", default = NULL,
    help = "desired output file"
  )
)

# TODO adjust description
description <- "Calculate consensus ... for selected BCs"

opt_parser <- OptionParser(
  usage = description,
  option_list = option_list
)
opt <- parse_args(opt_parser)

# Use these filepaths as input
input_file <- opt$input_file
output_file <- opt$output_file
bc_file <- opt$base_clusterings
# TODO adjust numbers in `n_bcs` and `n_clust`
# make sure that `n_clust` exists among the column names of `bc_file`
n_bcs <- ifelse(is.null(opt$n_bcs), 8, opt$n_bcs)
n_clust <- ifelse(is.null(opt$n_clusters), "7", opt$n_clusters)

# Read the selected base clusterings
label_df <- read.delim(input_file, stringsAsFactors = FALSE, row.names = 1, numerals="no.loss")
# bc_list stores all 'method_config_n_clust_label' entries matching n_clust
bc_list <- read.delim(bc_file, stringsAsFactors = FALSE, row.names = 1, numerals="no.loss", check.names=FALSE)[[as.character(n_clust)]]
bc_list <- bc_list[!is.na(bc_list)]
if (length(bc_list) < n_bcs){
  warning(sprintf("Not enough (%s) base clusterings(BCs) are available, use %s BCs instead.", n_bcs, length(bc_list)))
}
bc_list <- bc_list[1:min(n_bcs, length(bc_list))]

# Subset the label data to keep only the selected base clusterings
label_selected <- label_df[, bc_list]

# Make sure all the clusters are ranked 1 to n without jumping (SOTIP)
label_selected <- as.data.frame(lapply(label_selected, function(u){
    unique_labels <- sort(unique(u))
    if (all(unique_labels==seq_along(unique_labels))) {
        return(factor(u, levels = unique_labels))
    } else {
        # Count occurrences of each number
        freq <- table(u)
        rank_map <- rank(-freq, ties.method = "first") # Negative for descending order
        new_vec <- rank_map[as.character(u)]
        new_vec <- factor(as.numeric(new_vec))
        return(new_vec)
    }
}))
rownames(label_selected) <- rownames(label_df)

# Seed
seed <- opt$seed
set.seed(seed)
# TODO if the method requires the seed elsewhere please pass it on

## Your code goes here
# TODO
# label_selected: data frame with spots/cells as rows, base clusterings as columns
# Example:
#                     DRSC_default_10_label scanpy_default_10_label
# AAACAAGTATCTCCCA-1  1                     2
# AAACACCAATAACTGC-1  2                     1
# AAACAGAGCGACTCCT-1  1                     2


## Write output
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

write.table(output_df, file = output_file, sep = "\t", col.names = NA, quote = FALSE)
