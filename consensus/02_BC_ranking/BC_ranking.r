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
    c("-o", "--output_file"),
    type = "character", default = NULL,
    help = "desired output file"
  ),
  make_option(
    c("--selection_metrics"),
    type = "character", default = NULL,
    help = "file containing the metric information for BC selection"
  ),
  make_option(
    c("-mp", "--max_percentage"),
    type = "double", default = NULL,
    help = "maximal percentage of the largest class"
  )
)

description <- "Automatically select the base-clusterings based on different algorithms"

opt_parser <- OptionParser(
  usage = description,
  option_list = option_list
)
opt <- parse_args(opt_parser)

# Use these filepaths as input
input_file <- opt$input_file
output_file <- opt$output_file
smoothness_file <- opt$smoothness
ari_file <- opt$ari
max_percentage <- ifelse(is.null(opt$max_percentage), 0.8, opt$max_percentage)

##### Load files
label_df <- read.delim(input_file, stringsAsFactors = FALSE, row.names = 1, numerals="no.loss")

##### Filter out class-imbalanced case
label_df[, sapply(label_df, function(col) {
  max(table(col)) / length(col) <= max_percentage
})]

##### Separate label_df into different cluster number ones
n_clusters <- apply(label_df, 2, function(u){length(unique(u))})
label_lists <- split(names(n_clusters), n_clusters)

##### Select base-clsuterings based on algorithms

if (!is.null(ari_file)){
    # ARI df is a nxn dataframe with n = number of results instance, each value is the cross-ARI result
    selection_df <- read.delim(smoothness_file, stringsAsFactors = FALSE, row.names = 1, numerals="no.loss")
    s_mean <- colMeans(selection_df)
} else {
    if (is.null(smoothness_file)){
      file_searched <- list.files(path = dirname(input_file), pattern = "smoothness", full.names = TRUE)
      if (length(file_searched) == 0){
        warning("No smoothness entropy file found or defined.")
      }
      smoothness_file <- file_searched[1]
    }
    # Smoothness df is a one-column dataframe with row names refers to individual result instance
    selection_df <- read.delim(smoothness_file, stringsAsFactors = FALSE, row.names = 1, numerals="no.loss")
    # Reversed it as high entropy means less smooth
    s_mean <- - rowMeans(selection_df)
}

d_length <- max(lengths(label_lists))
s_bc_list <- sapply(label_lists, function(nclu_names){
    s_n <- s_mean[names(s_mean) %in% nclu_names]
    selected_names <- row.names(s_n)[order(s_n, decreasing = TRUE)]
    length(selected_names) <- d_length
    return(selected_names)
})

result_df <- as.data.frame(s_bc_list)

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

# Save the results
write.table(result_df, file = output_file, sep = "\t", col.names = NA, quote = FALSE)
