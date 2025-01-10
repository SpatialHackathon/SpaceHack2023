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
    c("--seed"),
    type = "integer", default = NULL,
    help = "seed for input"
  ),
  make_option(
    c("-c", "--columns"),
    type = "character", default = NULL,
    help = "List of column names for the input BCs"
  ), 
  make_option(
    c("-o", "--output_file"),
    type = "character", default = NULL,
    help = "desired output file"
  )
)

description <- "Calculate consensus for selected BCs"

opt_parser <- OptionParser(
  usage = description,
  option_list = option_list
)
opt <- parse_args(opt_parser)

# Use these filepaths as input
input_file <- opt$input_file
output_file <- opt$output_file
columns <- opt$columns
seed <- opt$seed

# Your code goes here
suppressPackageStartupMessages({
    library(diceR)
})

label_df <- read.delim(input_file, stringsAsFactors = FALSE, row.names = 1, numerals="no.loss")
label_selected <- label_df[, columns]

# Make sure all the clusters are ranked 1 to n without jumping (SOTIP)
label_selected <- apply(label_selected, 2, function(u){
    unique_labels <- sort(unique(u))
    if (all(unique_labels==seq_along(unique_labels))) {
        return(as.factor(u))
    } else {
        # Count occurrences of each number
        freq <- table(u)
        rank_map <- rank(-freq, ties.method = "first") # Negative for descending order
        new_vec <- rank_map[as.character(vec)]
        new_vec <- as.factor(as.numeric(new_vec))
        return(new_u)
    }
})

kmode_vec <- diceR:::k_modes(label_selected, is.relabelled = FALSE, seed = seed)
kmode_df <- data.frame(consensus_kmode=kmode_vec, row.names = row.names(label_selected))

# Save the results
write.table(kmode_df, file = output_file, sep = "\t", col.names = NA, quote = FALSE)
