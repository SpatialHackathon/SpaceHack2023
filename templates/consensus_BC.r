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
    c("-o", "--output_file"),
    type = "character", default = NULL,
    help = "desired output file"
  )
# TODO add additional arguments
)

# TODO adjust description
description <- "... to select base clusterings"

opt_parser <- OptionParser(
  usage = description,
  option_list = option_list
)
opt <- parse_args(opt_parser)

# Use these filepaths as input
input_file <- opt$input_file
output_file <- opt$output_file

##### Load files
label_df <- read.delim(input_file, stringsAsFactors = FALSE, row.names = 1, numerals="no.loss")

## Your code goes here
# TODO
# output_df: data frame with number of clusters as column headers, and clustering label names as values
# Example:
# 7                          8
# method1_default_7_label    method1_default_8_label
# method2_default_7_label    method2_default_8_label
# method3_default_7_label    method3_default_8_label


## Write output
dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

# Save the results
write.table(output_df, file = output_file, sep = "\t", col.names = NA, quote = FALSE)
