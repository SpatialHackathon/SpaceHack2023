#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Kim Vucinic; modified template and created script

suppressPackageStartupMessages(library(optparse))

# Arguments
option_list <- list(
  make_option(
    c("-l", "--labels"),
    type = "character", default = NULL,
    help = "Labels from domain clustering. Path to labels (as tsv)."
  ),
  make_option(
    c("--seed"),
    type = "integer", default = NULL,
    help = "Seed to use for random operations."
  ),
  make_option(
    c("-o", "--out_file"),
    type = "character", default = NULL,
    help = "Output file."
  )
)

# Description
description <- "Shuffling labels..."

opt_parser <- OptionParser(
  usage = description,
  option_list = option_list
)
opt <- parse_args(opt_parser)

# Use these filepaths as input
label_file <- opt$labels
seed <- opt$seed

## Your code goes here



## Write output
outfile <- file(opt$out_file)
dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)

write.table(df_shuffled, outfile, sep = "\t", col.names = NA, quote = "FALSE")