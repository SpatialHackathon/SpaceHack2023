#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Kim Vucinic; modified template and created script

suppressPackageStartupMessages(library(optparse))

# Arguments
option_list <- list(
  make_option(
    c("-c", "--coordinates"),
    type = "character", default = NULL,
    help = "Path to coordinates (as tsv)."
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
description <- "Shuffling coordinates in coordinates.tsv"

opt_parser <- OptionParser(
  usage = description,
  option_list = option_list
)
opt <- parse_args(opt_parser)

# Use these filepaths as input
coord_file <- opt$coordinates

# Seed
seed <- opt$seed
set.seed(seed)

## Your code goes here
df <- read.delim(coord_file, sep = "\t", row.names = 1)
if (any(!(c("x", "y") %in% colnames(df)))){
     stop("X and y coordinates are not present in the file. Check your file.")
}

# Randomize IDs, but keep the same order of IDs (not really necessary)
df_order <- rownames(df)
rownames(df) <- sample(rownames(df))
df_final <- df[order(match(rownames(df), df_order)),]

## Write output
outfile <- file(opt$out_file)
dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)

write.table(df_final, outfile, sep = "\t", col.names = NA, quote = "FALSE")