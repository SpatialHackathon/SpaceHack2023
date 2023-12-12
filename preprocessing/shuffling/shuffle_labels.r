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

# Seed
seed <- opt$seed
set.seed(seed)

## Your code goes here
df <- read.delim(label_file, sep = "\t", row.names = 1)
if (!("label" %in% colnames(df))){
     stop("Label column not present in the file. Check your file.")
}

# Randomize labels
df_randomized <- data.frame(label = sample(df$label))
rownames(df_randomized) <- rownames(df)

## Write output
outfile <- file(opt$out_file)
write.table(df_randomized, outfile, sep = "\t", col.names = NA, quote = FALSE)