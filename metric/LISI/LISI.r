#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Kirti Biharie; implemented LISI score

suppressPackageStartupMessages(library(optparse))

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

description <- "Calculate LISI Score"

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
library(lisi)
library(rjson)

if (is.na(opt$ground_truth)) {
  stop("Groundtruth labels needed to calculate the LISI Score")
}

if (is.na(opt$embedding)) {
  stop("Embeddings needed to calculate the LISI Score")
}

if (is.na(opt$config)) {
  stop("Config file not provided")
}

ground_truth <- read.delim(groundtruth_file, sep="\t", row.names=1)
embeddings <- read.delim(embedding_file, sep="\t", row.names=1)
config <- fromJSON(file=config_file)

common_index <- intersect(rownames(ground_truth), rownames(embeddings))
ground_truth <- ground_truth[common_index,,drop=FALSE]
embeddings <- embeddings[common_index,,drop=FALSE]

metric <- mean(compute_lisi(embeddings, ground_truth, "label", perplexity=config$perplexity)[,"label"])

## Write output
outfile <- file(opt$out_file)
dir.create(dirname(opt$out_file), showWarnings = FALSE, recursive = TRUE)

writeLines(format(metric, digits = 6, scientific = TRUE), outfile)
close(outfile)
