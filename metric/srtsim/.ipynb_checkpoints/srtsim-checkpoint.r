#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Lucie Pfeiferova; code for simulation matrix based on tissue/domain SRTsim package

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
description <- "Generate matrix(ces) from package SRTsim"

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
library('SRTsim')

coordinates <- read.delim(coord_file, sep = "\t", row.names = 1)
coordinates <- as.matrix(coordinates[rownames(colData), ])
matrix <- as(Matrix::t(Matrix::readMM(matrix_file)), "CsparseMatrix")
simSRT  <- createSRT(count_in= matrix,loc_in = coordinates)

simSRT1 <- srtsim_fit(simSRT,sim_schem="tissue")
## Generate synthetic data with estimated parameters
simSRT1 <- srtsim_count(simSRT1)

simSRT2 <- srtsim_fit(simSRT,sim_schem="domain")
## Generate synthetic data with estimated parameters
simSRT2 <- srtsim_count(simSRT2)




# from template
dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)

##write matrix
