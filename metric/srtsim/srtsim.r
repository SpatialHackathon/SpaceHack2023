#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Lucie Pfeiferova; code for simulation matrix based on 2 methods from SRTsim package, ref data needing

suppressPackageStartupMessages(library(optparse))

# TODO adjust description
option_list <- list(
  make_option(
    c("-c", "--coordinates"),
    type = "character", default = NULL,
    help = "Path to coordinates (as tsv)."
  ),
  make_option(
    c("-m", "--matrix"),
    type = "character", default = NA,
    help = "Path to (transformed) counts (as mtx)."
  ),
  make_option(
    c("-f", "--features"),
    type = "character", default = NULL,
    help = "Path to features (as tsv)."
  ),
  make_option(
    c("-o", "--observations"),
    type = "character", default = NULL,
    help = "Path to observations (as tsv)."
  ),
  make_option(
    c("-n", "--neighbors"),
    type = "character", default = NA,
    help = "Path to neighbor definitions. Square matrix (not necessarily symmetric) where each row contains the neighbors of this observation (as mtx)."
  ),
  make_option(
    c("-d", "--out_dir"),
    type = "character", default = NULL,
    help = "Output directory."
  ),
  make_option(
    c("--dim_red"),
    type = "character", default = NA,
    help = "Reduced dimensionality representation (e.g. PCA)."
  ),
  make_option(
    c("--image"),
    type = "character", default = NA,
    help = "Path to H&E staining."
  ),
  make_option(
    c("--n_clusters"),
    type = "integer", default = NULL,
    help = "Number of clusters to return."
  ),
  make_option(
    c("--technology"),
    type = "character", default = NULL,
    help = "The technology of the dataset (Visium, ST, imaging-based)."
  ),
  make_option(
    c("--seed"),
    type = "integer", default = NULL,
    help = "Seed to use for random operations."
  ),
  make_option(
    c("--config"),
    type = "character", default = NA,
    help = "Optional config file (json) used to pass additional parameters."
  )
)

# TODO adjust description
description <- "Generate matrix(ces) from package SRTsim - matrix and coordinates needed"

opt_parser <- OptionParser(
  usage = description,
  option_list = option_list
)
opt <- parse_args(opt_parser)

# Use these filepaths as input
label_file <- opt$labels



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

##write matrix(ces)
