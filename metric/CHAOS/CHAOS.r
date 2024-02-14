#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Kirti Biharie; implemented CHAOS score

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
    c("-x", "--coordinates"),
    type = "character", default = NA,
    help = "Coordinates of points in the spatial domain."
  ),
  make_option(
    c("-o", "--out_file"),
    type = "character", default = NULL,
    help = "Output file."
  )
)

description <- "Calculate CHAOS Score"

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
if (!is.na(opt$coordinates)) {
  coordinate_file <- opt$coordinates
}

## Your code goes here
library(pdist)

if (is.na(opt$coordinates)) {
  stop("Coordinates needed to calculate the CHAOS Score")
}

labels <- read.delim(label_file, sep="\t", row.names=1)
coordinates <- read.delim(coordinate_file, sep="\t", row.names=1)

common_index <- intersect(rownames(labels), rownames(coordinates))
labels <- labels[common_index,,drop=FALSE]
coordinates <- coordinates[common_index,,drop=FALSE]

# Copied from SpatialPCA (Lulu Shang, and Xiang Zhou (2022))

#' @title Calculate CHAOS score to measure clustering performance.
#' @description CHAOS score measures the spatial continuity of the detected spatial domains.
#' Lower CHAOS score indicates better spatial domian clustering performance.
#' @param clusterlabel Cluster labels.
#' @param location A n by k matrix of spatial locations.
#' @return A numeric value for CHAOS score.
fx_CHAOS = function(clusterlabel, location){
  matched_location=location
  NAs = which(is.na(clusterlabel))
  if(length(NAs>0)){
    clusterlabel=clusterlabel[-NAs]
    matched_location = matched_location[-NAs,]
  }
  matched_location = scale(matched_location)
  dist_val = rep(0,length(unique(clusterlabel)))
  count = 0
  for(k in unique(clusterlabel)){
    count = count + 1
    location_cluster = matched_location[which(clusterlabel == k),]
    if(length(location_cluster)==2){next}
    #require(parallel)
    results = lapply(1:dim(location_cluster)[1], fx_1NN, location_in=location_cluster)
    dist_val[count] = sum(unlist(results))
  }
  dist_val = na.omit(dist_val)
  return(sum(dist_val)/length(clusterlabel))
}

fx_1NN = function(i,location_in){
  line_i = rep(0,dim(location_in)[1])
  line_i = pdist(location_in[i,],location_in[-i,])@dist
  return(min(line_i))
}

metric <- fx_CHAOS(labels$label, coordinates)

## Write output
outfile <- file(opt$out_file)
dir.create(dirname(opt$out_file), showWarnings = FALSE, recursive = TRUE)

writeLines(format(metric, digits = 6, scientific = TRUE), outfile)
close(outfile)
