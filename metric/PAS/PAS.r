#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Kirti Biharie; implemented PAS score

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

description <- "Calculate Percentage of Abnormal Spots(PAS)"

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
  stop("Coordinates needed to calculate the PAS Score")
}

labels <- read.delim(label_file, sep="\t", row.names=1)
coordinates <- read.delim(coordinate_file, sep="\t", row.names=1)

common_index <- intersect(rownames(labels), rownames(coordinates))
labels <- labels[common_index,,drop=FALSE]
coordinates <- coordinates[common_index,,drop=FALSE]

# Copied from SpatialPCA (Lulu Shang, and Xiang Zhou (2022))

#' @title Calculate PAS score to measure clustering performance.
#' @description PAS score measures the randomness of the spots that located outside of the spatial region where it was clustered to.
#' Lower PAS score indicates better spatial domian clustering performance.
#' @param clusterlabel Cluster labels.
#' @param location A n by k matrix of spatial locations.
#' @return A numeric value for PAS score.
fx_PAS = function(clusterlabel, location){
  matched_location=location
  NAs = which(is.na(clusterlabel))
  if(length(NAs>0)){
    clusterlabel=clusterlabel[-NAs]
    matched_location = matched_location[-NAs,]
  }

  results = lapply(1:dim(matched_location)[1], fx_kNN, location_in=matched_location,k=10,cluster_in=clusterlabel)
  return(sum(unlist(results))/length(clusterlabel))
}

fx_kNN = function(i,location_in,k,cluster_in){
  line_i = rep(0,dim(location_in)[1])
  line_i = pdist(location_in[i,],location_in[-i,])@dist
  ind = order(line_i)[1:k]
  cluster_use = cluster_in[-i]
  if(sum(cluster_use[ind] != cluster_in[i])>(k/2)){
    return(1)
  }else{
    return(0)
  }
}

metric <- fx_PAS(labels$label, coordinates)

## Write output
outfile <- file(opt$out_file)
dir.create(dirname(opt$out_file), showWarnings = FALSE, recursive = TRUE)

writeLines(format(metric, digits = 6, scientific = TRUE), outfile)
close(outfile)
