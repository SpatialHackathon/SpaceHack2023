#!/usr/bin/env Rscript

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Mark D. Robinson; coded the domain-specific F1

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
description <- "Calculate domain-specific F1 score (returns JSON with vector: F1 for each true domain)"

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
library(clue)
library(jsonlite)

# # for testing - start
# label_file <- "results/libd_dlpfc/Br5595_151670/SpaGCN/domains.tsv"
# groundtruth_file <- "data/libd_dlpfc/Br5595_151670/labels.tsv"
# outfile <- "./domain-specific-f1.json"
# # for testing - stop

domains <- read.delim(label_file, sep="\t", row.names = 1)
groundtruth <- read.delim(groundtruth_file, sep="\t", row.names = 1)

rn <- intersect(rownames(domains), rownames(groundtruth))

# subset to common set
domains <- domains[rn,,drop = FALSE]
groundtruth <- groundtruth[rn,,drop = FALSE]

n_gt_labels <- length(unique(groundtruth$label))
n_clust_labels <- length(unique(domains$label))


tb <- table(domains$label, groundtruth$label)

# hungarian solver will fail, pad with 0 columns
if( n_clust_labels > n_gt_labels ) {
  n_cols <- n_clust_labels - n_gt_labels
  for(i in seq_len(n_cols)) {
    tb <- cbind(tb,0)
    colnames(tb)[ncol(tb)] <- paste0("dummy",i)
  }
}

print(tb)

# for the case that there are more clusters than true labels, need to 
# add dummy variables

hg <- clue::solve_LSAP(tb, maximum = TRUE)
sa <- seq_along(hg)
mapping <- data.frame(cluster = rownames(tb)[sa], 
                      label = colnames(tb)[hg],
                      tp = tb[cbind(seq_along(hg), hg)])
rownames(mapping) <- mapping$label

tt <- table(groundtruth$label)
mapping[names(tt),"domain_size"] <- tt

tl <- table(domains$label)
m <- match(mapping$cluster, names(tl))
mapping$cluster_size <- as.numeric(tl[m])

mapping$recall <- with(mapping, tp/domain_size)
mapping$precision <- with(mapping, tp/cluster_size)
mapping$f1 <- with(mapping, 2/(1/recall+1/precision))

# remove dummy labels
(mapping <- mapping[!grepl("dummy", mapping$label),])

## Write output
outfile <- file(opt$out_file)
write_json(unname(mapping[,c("label","f1")]),outfile)

# from template
# dir.create(dirname(outfile), showWarnings = FALSE, recursive = TRUE)
# 
# writeLines(format(metric, digits = 6, scientific = TRUE), outfile)
# close(outfile)
