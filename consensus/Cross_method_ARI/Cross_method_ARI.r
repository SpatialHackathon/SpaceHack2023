#!/usr/bin/env Rscript

# Author_and_contribution: Jieran Sun & Mark Robinson; Create the script

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
)

description <- "Return cross-method ARI"

opt_parser <- OptionParser(
  usage = description,
  option_list = option_list
)
opt <- parse_args(opt_parser)

# Use these filepaths as input
input_file <- opt$input_file
output_file <- opt$output_file

set.seed(seed)

suppressPackageStartupMessages({
  library(mclust)
})

##### Define ARI function
calc_aris <- function(m, flavour="ARI") {
  a <- diag(ncol(m))
  for(i in 1:(ncol(m)-1))
    for(j in 2:ncol(m)) {
      if(flavour=="ARI") {
        require(mclust)
        a[i,j] <- a[j,i] <- mclust::adjustedRandIndex(m[,i], m[,j])
      } else if(flavour=="sARI") {
        # TODO implement spatial ARI
      }
    }
  rownames(a) <- colnames(a) <- colnames(m)
  a
}

##### Load files
label_df <- read.delim(input_file, stringsAsFactors = FALSE, row.names = 1, numerals="no.loss")

# Calculate cross-method ARI
ari_mat <- calc_aris(label_df)

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)

# Save the results
write.table(ari_mat, file = output_file, sep = "\t", col.names = NA, quote = FALSE)
