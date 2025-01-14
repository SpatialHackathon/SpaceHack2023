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
    help = "desired output file (spot-wise cross-method entropy)"
  ),
  make_option(
    c("--BC_ranking"),
    type = "character", default = NULL,
    help = "Optional. Define which BCs are used as input for CME calculation. Can either be path to a tsv file or a string of BC instance separated by ','"
  ),
  make_option(
    c("-n", "--n_clusters"),
    type = "integer", default = NULL,
    help = "Desired cluster number. If not NULL, the input file will be filtered based on the number of clusters"
  ),
  make_option(
    c("--n_bcs"),
    type = "integer", default = 8,
    help = "Number of base_clusterings to select as input. Default is 8"
  ),
  make_option(
    c("--BC_output"),
    type = "character", default = NULL,
    help = "If defined, it will output a this file showing which BC instances are used for calculating the CME."
  )
)

description <- "Calculate cross-method entropy at each point. All input columns must have the same number of clusters."

opt_parser <- OptionParser(
  usage = description,
  option_list = option_list
)
opt <- parse_args(opt_parser)

# Use these filepaths as input
input_file <- opt$input_file
output_file <- opt$output_file
bc_ranking <- opt$BC_ranking
n_clust <- opt$n_clusters
n_bcs <- opt$n_bcs

suppressPackageStartupMessages({
  library(dbscan)
  library(clue)
})

##### Define function
calc_entropy <- function(u) {
  p <- u[u>0]
  p <- p/sum(p)
  -sum(p*log(p))
}

align_classes <- function(d, ref, columns=NULL) {
  suppressPackageStartupMessages(require(clue))
  bcs <- as.data.frame(d)

  # Adding a columns filtering here
  if (!is.null(columns)){
    all_cols <- colnames(bcs)
    bcs <- bcs[,union(columns, ref)]
  }

  for(j in 1:ncol(bcs)){
    bcs[,j] <- factor(bcs[,j], levels = as.factor(as.numeric(unique(bcs[,j]))))
  }
  levs <- lapply(bcs, levels)
  n_levs <- sapply(levs, length)
  
  cols_to_change <- setdiff(colnames(bcs), ref)
  
  for(i in cols_to_change) {
    if(n_levs[i] != n_levs[ref]) {
      message(paste0(i," has ",n_levs[i],
                     " levels; reference (", ref, 
                     ") has ",n_levs[ref], " (not modifying)"))
      next
    }
    hung <- clue::solve_LSAP(table(bcs[,ref], bcs[,i]), 
                             maximum = TRUE)
    lookup <- cbind(seq_along(hung), levs[[i]][hung])
    levels(bcs[,i]) <- lookup[order(as.integer(lookup[,2])),1]
    bcs[,i] <- as.factor(as.character(bcs[,i]))
  }

  # Adding back the remaining columns 
  if (!is.null(columns)){
    bcs <- cbind(bcs, as.data.frame(d)[, setdiff(all_cols, colnames(bcs))])
  }
  
  return(bcs)
}

cross_meth_entropy <- function(df){
  cm_ent_list <- lapply(colnames(df), function(c){
    # align the classes with the reference
    df <- align_classes(d=df, ref=c)
    # Extract existing classes
    lv <- levels(df[[c]])
    # For each spot(row), return the frequency of each class
    ps <- apply(df, 1,
                function(u) table(factor(u, levels = lv))) 
    # since apply return the results by cbind, transpose it such that rows refers to the spots and columns means the cluster classes
    ps <- t(ps)
    # Calculate the entropy spot-wise
    ent <- apply(ps, 1, calc_entropy)
    # attach it with same name
    names(ent) <- row.names(df)
  })
  # Take the mean of the entropy for each spot
  mean_ent <- rowMeans(as.data.frame(cm_ent_list))
  return(mean_ent)
}

##### Load files
label_df <- read.delim(input_file, stringsAsFactors = FALSE, row.names = 1, numerals="no.loss")

# Calculate point-wise entropy with its neighbor
if (!is.null(bc_ranking)){
  if (file.exists(bc_ranking)){
    bc_df <- read.delim(bc_ranking, stringsAsFactors = FALSE, row.names = 1, numerals="no.loss")

    if (!is.null(n_clust) && as.character(n_clust) %in% colnames(bc_df)){
      columns <- bc_df[[as.character(n_clust)]]
    } else {
      waning(sprintf("%s is not in the column names of BC_ranking fiel, use the first column instead.", as.character(n_clust)))
      columns <- bc_df[,1]
    }

    columns <- columns[!is.na(columns)]
    columns <- columns[1:min(n_bcs, length(columns))]
  } else {
    columns <- unlist(strsplit(opt$bc_ranking, ",", fixed=TRUE))
  }

  result_ent <- cross_meth_entropy(label_df[, columns])
} else {
  # Split name based on cluster numbers
  n_clusters <- apply(label_df, 2, function(u){length(unique(u))})
  label_lists <- split(names(n_clusters), n_clusters)

  if (!is.null(n_clust) && as.character(n_clust) %in% names(label_lists)){
    cn <- as.character(n_clust)
    columns <-  label_lists[[cn]]
    columns <- columns[1:min(n_bcs, length(columns))]
    result_ent <- cross_meth_entropy(label_df[, columns])
  } else {
    columns <- label_lists
    result_ent <- sapply(columns, function(cn_names){
      cn_names <- cn_names[1:min(n_bcs, length(cn_names))]
      ent <- cross_meth_entropy(label_df[, cn_names])
      return(ent)
    })
  }
}

# Save the results
result_df <- as.data.frame(result_list)
if (ncol(result_df)==1){colnames(result_df)<-c("entropy")}

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
write.table(result_df, file = output_file, sep = "\t", col.names = NA, quote = FALSE)

# add meta-data file
if (!is.null(opt$BC_output)){

  if (class(columns)=="list"){
    max_length <- max(lengths(columns))

    column_df <- do.call(cbind, lapply(columns, function(c){
      length(c) <- max_length
      return(c)
    }))

  } else {
    column_df <- data.frame(columns = columns)
  }

  dir.create(dirname(opt$BC_output), showWarnings = FALSE, recursive = TRUE)
  write.table(column_df, file = opt$BC_output, sep = "\t", col.names = NA, quote = FALSE)
}