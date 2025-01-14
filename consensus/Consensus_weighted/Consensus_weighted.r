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
    c("--seed"),
    type = "integer", default = NULL,
    help = "seed for input"
  ),
  make_option(
    c("-b", "--base_clusterings"),
    type = "character", default = NULL,
    help = "Path to base-clustering ranking file"
  ), 
  make_option(
    c("--n_clusters"),
    type = "character", default = NULL,
    help = "Desired number of clusters in the consensus output"
  ),
  make_option(
    c("--n_bcs"),
    type = "integer", default = NULL,
    help = "Desired number of base clustering results feed into the algorithm"
  ),
  make_option(
    c("-o", "--output_file"),
    type = "character", default = NULL,
    help = "desired output file"
  ),
  make_option(
    c("--lambda"),
    type = "numeric", default = NULL,
    help = "regulation term coefficient"
  ),
  make_option(
    c("--jar_file"),
    type = "character", default = NULL,
    help = "path to the jar file to do leiden clustering"
  )
)

description <- "Calculate consensus for selected BCs"

opt_parser <- OptionParser(
  usage = description,
  option_list = option_list
)
opt <- parse_args(opt_parser)

# Use these filepaths as input
input_file <- opt$input_file
output_file <- opt$output_file
bc_file <- opt$base_clusterings
n_bcs <- ifelse(is.null(opt$n_bcs), 8, opt$n_bcs)
n_clust <- ifelse(is.null(opt$n_clusters), "7", opt$n_clusters)
seed <- opt$seed
jar_file <- opt$jar_file
lambda <- opt$lambda

set.seed(seed)

suppressPackageStartupMessages({
    library(Matrix)
    library(dplyr)
    library(future.apply)
})

###################### Define functions ######################
#' The adaptive weighted ensemble-based learning method to integrate the multiple binary spots similarity matrix
#'
#'
#' @importFrom parallel makeCluster stopCluster parApply
#' @importFrom abind abind
#'
#' @param Results.clustering a list contains all the results of individual similarity matrix. The elements of list is a matrix, spots * spots.
#' @param lambda hyper-parameter constrain the weight of individual methods for ensemble. If the parameter is set to NULL, then, we will adopt the value in our algorithm.
#' @param prob.quantile numeric of probabilities with values in [0,1]. Default setting is 0.5.
#' @param niter a positive integer represents the maximum number of updating algorithm. Default setting is 100.
#' @param epsilon a parameter represents the stop criterion.
#'
#' @return a list contains a matrix of the ensemble similarity of spots and a vector of the weight assigned to base results.
#'
#'@export

solve_ensemble <- function(Results.clustering, 
                          lambda = NULL, 
                          prob.quantile = 0.5,
                          niter = 100, 
                          epsilon = 1e-5,
                          verbose = FALSE){
  suppressPackageStartupMessages(require(future.apply))
  
  plan(multisession, workers = parallelly::availableCores() - 1)

  options(digits = 7)
  # Results.clustering <- Results.clustering.all[[1]]
  num.methods <- length(Results.clustering)
  num.spots <- nrow(Results.clustering[[1]])
  num.cell.type <- ncol(Results.clustering[[1]])

  ## initialization V by the mean of individual values
  w <- c(rep(1/num.methods, num.methods))
  H <-  Reduce("+", Map("*", Results.clustering, w))

  if(is.null(lambda)){
    cat("We will adpote a value for lambda in our algorithm...", "\n")
  }

  k <- 1

  while (k <= niter) {
    st <- proc.time()
    if(k == 1){
      loss_all_temp <- 0
      # Generate the first loss value 
      temp2 <-  future_sapply(Results.clustering, JSD_Matrix, Y = H)
      # Empricial estimation of lambda in the paper
      if(is.null(lambda)){
        lambda <- quantile(temp2, probs = prob.quantile)
      }
    }else{
      loss_all_temp <- loss_all
    }
    ##### update w
    temp2 <-  future_sapply(Results.clustering, JSD_Matrix, Y = H)
    w <- exp(-temp2/lambda)/sum(exp(-temp2/lambda))
    ##### update H
    H <-  Reduce("+", Map("*", Results.clustering, w))

    # Objective function loss
    loss_main <- sum(sapply(Results.clustering, JSD_Matrix, Y = H) * w)
    loss_entropy <- sum(w * log(w))
    loss_all <- loss_main + lambda * loss_entropy

    if(k == niter){
      cat("The method maybe not convergens, the algorithm need an larger max_epoches!", "\n")}

    # Stopping criteria
    diff_iter <- abs(loss_all - loss_all_temp)
    delta <- proc.time() - st 
    if (verbose){
      cat("iter: ", k, "loss_main: ", loss_main, "loss_entropy: ", loss_entropy,
          "loss_all: ", loss_all, "lambda: ", lambda, "diff:", 
          diff_iter, "epoch_time:", delta[3], "\n")
    }

    if(diff_iter < epsilon | k >= niter){
      break
      }
    k <- k + 1

  }
  colnames(H) <- colnames(Results.clustering[[1]])
  return(list(H = H, w = w))

}

JSD_Matrix <- function(X, Y, epi=1e-10){
  # Flatten the matrix, get probability, add a pseudo-count
  x <- pmax(as.vector(X)/sum(X), epi)
  y <- pmax(as.vector(Y)/sum(Y), epi)
  m <- (x + y)/2

  # Make the JSD bounded by 1 by using log2
  JSD <- sqrt(0.5 * (sum(x * log2(x / m)) +
                     sum(y * log2(y / m))))
  return(JSD)
}

#' Get binary similarity matrix from a cluster vector,
get_binary_matrix <- function(cluster_vector){
  suppressPackageStartupMessages(require(Matrix))
  N <- length(cluster_vector)

  pairs <- which(outer(cluster_vector, cluster_vector, FUN = "=="), arr.ind = TRUE)
  # Remove diagonal entries
  pairs <- pairs[pairs[,1] != pairs[,2], ]

  # Create a sparse matrix using the pairs of indices
  bm <- Matrix::sparseMatrix(
    i = pairs[,1],
    j = pairs[,2],
    x = rep(1, nrow(pairs)),
    dims = c(N, N)
  )

  return(bm)
}

#' Get cluster label from Leiden clustering of the consensus binary matrix.
#' Using also binary_search function to ensure the proper number of clusters
get_cluster_label <- function(binary_matrix,
                              n_clust_target,
                              resolution_update = 2,
                              resolution_init = 0.5,
                              num_rs = 100,
                              tolerance = 1e-5,
                              verbose=FALSE,
                              min_target=NULL, 
                              jar_file=NULL){
  suppressPackageStartupMessages(require(igraph))

  graph <- igraph::graph_from_adjacency_matrix(binary_matrix,
                                       mode = "upper",
                                       weighted = TRUE,
                                       diag = FALSE)
  write.table(as_data_frame(graph, what="edges") %>% sort("from"), 
              "edgeList.txt", sep="\t", 
              col.names=FALSE, row.names=FALSE)
  if (verbose){cat("Weighted neighborhood graph created... \n")}

  # Initialize boundaries
  lb <- rb <- NULL
  n_clust <- -1
  if (is.null(min_target)){
    min_target <- nrow(binary_matrix)/(10*n_clust_target)
  }
  print(min_target)

  get_clusters <- function(graph, resolution, mt=min_target, jar_file=NULL){
    if (!is.null(jar_file) && file.exists(jar_file)){
        system("cp %s networkanalysis-1.3.0.jar", jar_file)
    } else{
        if (!file.exists("networkanalysis-1.3.0.jar")){
        system("wget https://repo1.maven.org/maven2/nl/cwts/networkanalysis/1.3.0/networkanalysis-1.3.0.jar")
        }
    }
    system(sprintf("java -cp networkanalysis-1.3.0.jar nl.cwts.networkanalysis.run.RunNetworkClustering -w -m %1.0f -r %s -o cluster.txt edgeList.txt", mt, resolution))

    results <- read.table("cluster.txt", header=TRUE, row.names=NULL)[,2]

    return(results)
  }

  res <-  resolution_init
  result <- get_clusters(graph, resolution = res, jar_file = jar_file)
  # Adjust cluster_ids extraction per method
  n_clust <- length(unique(result))
  if (verbose){cat(sprintf("Boundary search starts..res = %s, n_clust=%s \n", res, n_clust))}
  if (n_clust > n_clust_target) {
    while (n_clust > n_clust_target && res > 1e-5) {
      rb <- res
      res <- res / resolution_update
      result <- get_clusters(graph, resolution = res)
      n_clust <- length(unique(result))
      if (verbose){cat(sprintf("Boundary search..lb = %s, rb = %s, n_clust=%s \n", res, rb, n_clust))}
    }
    lb <- res
  } else if (n_clust < n_clust_target) {
    while (n_clust < n_clust_target) {
      lb <- res
      res <- res * resolution_update
      result <- get_clusters(graph, resolution = res)
      n_clust <- length(unique(result))
      if (verbose){cat(sprintf("Boundary search..lb = %s, rb = %s, n_clust=%s \n", lb, res, n_clust))}
    }
    rb <- res
  }
  if (n_clust == n_clust_target) {lb = rb = res}

  i <- 0
  if (verbose){cat(sprintf("Boundary search done. lb = %s, rb = %s, res = %s, n_clust=%s \n", lb, rb, res, n_clust))}

  while ((rb - lb > tolerance || lb == rb) && i < num_rs) {
    mid <- sqrt(lb * rb)
    # message("Resolution: ", mid)
    result <- get_clusters(graph, resolution = mid)
    n_clust <- length(unique(result))
    min_clust_size <- min(table(result))
    if (verbose){cat(sprintf("iter: %s, res: %s, n_clust: %s, min_clust_size: %s \n", i,  mid, n_clust, min_clust_size))} # nolint
    if (n_clust == n_clust_target && min_clust_size >= min_target) break
    if (n_clust > n_clust_target) {
      rb <- mid
    } else if (n_clust < n_clust_target){
      lb <- mid
    } else if (n_clust == n_clust_target && min_clust_size < min_target){
      rb <- mid
      if (rb == lb){lb <- lb*0.9}
    }
    i <- i + 1
  }

  # Warning if target not met
  if (n_clust != n_clust_target) {
    warning(sprintf("Warning: n_clust = %d not found in binary search, return best approximation with res = %f and n_clust = %d. (rb = %f, lb = %f, i = %d)", n_clust_target, mid, n_clust, rb, lb, i))
  }

  cluster_vector <- result
  file.remove("cluster.txt", "edgeList.txt", "networkanalysis-1.3.0.jar")
  return(cluster_vector)
}

###################### Consensus calling begin ######################

label_df <- read.delim(input_file, stringsAsFactors = FALSE, row.names = 1, numerals="no.loss")
bc_list <- read.delim(bc_file, stringsAsFactors = FALSE, row.names = 1, numerals="no.loss")[[as.character(n_clust)]]
bc_list <- bc_list[!is.na(bc_list)]

if (length(bc_list) < n_bcs){
  warning(sprintf("Not enough (%s) base clusterings(BCs) are available, use %s BCs instead.", n_bcs, length(bc_list)))
}
bc_list <- bc_list[1:min(n_bcs, length(bc_list))]

label_selected <- label_df[, bc_list]


# Make sure all the clusters are ranked 1 to n without jumping (SOTIP)
label_selected <- apply(label_selected, 2, function(u){
    unique_labels <- sort(unique(u))
    if (all(unique_labels==seq_along(unique_labels))) {
        return(as.factor(u))
    } else {
        # Count occurrences of each number
        freq <- table(u)
        rank_map <- rank(-freq, ties.method = "first") # Negative for descending order
        new_vec <- rank_map[as.character(vec)]
        new_vec <- as.factor(as.numeric(new_vec))
        return(new_u)
    }
})
# Get binary matrix for individual clustering
binary_matrices <- lapply(label_selected, get_binary_matrix)

# Get consensus clustering from KL EnSDD method
ensemble_result <- solve_ensemble(Results.clustering = binary_matrices,
                                  lambda=lambda,
                                  verbose=FALSE)

# Resulting binary matrix
binary_ensemble <- ensemble_result$H
weight_vector <- ensemble_result$w

if (is.null(nclust)){
    nclust <- max(label_selected)
}

# Ensemble labelling
ensemble_label <- get_cluster_label(binary_ensemble,
                                    n_clust_target=nclust,
                                    jar_file = jar_file,
                                    verbose=FALSE)
ensemble_df <- data.frame(consensus_weighted=ensemble_label, row.names = row.names(label_selected))

dir.create(dirname(output_file), showWarnings = FALSE, recursive = TRUE)
# Save the results
write.table(ensemble_df, file = output_file, sep = "\t", col.names = NA, quote = FALSE)
