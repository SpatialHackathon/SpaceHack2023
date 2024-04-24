binary_search <- function(
    spe,
    do_clustering,
    extract_nclust,
    n_clust_target,
    resolution_update = 2,
    resolution_init = 1,
    resolution_boundaries=NULL,
    num_rs = 100,
    tolerance = 1e-3,
    ...) {

  # Initialize boundaries
  lb <- rb <- NULL
  n_clust <- -1

  if (!is.null(resolution_boundaries)){
    lb <- resolution_boundaries[1]
    rb <- resolution_boundaries[2]
  } else {
    res <-  resolution_init
    result <- do_clustering(spe, resolution = res, ...)
    # Adjust cluster_ids extraction per method
    n_clust <- extract_nclust(result)
    if (n_clust > n_clust_target) {
      while (n_clust > n_clust_target && res > 1e-5) {
        rb <- res
        res <- res/resolution_update
        result <- do_clustering(spe, resolution = res, ...)
        n_clust <- extract_nclust(result)
      }
      lb <- res
    } else if (n_clust < n_clust_target) {
      while (n_clust < n_clust_target) {
        lb <- res 
        res <- res*resolution_update
        result <- do_clustering(spe, resolution = res, ...)
        n_clust <- extract_nclust(result)
      }
      rb <- res
    }
    if (n_clust == n_clust_target) {lb = rb = res }
  }

  i <- 0
  while ((rb - lb > tolerance || lb == rb) && i < num_rs) {
    mid <- sqrt(lb * rb)
    message("Resolution: ", mid)
    result <- do_clustering(spe, resolution = mid, ...)
    n_clust <- extract_nclust(result)
    if (n_clust == n_clust_target || lb == rb) break
    if (n_clust > n_clust_target) {
      rb <- mid
    } else {
      lb <- mid
    }
    i <- i + 1
  }

  # Warning if target not met
  if (n_clust != n_clust_target) {
    warning(sprintf("Warning: n_clust = %d not found in binary search, return best approximation with res = %f and n_clust = %d. (rb = %f, lb = %f, i = %d)", n_clust_target, mid, n_clust, rb, lb, i))
  }
  return(result)
}