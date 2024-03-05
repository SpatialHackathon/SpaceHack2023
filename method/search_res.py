def binary_search(
    adata, n_clust_target,
    method="louvain",
    resolution_boundaries=None,
    resolution_init=1,
    resolution_update=2,
    num_rs=1e2,
    tolerance = 1e-3,
    seed = 2023,
):
    """
    Uses binary search to find the resolution parameter that results in the target number of clusters.

    Parameters
    ------------
    adata (Anndata)
        AnnData object for clustering. Should have neighbor graphs already.

    target_n_clusters (int)
        The desired number of clusters.

    method (str)
        Cluster function of choice. Can be "louvain" or "leiden".

    resolution_boundary (list)
        a list defining the boundary for resolution search:

            - If `None`, the function will find a rough boundary that contains the target resolution using `resolution_init` and `resolution_update`;
            - If defined, follows [`left_boundary`, `right_boundary`];

    resolution_init (float)
        initial resolution to start the search with.

    resolution_update(float)
        a positive number for the scale of coarse boundary definition. Only used when `resolution_boundary = None`.

    num_rs (int):
        Highest number of iteration for search.

    tolerance (float):
        Smallest gap between search boundaries. search will stop if the interval is smaller than tolerance.

    seed (int):
        seed for clustering.

    Returns
    ------------
    y: a pandas dataframe with one clumn denoting the clustering results from best resolution and with barcode as index.
    """
    import scanpy as sc
    import numpy as np
    import warnings

    y = None

    def do_clustering(res):
        getattr(sc.tl, method)(adata, resolution=res, random_state=seed)
        y = adata.obs[[method]].astype(int)
        n_clust = len(np.unique(y))
        return y, n_clust

    lb = rb = None
    n_clust = -1
    if resolution_boundaries is not None:
        lb, rb = resolution_boundaries
    else:
        res = resolution_init
        y, n_clust = do_clustering(res)
        # coarse search for the boundary containing n_clust_target
        if n_clust > n_clust_target:
            while n_clust > n_clust_target and res > 1e-4:
                rb = res
                res /= resolution_update
                y, n_clust = do_clustering(res)
            lb = res
        elif n_clust < n_clust_target:
            while n_clust < n_clust_target:
                lb = res
                res *= resolution_update
                y, n_clust = do_clustering(res)
            rb = res
        if n_clust == n_clust_target: lb = rb = res

    i = 0
    while (rb - lb > tolerance or lb == rb) and i < num_rs:
        mid = (lb * rb) ** .5
        y, n_clust = do_clustering(mid)
        if n_clust == n_clust_target or lb == rb: break
        if n_clust > n_clust_target: rb = mid
        else: lb = mid
        i += 1

    # Check if the situation is met
    if n_clust != n_clust_target:
        warnings.warn(f"Warning: n_clust = {n_clust_target} not found in binary search, \
        return best proximation with res = {mid} and \
        n_clust = {n_clust}. (rb = {rb}, lb = {lb}, i = {i})")

    return y