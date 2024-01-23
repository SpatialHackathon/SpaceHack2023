#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Jieran Sun; implemented method cellCharter

import argparse

# TODO adjust description
parser = argparse.ArgumentParser(description="Method CellCharter: https://doi.org/10.1038/s41588-023-01588-4")

parser.add_argument(
    "-c", "--coordinates", help="Path to coordinates (as tsv).", required=True
)
parser.add_argument(
    "-m", "--matrix", help="Path to (transformed) counts (as mtx).", required=False
)
parser.add_argument(
    "-f", "--features", help="Path to features (as tsv).", required=True
)
parser.add_argument(
    "-o", "--observations", help="Path to observations (as tsv).", required=True
)
parser.add_argument(
    "-n",
    "--neighbors",
    help="Path to neighbor definitions. Square matrix (not necessarily symmetric) where each row contains the neighbors of this observation (as mtx).",
    required=False,
)
parser.add_argument("-d", "--out_dir", help="Output directory.", required=True)
parser.add_argument(
    "--dim_red",
    help="Reduced dimensionality representation (e.g. PCA).",
    required=False,
)
parser.add_argument("--image", help="Path to H&E staining.", required=False)
parser.add_argument(
    "--n_clusters", help="Number of clusters to return.", required=True, type=int
)
parser.add_argument(
    "--technology",
    help="The technology of the dataset (Visium, ST, imaging-based).",
    required=True,
)
parser.add_argument(
    "--seed", help="Seed to use for random operations.", required=True, type=int
)
parser.add_argument(
    "--config",
    help="Optional config file (json) used to pass additional parameters.",
    required=False,
)

################# Get args #######################

args = parser.parse_args()

n_clusters = args.n_clusters
technology = args.technology
seed = args.seed
# Def config
import json
with open(args.config, "r") as f:
    config = json.load(f)

# Set up output files
from pathlib import Path
out_dir = Path(args.out_dir)
label_file = out_dir / "domains.tsv"
embedding_file = out_dir / "embedding.tsv"

################# Define functions #################
def get_anndata(args):
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy as sp
    from PIL import Image

    observations = pd.read_table(args.observations, index_col=0)
    features = pd.read_table(args.features, index_col=0)

    coordinates = (
        pd.read_table(args.coordinates, index_col=0)
        .loc[observations.index, :]
        .to_numpy()
    )

    adata = ad.AnnData(obs=observations, 
                       var=features, 
                       obsm={"spatial": coordinates})

    if args.matrix is not None:
        X = sp.io.mmread(args.matrix)
        if sp.sparse.issparse(X):
            X = X.tocsr()
        adata.X = X

    if args.neighbors is not None:
        adata.obsp["spatial_connectivities"] = sp.io.mmread(args.neighbors).T.tocsr()

    # Filter by selected samples/features
    if "selected" in adata.obs.columns:
        adata = adata[observations["selected"].astype(bool), :]
    if "selected" in adata.var.columns:
        adata = adata[:, features["selected"].astype(bool)]

    if args.dim_red is not None:
        adata.obsm["reduced_dimensions"] = (
            pd.read_table(args.dim_red, index_col=0).loc[adata.obs_names].to_numpy()
        )

    if args.image is not None:
        adata.uns["image"] = np.array(Image.open(args.image))
    else:
        adata.uns["image"] = None

    return adata

## Your code goes here
import cellcharter as cc
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import euclidean_distances
from scipy.sparse import csr_matrix
from itertools import chain
import random

random.seed(seed)

# Load data
adata = get_anndata(args)
adata.var_names_make_unique()

# Only calculates distance and remove edges if it's not Visium(grid structure)
if technology != "Visium":
    
    coords = adata.obsm["spatial"]
    N = coords.shape[0]
    Adj = adata.obsp["spatial_connectivities"] 
    
    # Calculate euclidean distances
    dists = np.array(list(chain(*(
        euclidean_distances(coords[Adj.indices[Adj.indptr[i] : Adj.indptr[i + 1]], :], 
                            coords[np.newaxis, i, :])
        for i in range(N)
        if len(Adj.indices[Adj.indptr[i] : Adj.indptr[i + 1]])
    )))).squeeze()
    
    Dst = csr_matrix((dists, Adj.indices, Adj.indptr), shape=(N, N))
    Dst.setdiag(0.0)
    
    # Fitting into the respective data slot
    adata.obsp["spatial_distances"] = Dst
    adata.uns["spatial_neighbors"] = {
            "connectivities_key": "spatial_connectivities",
            "distances_key": "spatial_distances",
            "params": {"coord_type": "generic", 
                       "delaunay": True, 
                       "transform": None,
                       "technology": technology},
        }
    
    # Remove links between cells at a distance bigger than 99% of all positive distances.
    cc.gr.remove_long_links(adata, distance_percentile = 0.99)

else:
    adata.uns["spatial_neighbors"] = {
        "connectivities_key": "spatial_connectivities",
        "distances_key": None,
        "params": {"coord_type": "grid", 
                   "delaunay": False, 
                   "transform": None,
                   "technology": technology},
    }


# Aggregate neighborhoods
cc.gr.aggregate_neighbors(adata, 
                          n_layers=config["n_layers"], 
                          aggregations=config["aggregations"],
                          use_rep='reduced_dimensions', 
                          out_key='X_cellcharter', 
                          sample_key=None)

# clustering
model = cc.tl.Cluster(
    n_clusters=n_clusters, 
    random_state=seed,
    convergence_tolerance=config["convergence_tolerance"], 
    covariance_regularization=config["covariance_regularization"]
    # if running on gpu
    #trainer_params=dict(accelerator='gpu', devices=1)
    )

model.fit(adata, use_rep='X_cellcharter')

# label prediction
label_df = pd.DataFrame({"label" : model.predict(adata, use_rep='X_cellcharter')}) 
label_df.index = adata.obs.index

embedding_df = pd.DataFrame(adata.obsm["X_cellcharter"])
embedding_df.index = adata.obs_names

## Write output
out_dir.mkdir(parents=True, exist_ok=True)

label_df.to_csv(label_file, sep="\t", index_label="")
embedding_df.to_csv(embedding_file, sep="\t", index_label="")
