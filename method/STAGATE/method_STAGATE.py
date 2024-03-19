#!/usr/bin/env python

# Author_and_contribution: 
# Niklas Mueller-Boetticher; created template
# Peiying Cai; implemented method

import argparse


parser = argparse.ArgumentParser(
    description="""STAGATE (https://www.nature.com/articles/s41467-022-29439-6)
"""
)

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
    help="The technology of the dataset (Visium, ST, ...).",
    required=False,
)
parser.add_argument(
    "--seed", help="Seed to use for random operations.", required=True, type=int
)
parser.add_argument(
    "--config",
    help="Optional config file used to pass additional parameters.",
    required=False,
)

args = parser.parse_args()

from pathlib import Path

out_dir = Path(args.out_dir)

# Output files
label_file = out_dir / "domains.tsv"
embedding_file = out_dir / "embedding.tsv"

n_clusters = args.n_clusters
seed = args.seed


## Your code goes here
import json

with open(args.config, "r") as f:
    config = json.load(f)

import random
import numpy as np
import pandas as pd
import scanpy as sc
# Use tensorflow or pyG
import tensorflow.compat.v1 as tf
tf.compat.v1.disable_eager_execution()
# the location of R (used for the mclust clustering)
import scipy as sp

# import res-n_clust tuning function
import sys
from pathlib import Path

# Add the parent directory of the current file to sys.path
method_dir = Path(__file__).resolve().parent.parent  # Navigate two levels up
sys.path.append(str(method_dir))

from search_res import binary_search

def get_anndata(args):
    import anndata as ad
    X = sp.io.mmread(args.matrix)
    if sp.sparse.issparse(X):
        X = sp.sparse.csr_matrix(X, dtype= 'float32')
    observations = pd.read_table(args.observations, index_col=0)
    features = pd.read_table(args.features, index_col=0)
    # Filter
    if "selected" in observations.columns:
        X = X[observations["selected"].to_numpy().nonzero()[0], :]
        observations = observations.loc[lambda df: df["selected"]]
    if "selected" in features.columns:
        X = X[:, features["selected"].to_numpy().nonzero()[0]]
        features = features.loc[lambda df: df["selected"]]
    coordinates = (
        pd.read_table(args.coordinates, index_col=0)
        .loc[observations.index, :]
        .to_numpy()
    )
    adata = ad.AnnData(
        X=X, obs=observations, var=features, obsm={"spatial": coordinates}
    )
    return adata


adata = get_anndata(args)
adata.var_names_make_unique()

# Set seed
random.seed(seed)
tf.random.get_seed(seed)
np.random.seed(seed)

import STAGATE as sg
from scipy.spatial.distance import euclidean
# To skip the neighbor computation step in STAGATE
if args.neighbors is not None:
    neighbors = sp.io.mmread(args.neighbors).T.tocoo()
    row_indices, col_indices = neighbors.row, neighbors.col
    coor = adata.obsm['spatial']
    cell_ids = adata.obs_names
    Spatial_Net = pd.DataFrame({'Cell1': cell_ids[row_indices],
                                'Cell2': cell_ids[col_indices],
                                'Distance': [euclidean(coor[i], coor[j]) for i, j in zip(row_indices, col_indices)]})
    Spatial_Net = Spatial_Net.loc[Spatial_Net['Distance']>0,]
    adata.uns['Spatial_Net'] = Spatial_Net
elif config["model"] == 'Radius':
    sg.Cal_Spatial_Net(adata, rad_cutoff=config["rad_cutoff"], model=config["model"])
elif config["model"] == 'KNN':
    sg.Cal_Spatial_Net(adata, k_cutoff=config["k_cutoff"], model=config["model"])

sg.Stats_Spatial_Net(adata)

# Run
# alpha: the weight of the cell type-aware spatial network
# pre-resolution: the resolution parameter of pre-clustering
# with cell type-aware module:
# adata = sg.train_STAGATE(adata, alpha=0.5, pre_resolution=0.2,
#                              n_epochs=1000, save_attention=True)
adata = sg.train_STAGATE(adata, alpha=0, random_seed=seed)
sc.pp.neighbors(adata, use_rep='STAGATE')
sc.tl.umap(adata)

if config["method"] == "mclust":
    adata = sg.mclust_R(adata, used_obsm='STAGATE', num_cluster=n_clusters, random_seed=seed)
    label_df = adata.obs[["mclust"]]
elif config["method"] == "louvain":
    label_df = binary_search(adata, n_clust_target=n_clusters, method="louvain", seed = seed)
    #sc.tl.louvain(adata, resolution=config["res"])
    #label_df = adata.obs[["louvain"]]

## Write output
out_dir.mkdir(parents=True, exist_ok=True)

label_df.columns = ["label"]
label_df.to_csv(label_file, sep="\t", index_label="")

# Convert the NumPy array to Pandas DataFrame with row names
embedding_df = pd.DataFrame(adata.obsm['STAGATE'], index=adata.obs_names)
if embedding_df is not None:
    embedding_df.to_csv(embedding_file, sep="\t", index_label="")
