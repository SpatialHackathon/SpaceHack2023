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
    "--n_genes", help="Number of genes to use.", required=False, type=int
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

# ... or AnnData if you want
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
        .iloc[:, 0:2]
        .to_numpy()
    )

    adata = ad.AnnData(obs=observations, var=features, obsm={"spatial": coordinates})

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
        adata.uns["image"] = np.array(Image.open(args.img))
    else:
        adata.uns["image"] = None

    return adata

adata = get_anndata(args)
adata.var_names_make_unique()

# Set seed
random.seed(seed)
tf.random.get_seed(seed)
np.random.seed(seed)

import STAGATE as sg
from scipy.spatial.distance import euclidean

with open(args.config, "r") as f:
    config = json.load(f)

if args.n_genes is not None:
    n_genes = args.n_genes
else:
    n_genes = config["n_genes"] # default setting: 3000

method = config["method"]
model = config["model"]
rad_cutoff = config["rad_cutoff"] if model == "Radius" else None
k_cutoff = config["k_cutoff"] if model == "KNN" else None
min_cells = config["min_cells"]

sc.pp.filter_genes(adata, min_cells=min_cells)
#Normalization
if adata.n_vars > n_genes:
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=n_genes)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sg.Cal_Spatial_Net(adata, rad_cutoff=rad_cutoff, k_cutoff = k_cutoff, model = model)

sg.Stats_Spatial_Net(adata)

# Run

adata = sg.train_STAGATE(adata, random_seed=seed)
sc.pp.neighbors(adata, use_rep='STAGATE')
if method == "mclust":
    adata = sg.mclust_R(adata, used_obsm='STAGATE', num_cluster=n_clusters, random_seed=seed)
    label_df = adata.obs[["mclust"]]
else:
    label_df = binary_search(adata=adata, n_clust_target=n_clusters, method=method, seed=seed)

## Write output
out_dir.mkdir(parents=True, exist_ok=True)
label_df.columns = ["label"]
label_df.to_csv(label_file, sep="\t", index_label="")

# Convert the NumPy array to Pandas DataFrame with row names
embedding_df = pd.DataFrame(adata.obsm['STAGATE'], index=adata.obs_names)
if embedding_df is not None:
    embedding_df.to_csv(embedding_file, sep="\t", index_label="")
