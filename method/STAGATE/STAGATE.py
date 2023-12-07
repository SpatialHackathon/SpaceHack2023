#!/usr/bin/env python

# Author_and_contribution: 
# Niklas Mueller-Boetticher; created template
# Peiying Cai; implemeted method

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
parser.add_argument("--R_dir", help="Path to R.", required=False)
parser.add_argument("--rpy2_dir", help="Path to the side package: rpy2.", required=False)
# Example: 'D:\Program Files\R\R-4.3.1' and 'D:\ProgramData\Anaconda3\Lib\site-packages\rpy2'

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

import os
if config["method"] in ['mclust']:
    if args.R_dir is None or args.rpy2_dir is None:
        raise Exception(
            f"The location of R and rpy2 is necessary for the mclust algorithm. Please specify the paths to the installed R and rpy2."
            )
    else:
        os.environ['R_HOME'] = args.R_dir
        os.environ['R_USER'] = args.rpy2_dir


import random
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import sys
import STAGATE
import torch
# Use tensorflow or pyG
import tensorflow.compat.v1 as tf
tf.compat.v1.disable_eager_execution()
# the location of R (used for the mclust clustering)

def get_anndata(args):
    import anndata as ad
    import scipy as sp

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
        X=X, obs=observations, var=features, obsm={"spatial_pixel": coordinates}
    )

    return adata


adata = get_anndata(args)
adata.var_names_make_unique()

# Set seed
random.seed(seed)
torch.manual_seed(seed)
np.random.seed(seed)

# Construct the spatial network
STAGATE.Cal_Spatial_Net(adata, rad_cutoff=config["rad"])
STAGATE.Stats_Spatial_Net(adata)

# Run
# alpha: the weight of the cell type-aware spatial network
# pre-resolution: the resolution parameter of pre-clustering
# with cell type-aware module:
# adata = STAGATE.train_STAGATE(adata, alpha=0.5, pre_resolution=0.2,
#                              n_epochs=1000, save_attention=True)
adata = STAGATE.train_STAGATE(adata, alpha=0, random_seed=seed)
sc.pp.neighbors(adata, use_rep='STAGATE')
sc.tl.umap(adata)

if config["method"] == "mclust":
    import os
    adata = STAGATE.mclust_R(adata, used_obsm='STAGATE', num_cluster=n_clusters, random_seed=seed)
    label_df = adata.obs['mclust']
elif config["method"] == "louvain":
    sc.tl.louvain(adata, resolution=config["res"])
    label_df = adata.obs['louvain']

## Write output
out_dir.mkdir(parents=True, exist_ok=True)

label_df.columns = ["label"]
label_df.to_csv(label_file, sep="\t", index_label="")

embedding_df = adata.obsm['STAGATE']
if embedding_df is not None:
    embedding_df.to_csv(embedding_file, sep="\t", index_label="")