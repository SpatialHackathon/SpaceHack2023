#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created script

import argparse

parser = argparse.ArgumentParser(
    description="""SpaGCN (https://www.nature.com/articles/s41592-021-01255-8)
Requirements: Visium or ST. Images can be used."""
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
    required=True,
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
technology = args.technology
seed = args.seed


## Your code goes here
import json

with open(args.config, "r") as f:
    config = json.load(f)

if technology not in ["Visium", "ST"]:
    config["refine"] = False

# TODO how to determine beta
beta = 49


import random

import numpy as np
import pandas as pd
import scanpy as sc
import SpaGCN as spg
import torch


def get_anndata(args):
    import anndata as ad
    import scipy as sp
    from PIL import Image

    X = sp.io.mmread(args.matrix)
    if sp.sparse.issparse(X):
        X = X.tocsr()

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

    if args.image is not None:
        Image.MAX_IMAGE_PIXELS = None
        adata.uns["image"] = np.array(Image.open(args.image))
    else:
        adata.uns["image"] = None

    return adata

# Set seed
random.seed(seed)
torch.manual_seed(seed)
np.random.seed(seed)

# reference: https://github.com/jianhuupenn/SpaGCN/blob/master/tutorial/tutorial.md#5-spatial-domain-detection-using-spagcn

adata = get_anndata(args)
adata.var_names_make_unique()
# adata.write_h5ad("adata.h5ad")

spg.prefilter_genes(adata, min_cells=3)
spg.prefilter_specialgenes(adata)

# Normalize and take log for UMI
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

if technology in ["Visium", "ST"]:
    if adata.uns["image"] is not None:
        adj = spg.calculate_adj_matrix(
            adata.obs["col"],
            adata.obs["row"],
            adata.obsm["spatial_pixel"][:, 0],
            adata.obsm["spatial_pixel"][:, 1],
            image=adata.uns["image"],
            alpha=config["alpha"],
            beta=beta,
            histology=True,
        )
    else:
        adj = spg.calculate_adj_matrix(
            adata.obs["col"], adata.obs["row"], histology=False
        )
else:
    adj = spg.calculate_adj_matrix(
        adata.obsm["spatial_pixel"][:, 0],
        adata.obsm["spatial_pixel"][:, 1],
        histology=False,
    )


clf = spg.SpaGCN()

# Find the l value given p
l = spg.search_l(config["p"], adj)
clf.set_l(l)


# Run
if config["method"] == "louvain":
    # Search for suitable resolution
    res = spg.search_res(
        adata, adj, l, n_clusters, r_seed=seed, t_seed=seed, n_seed=seed
    )
    clf.train(
        adata,
        adj,
        init_spa=True,
        init=config["method"],
        res=res,
        n_neighbors=config["n_neighbors"],
        num_pcs=config["n_pcs"],
    )
else:
    clf.train(
        adata,
        adj,
        init_spa=True,
        init=config["method"],
        n_clusters=n_clusters,
        num_pcs=config["n_pcs"],
    )
y_pred, prob = clf.predict()

adata.obs["cluster"] = pd.Series(y_pred, index=adata.obs_names, dtype="category")

if technology in ["Visium", "ST"] and config["refine"]:
    # Do cluster refinement(optional)
    adj_2d = spg.calculate_adj_matrix(
        adata.obs["col"], adata.obs["row"], histology=False
    )

    shape = "hexagon" if technology == "Visium" else "square"
    refined_pred = spg.refine(
        sample_id=adata.obs_names.tolist(),
        pred=adata.obs["cluster"].tolist(),
        dis=adj_2d,
        shape=shape,
    )

    adata.obs["refined_cluster"] = pd.Series(
        refined_pred, index=adata.obs_names, dtype="category"
    )

    label_df = adata.obs[["refined_cluster"]]

else:
    label_df = adata.obs[["cluster"]]


## Write output
out_dir.mkdir(parents=True, exist_ok=True)

label_df.columns = ["label"]
label_df.to_csv(label_file, sep="\t", index_label="")
