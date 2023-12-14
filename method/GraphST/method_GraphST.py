#!/usr/bin/env python

# Author_and_contribution: 
# Niklas Mueller-Boetticher; created template
# Peiying Cai; implemented method

import argparse
from ast import If

parser = argparse.ArgumentParser(
    description="""GraphST (https://www.nature.com/articles/s41467-023-36796-3)
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
seed = args.seed

# Map technology spelling
def map_technology(input_technology):
    technology_mapping = {
        "Visium": "10X",
        "Stereo-seq": "Stereo",
        "Slide-seq": "Slide"
    }
    return technology_mapping.get(input_technology, None)

technology = map_technology(str(args.technology))
if technology is None:
    raise Exception(
        f"Invalid technology. GraphST only supports 10X Visium, Stereo-seq, and Slide-seq/Slide-seqV2 not {args.technology}. "
        )
    
# TODO use the default settings or place it into the config.
# start=0.1
# end=3.0
# increment=0.01

import json

with open(args.config, "r") as f:
    config = json.load(f)

import random

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn import metrics
import torch
from GraphST import GraphST

# Run device, by default, the package is implemented on 'cpu'. The author recommend using GPU.
device = torch.device('cuda:1' if torch.cuda.is_available() else 'cpu')

def get_anndata(args):
    import anndata as ad
    import scipy as sp

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
        X=X, obs=observations, var=features, obsm={"spatial": coordinates}
    )

    return adata


adata = get_anndata(args)
adata.var_names_make_unique()

# Set seed
random.seed(seed)
torch.manual_seed(seed)
np.random.seed(seed)

# Train
# define model
model = GraphST.GraphST(adata, datatype=technology, device=device)

# train model
adata = model.train()

# Refinement:
# This step is not recommended for ST data with fine-grained domains, Stereo-seq, and Slide-seqV2
# Radius: speficy the number of neighbors considered during refinement
# If not refinement is intended, set it as the default setting in GraphST, but it is not going to be used. 
if config["refine"]:
    radius = config["radius"]
else:
    radius = 50
    

# Run
from GraphST.utils import clustering

clustering(adata, 
           n_clusters=n_clusters, 
           radius=radius, 
           method=config["method"], 
           refinement=config["refine"]
           )

label_df = adata.obs[["domain"]]

## Write output
out_dir.mkdir(parents=True, exist_ok=True)

label_df.columns = ["label"]
label_df.to_csv(label_file, sep="\t", index_label="")
