#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Zaira Seferbekova; wrote the code for SpaceFlow

import argparse

# TODO adjust description
parser = argparse.ArgumentParser(description="Method ...")

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

args = parser.parse_args()

from pathlib import Path

out_dir = Path(args.out_dir)

# Output files
label_file = out_dir / "domains.tsv"
embedding_file = out_dir / "embedding.tsv"

n_clusters = args.n_clusters
technology = args.technology
seed = args.seed

# Load config file
import json

with open(args.config, "r") as f:
    config = json.load(f)
    
def get_anndata(args):
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy.io import mmread
    from scipy.sparse import issparse
    from PIL import Image

    observations = pd.read_table(args.observations, index_col=0)
    features = pd.read_table(args.features, index_col=0)

    coordinates = (
        pd.read_table(args.coordinates, index_col=0)
        .loc[observations.index, :]
        .to_numpy()
    )

    adata = ad.AnnData(obs=observations, var=features, obsm={"spatial": coordinates})

    if args.matrix is not None:
        X = mmread(args.matrix)
        if issparse(X):
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


# Set the seed
import random
random.seed(seed)

import scanpy as sc
from SpaceFlow import SpaceFlow
import pandas as pd
import warnings
import torch

use_cuda = torch.cuda.is_available()
device = 1 if use_cuda else 0

# Filter genes expressed in < 3 cells
sc.pp.filter_genes(adata, min_cells=3)

# Create a SpaceFlow object and preprocess it
sf = SpaceFlow.SpaceFlow(adata=adata)
sf.preprocessing_data(n_top_genes=3000)

# Train the network
sf.train(spatial_regularization_strength=0.1, z_dim=50, lr=1e-3, epochs=1000, max_patience=50, min_stop=100, random_seed=seed, gpu=device, regularization_acceleration=True, edge_subset_sz=1000000, embedding_save_filepath=out_dir)

if config is not None:
    res = int(config['res'])
    nn = int(config['n_neighbours'])
else:
    res = 0.5
    nn = 15

# Raise a warning that clustering is based on resolution and not n_clusters
warnings.warn("The `n_clusters` parameter was not used; config['res'] used instead.")

# Segment the domains given the resolution
sf.segmentation(domain_label_save_filepath=out_dir, n_neighbors=nn, resolution=res)

label_df = pd.DataFrame(sf.domains)  # DataFrame with index (cell-id/barcode) and 1 column (label)
embedding_df = pd.DataFrame(sf.embedding, index=adata.obs_names) # DataFrame with index (cell-id/barcode) and n columns

## Write output
out_dir.mkdir(parents=True, exist_ok=True)

label_df.columns = ["label"]
label_df.to_csv(label_file, sep="\t", index_label="")

if embedding_df is not None:
    embedding_df.to_csv(embedding_file, sep="\t", index_label="")
