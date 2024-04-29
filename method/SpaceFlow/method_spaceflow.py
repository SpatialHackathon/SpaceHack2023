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
    "--n_genes", help="Number of genes to use.", required=False, type=int
)
parser.add_argument(
    "--n_pcs", help="Number of PCs to use.", required=False, type=int
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

if args.n_genes is not None:
    n_genes = args.n_genes
else:
    n_genes = config["n_genes"] 

if args.n_pcs is not None:
    n_pcs = args.n_pcs
else:
    n_pcs = config["n_pcs"] 
    
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
import anndata as ad
from SpaceFlow import SpaceFlow
import pandas as pd
import warnings
import torch

# import res-n_clust tuning function
import sys
from pathlib import Path

# Add the parent directory of the current file to sys.path
method_dir = Path(__file__).resolve().parent.parent  # Navigate two levels up
sys.path.append(str(method_dir))

from search_res import binary_search

use_cuda = torch.cuda.is_available()
device = 1 if use_cuda else 0

# adata.write_h5ad("adata.h5ad")
n_vars = config["n_vars"]
nn = config["n_neighbours"]

# Create a SpaceFlow object 
sc.pp.filter_genes(adata, min_cells=1)
sf = SpaceFlow.SpaceFlow(adata=adata)

# modified from sf.preprocessing_data (https://github.com/hongleir/SpaceFlow/blob/master/SpaceFlow/SpaceFlow.py#L110)
def preprocessing_data(self, n_top_genes=None, n_comps = 50, n_neighbors=10):
    adata = self.adata
    if not adata:
        print("No annData object found, please run SpaceFlow.SpaceFlow(expr_data, spatial_locs) first!")
        return
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor='cell_ranger', subset=True)
    sc.pp.pca(adata, n_comps)
    spatial_locs = adata.obsm['spatial']
    spatial_graph = self.graph_alpha(spatial_locs, n_neighbors=n_neighbors)
    self.adata_preprocessed = adata
    self.spatial_graph = spatial_graph

preprocessing_data(n_top_genes = min(adata.n_vars, n_genes), n_comps = n_pcs)

# Train the network
sf.train(spatial_regularization_strength=0.1, z_dim=50, lr=1e-3, epochs=1000, 
         max_patience=50, min_stop=100, random_seed=seed, gpu=device, regularization_acceleration=True, 
         edge_subset_sz=1000000, embedding_save_filepath=embedding_file)

# Raise a warning that clustering is based on resolution and not n_clusters
# warnings.warn("The `n_clusters` parameter was not used; config['res'] used instead.")

# Segment the domains given the resolution
embedding_adata = ad.AnnData(sf.embedding)
sc.pp.neighbors(embedding_adata, n_neighbors=nn, use_rep="X")
label_df = binary_search(embedding_adata, n_clust_target=n_clusters, method="leiden", seed = seed)
# sf.segmentation(domain_label_save_filepath=label_file, n_neighbors=nn, resolution=res)

# label_df = pd.DataFrame(sf.domains)  # DataFrame with index (cell-id/barcode) and 1 column (label)
embedding_df = pd.DataFrame(sf.embedding, index=adata.obs_names) # DataFrame with index (cell-id/barcode) and n columns

## Write output
out_dir.mkdir(parents=True, exist_ok=True)

label_df.columns = ["label"]
label_df.index = adata.obs.index
label_df.to_csv(label_file, sep="\t", index_label="")

if embedding_df is not None:
    embedding_df.to_csv(embedding_file, sep="\t", index_label="")