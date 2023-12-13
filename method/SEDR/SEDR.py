#!/usr/bin/env python

# Author_and_contribution: 
# Niklas Mueller-Boetticher; created template
# Jieran Sun; Implement SEDR method

import argparse

# TODO adjust description
parser = argparse.ArgumentParser(
    description="SEDR – Unsupervised spatially embedded deep representation of spatial transcriptomics. See https://doi.org/10.1101/2021.06.15.448542 for details")

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

## Session for code

# anndata input
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

    if args.img is not None:
        adata.uns["image"] = np.array(Image.open(args.img))
    else:
        adata.uns["image"] = None

    return adata

# Import packages for SEDR
import scanpy as sc
import numpy as np
import scipy as sp
import torch
import SEDR

import os
import warnings
warnings.filterwarnings('ignore')

# Get args
args = parser.parse_args()

# Def config
import json
with open(args.config, "r") as f:
    config = json.load(f)

# Set up output files
from pathlib import Path
out_dir = Path(args.out_dir)
label_file = out_dir / "domains.tsv"
embedding_file = out_dir / "embedding.tsv"

# Def vars
n_clusters = args.n_clusters
technology = args.technology
seed = args.seed

# Set SEED for SEDR
import random
random.seed(seed)
SEDR.fix_seed(seed)

# Load data
adata = get_anndata(args)
adata.var_names_make_unique()

# if dim-red is not provided, use default PCA dimRed
if "reduced_dimensions" not in adata.obsm_keys():
    from sklearn.decomposition import PCA 
    adata_X = PCA(n_components=200, 
                  random_state=seed).fit_transform(adata.X)
    adata.obsm['reduced_dimensions'] = adata_X

# Constructing neighborhood graphs if neighbors not provided
if "spatial_connectivities" in adata.obsp_keys():
    # Copy from source code in order for customization
    adj_m1 = adata.obsp["spatial_connectivities"]
    adj_m1 = sp.coo_matrix(adj_m1)

    # Store original adjacency matrix (without diagonal entries) for later
    adj_m1 = adj_m1 - sp.dia_matrix((adj_m1.diagonal()[np.newaxis, :], [0]), shape=adj_m1.shape)
    adj_m1.eliminate_zeros()

    # Some preprocessing
    adj_norm_m1 = SEDR.preprocess_graph(adj_m1)
    adj_m1 = adj_m1 + sp.eye(adj_m1.shape[0])

    adj_m1 = adj_m1.tocoo()
    shape = adj_m1.shape
    values = adj_m1.data
    indices = np.stack([adj_m1.row, adj_m1.col])
    adj_label_m1 = torch.sparse_coo_tensor(indices, values, shape)

    norm_m1 = adj_m1.shape[0] * adj_m1.shape[0] / float((adj_m1.shape[0] * adj_m1.shape[0] - adj_m1.sum()) * 2)

    graph_dict = {
        "adj_norm": adj_norm_m1,
        "adj_label": adj_label_m1.coalesce(),
        "norm_value": norm_m1
    }
    
else: 
    graph_dict = SEDR.graph_construction(adata, 
                                         n=12, 
                                         dmax=50, 
                                         mode=config['graph_mode'])
    

# Training SEDR
# device: using cpu or gpu (if avaliable)
# using_dec: boolean, whether to use the unsupervised deep embedded clustering (DEC) method to improve clustering results 
sedr_net = SEDR.Sedr(adata.obsm['reduced_dimensions'], 
                     graph_dict, 
                     mode='clustering', 
                     device=config["device"])

if config["using_dec"]:
    sedr_net.train_with_dec(N=1)
else:
    sedr_net.train_without_dec(N=1)
sedr_feat, _, _, _ = sedr_net.process()
# latent embedding
adata.obsm['SEDR'] = sedr_feat

# Clustering 
match config['cluster_method']:
    case "mclust":
        adata = SEDR.mclust_R(adata, 
                              n_clusters = n_clusters, 
                              use_rep='SEDR', 
                              key_added='SEDR', 
                              random_seed=seed
                             )
    case "louvain":
         adata = SEDR.louvain(adata, 
                              n_clusters = n_clusters, 
                              use_rep='SEDR', 
                              key_added='SEDR', 
                              random_seed=seed
                             )
    case "leiden":
         adata = SEDR.leiden(adata, 
                              n_clusters = n_clusters, 
                              use_rep='SEDR', 
                              key_added='SEDR', 
                              random_seed=seed
                             )
        
# Output dataframes
label_df = adata.obs[["SEDR"]]
embedding_df = adata.obsm['SEDR']

## Write output
out_dir.mkdir(parents=True, exist_ok=True)
label_df.columns = ["label"]
label_df.to_csv(label_file, sep="\t", index_label="")

if embedding_df is not None:
    embedding_df.to_csv(embedding_file, sep="\t", index_label="")
