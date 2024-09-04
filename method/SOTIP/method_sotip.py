#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Zaira Seferbekova; wrote the code for SOTIP

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

# Load config file
import json

with open(args.config, "r") as f:
    config = json.load(f)
    
# Output files
label_file = out_dir / "domains.tsv"
embedding_file = out_dir / "embedding.tsv"

n_clusters = args.n_clusters
technology = args.technology
seed = args.seed

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

from sotip import *
import scanpy as sc

if "n_pcs" not in config.keys():
    config["n_pcs"] = 50 # default
  
n_pcs = args.n_pcs if args.n_pcs is not None else config["n_pcs"]

# Process the data with scanpy routine
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata,n_comps=n_pcs)
sc.pp.neighbors(adata, use_rep='X_pca')
sc.tl.umap(adata)

# Find resolution for the given n_clusters
if "res" in config.keys():
    res_recom = config["res"]
else:
    res_recom = search_res(adata, n_clusters, start=0.1, step=0.1, tol=5e-3, max_run=10)

sc.tl.leiden(adata,resolution=res_recom)
leiden_n_clusters = len(adata.obs['leiden'].cat.categories)
# If the number of clusters is less than expected, log a message and stop the execution
import sys
# import logging
# logging.error
if leiden_n_clusters < n_clusters:
    print(f"Number of clusters obtained ({leiden_n_clusters}) is less than the expected number ({n_clusters}). "
          f"Consider setting higher resolution values in config['res'] to get more clusters.")
    sys.exit(1)

# ME size 
knn = 10
n_neighbors = config['n_neighbours']
    
# Order of cluster label for ME representation (adata.obsm['ME'])
ME_var_names_np_unique = np.array(adata.obs['leiden'].cat.categories) 

# Add a ME obsm for adata
MED(adata, use_cls='leiden', nn=knn, copy=False, ME_var_names_np_unique=ME_var_names_np_unique, spatial_var='spatial') 

# Compute the topological structure with paga
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, add_pos=True, show=False)

# Use the connectivities between cell clusters to guide the graph distance computation
gd = get_ground_distance(adata,method='paga_guided_umap',cls_key='leiden')

# Add a X_ME_EMD_mat obsm to adata_phEMD
adata_phEMD = MED_phEMD_mp(
    adata.copy(),         # the used anndata
    GD_method='paga_guided_umap',  # use CGMGD as ground distance
    MED_knn=knn,          # ME size is set consistently as 10
    CT_obs='leiden',       # use leiden cluster label
    ifspatialplot=False,  # do not show intermediate result plot
    OT_method='pyemd',    # use pyemd to compute EMD as ME distance
    ME_precompyted=True,  # use precomputed ME representation (already computed in step 1)
    GD_precomputed=True,  # use precomputed ground distance (already computed in step 2)
)

adata.obsp['ME_EMD_mat'] = adata_phEMD.obsm['X_ME_EMD_mat']

# Compute the MEG, each node is a ME, edge is the connectivity
# sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="reduced_dimensions") # related to issues#208
sc.pp.neighbors(adata, n_neighbors=n_neighbors)
knn_indices, knn_dists, forest = sc.neighbors.compute_neighbors_umap(adata_phEMD.obsm['X_ME_EMD_mat'], n_neighbors=n_neighbors, metric='precomputed')
adata.obsp['distances'], adata.obsp['connectivities'] = sc.neighbors._compute_connectivities_umap(
    knn_indices,
    knn_dists,
    adata.shape[0],
    n_neighbors, 
)

# Set the ME graph's associated information (connectivity matrix, distance matrix) to neighbors_EMD
adata.uns['neighbors_EMD'] = adata.uns['neighbors'].copy()

# Use the computed MEG as input of umap and leiden clustering
sc.tl.umap(adata,neighbors_key='neighbors_EMD')
sc.tl.leiden(adata,neighbors_key='neighbors_EMD', key_added='leiden_EMD')

#use_rep = 'leiden_EMD'
#if len(adata.obs[use_rep].cat.categories) > n_clusters:
# Merge regions according to MEG connectivities
sc.tl.paga(adata, groups='leiden_EMD', neighbors_key='neighbors_EMD')
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scanpy as sc
import palettable
import anndata as ad
import scanpy.external as sce
import sys
from scipy.stats import *
from sklearn.metrics import *

# cls_key = 'leiden_EMD'
# neighbor_key = 'neighbors_EMD'
# thresh = 0.5

merge_cls_paga(adata, thresh=0, min_cls=n_clusters, paga_plot=False)
use_rep = 'leiden_EMD_merge'

embedding_df = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names) # optional, DataFrame with index (cell-id/barcode) and n columns
label_df = adata.obs[[use_rep]]  # DataFrame with index (cell-id/barcode) and 1 column (label)

## Write output
out_dir.mkdir(parents=True, exist_ok=True)

label_df.columns = ["label"]
label_df.to_csv(label_file, sep="\t", index_label="")

if embedding_df is not None:
    embedding_df.to_csv(embedding_file, sep="\t", index_label="")
