#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Zaira Seferbekova; wrote the code for SCAN-IT

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

    if args.image is not None:
        adata.uns["image"] = np.array(Image.open(args.img))
    else:
        adata.uns["image"] = None

    return adata


adata = get_anndata(args)


# Set the seed
import random
random.seed(seed)

import scanit
import scanpy as sc
import pandas as pd

# Construct the spatial graph
scanit.tl.spatial_graph(adata, method='alpha shape', alpha_n_layer=2, knn_n_neighbors=5)

# Generate low dimentional embedding (saved to X_scanit)
scanit.tl.spatial_representation(adata, n_h=30, n_epoch=2000, lr=0.001, device='cpu', n_consensus=1, projection='mds', python_seed=seed, torch_seed=seed, numpy_seed=seed)

# Construct a NN graph based on the embedding
sc.pp.neighbors(adata, use_rep='X_scanit', n_neighbors=15)

# Run clustering and adjust resolution until the n_clusters is reached
res=0.5
step=0.1
old_num = -1  # Set an initial value for old_num

while old_num != n_clusters:
    step_sign = 1 if (old_num < n_clusters) else -1
    
    sc.tl.leiden(adata, resolution=res)
    new_num = adata.obs['leiden'].nunique()
    print("Res = ", res, "Num of clusters = ", new_num)
    
    if new_num == n_clusters:
        print("Recommended res = ", str(res))
        break  # Break the loop if the desired number of clusters is reached
    
    if step_sign == 1:
        res += step
    else:
        res -= step
        step /= 2  # Reduce step size when changing direction
    
    old_num = new_num

label_df = adata.obs[['leiden']]  # DataFrame with index (cell-id/barcode) and 1 column (label)
embedding_df = pd.DataFrame(adata.obsm['X_scanit'], index=adata.obs_names)  # DataFrame with index (cell-id/barcode) and n columns

## Write output
out_dir.mkdir(parents=True, exist_ok=True)

label_df.columns = ["label"]
label_df.to_csv(label_file, sep="\t", index_label="")

if embedding_df is not None:
    embedding_df.to_csv(embedding_file, sep="\t", index_label="")
