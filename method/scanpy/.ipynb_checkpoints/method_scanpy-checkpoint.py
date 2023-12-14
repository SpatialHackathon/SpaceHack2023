#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Liya Zaygerman

import argparse

# TODO adjust description
parser = argparse.ArgumentParser(description="Method scanpy")

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
# if additional output files are required write it also to out_dir

# Use these filepaths as input ...
coord_file = args.coordinates
feature_file = args.features
observation_file = args.observations

if args.neighbors is not None:
    neighbors_file = args.neighbors
if args.matrix is not None:
    matrix_file = args.matrix
if args.dim_red is not None:
    dimred_file = args.dim_red
if args.image is not None:
    image_file = args.image
if args.config is not None:
    config_file = args.config

n_clusters = args.n_clusters
technology = args.technology
seed = args.seed

import anndata as ad
import numpy as np
import pandas as pd
import scipy as sp
from PIL import Image

# ... or AnnData if you want
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
        adata.uns["image"] = np.array(Image.open(args.image))
    else:
        adata.uns["image"] = None

    return adata


adata = get_anndata(args)
adata.var_names_make_unique()


# TODO set the seed, if the method requires the seed elsewhere please pass it on
import random

random.seed(seed)
# np.random.seed(seed)
# torch.manual_seed(seed)

## Your code goes here
import scanpy as sc
import json

# get the json config
with open (args.config, "r") as c:
    config = json.load(c)

# scanpy settings
# sc.settings.set_figure_params(scanpy=True, dpi=80, dpi_save=500, frameon=True, vector_friendly=True, fontsize=14)

if __name__ == "__main__":
    # throw warning about not using the num_clusters as a parameter, because scanpy uses leiden or louvain and needs the resolution parameter, which is defined in the config.json. There exists a good function to perform an extensive search for the right resolution parameter to define the desired num_clusters, but we still have licensing issues. For more, see the SpaceHack2.0 GitHub issue #139
    import warnings
    
    warnings.warn("Scanpy uses leiden/louvain for clustering, which relies on the resolution parameter in the config file. The parameter num_clusters will be ignored.", UserWarning)
    
    # scanpy starts here
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_neighbors = config["n_neighbors"], n_pcs=config["n_pcs"], random_state=seed)

    #two options - leiden or loivain
    if config['clustering'] == "louvain":
        sc.tl.louvain(adata, resolution=config["resolution"], random_state=seed, key_added='louvain')
    elif config['clustering'] == "leiden":
        sc.tl.leiden(adata, resolution=config["resolution"], random_state=seed)
    else: 
        print("No clustering method defined or your method is not available, performing leiden")
        sc.tl.leiden(adata, resolution=config["resolution"], random_state=seed)
    
    label_df = adata.obs["leiden"]


# embedding_df = None  # optional, DataFrame with index (cell-id/barcode) and n columns


## Write output
out_dir.mkdir(parents=True, exist_ok=True)

label_df.columns = ["label"]
label_df.to_csv(label_file, sep="\t", index_label="")

# if embedding_df is not None:
#     embedding_df.to_csv(embedding_file, sep="\t", index_label="")
