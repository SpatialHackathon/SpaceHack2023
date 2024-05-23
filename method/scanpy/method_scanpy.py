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
parser.add_argument(
    "--n_pcs",
    help="Optional specific number of prinicpal components to use.",
    required=False,
)
parser.add_argument(
    "--n_genes",
    help="Optional specific number of variable features to take.",
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

# import res-n_clust tuning function
import sys
from pathlib import Path

# Add the parent directory of the current file to sys.path
method_dir = Path(__file__).resolve().parent.parent  # Navigate two levels up
sys.path.append(str(method_dir))

from search_res import binary_search

# get the json config
with open (args.config, "r") as c:
    config = json.load(c)

# scanpy settings
# sc.settings.set_figure_params(scanpy=True, dpi=80, dpi_save=500, frameon=True, vector_friendly=True, fontsize=14)


    # throw warning about not using the num_clusters as a parameter, because scanpy uses leiden or louvain and needs the resolution parameter, which is defined in the config.json. There exists a good function to perform an extensive search for the right resolution parameter to define the desired num_clusters, but we still have licensing issues. For more, see the SpaceHack2.0 GitHub issue #139
import warnings

# warnings.warn("Scanpy uses leiden/louvain for clustering, which relies on the resolution parameter in the config file. The parameter num_clusters will be ignored.", UserWarning)

if "n_genes" not in config.keys():
    config["n_genes"] = None
if "n_pcs" not in config.keys():
    config["n_pcs"] = None

n_pcs = args.n_pcs if args.n_pcs is not None else config["n_pcs"]
n_genes = args.n_genes if args.n_genes is not None else config["n_genes"]
directed = None if "directed" not in config.keys() else config["directed"]
n_iterations = -1 if "n_iterations" not in config.keys() else config["n_iterations"]
flavor = "leidenalg" if "flavor" not in config.keys() else config["flavor"]

# scanpy starts here
if not args.dim_red:
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    if adata.n_vars > 2000:
        sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=n_genes)
    sc.pp.scale(adata)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_neighbors=config["n_neighbors"], n_pcs=n_pcs, random_state=seed)
else:
    sc.pp.neighbors(adata, n_neighbors=config["n_neighbors"], use_rep="reduced_dimensions")

#two options - leiden or loivain
if config['clustering'] not in ["louvain", "leiden"]:
    print("No clustering method defined or your method is not available, performing leiden")
    label_df = binary_search(adata, n_clust_target=n_clusters, method="leiden", seed = seed)
else:
    label_df = binary_search(adata, n_clust_target=n_clusters, 
                             method=config['clustering'], 
                             seed = seed, 
                             directed = directed,
                             n_iterations = n_iterations,
                             flavor = flavor,)

# sc.tl.leiden(adata, resolution=config["resolution"], random_state=seed)

# label_df = adata.obs[["leiden"]]


# embedding_df = None  # optional, DataFrame with index (cell-id/barcode) and n columns


## Write output
out_dir.mkdir(parents=True, exist_ok=True)

label_df.columns = ["label"]
label_df.to_csv(label_file, sep="\t", index_label="")

# if embedding_df is not None:
#     embedding_df.to_csv(embedding_file, sep="\t", index_label="")
