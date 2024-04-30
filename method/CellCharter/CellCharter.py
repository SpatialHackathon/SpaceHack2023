#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Jieran Sun; implemented method cellCharter

import argparse

# TODO adjust description
parser = argparse.ArgumentParser(description="Method CellCharter: https://doi.org/10.1038/s41588-023-01588-4")

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

################# Get args #######################

args = parser.parse_args()

n_clusters = args.n_clusters
technology = args.technology
seed = args.seed
# Def config
import json
with open(args.config, "r") as f:
    config = json.load(f)

# Set up output files
from pathlib import Path
out_dir = Path(args.out_dir)
label_file = out_dir / "domains.tsv"
embedding_file = out_dir / "embedding.tsv"

################# Define functions #################
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
        .iloc[:,0:2]
        .loc[observations.index, :]
        .to_numpy()
    )

    adata = ad.AnnData(obs=observations, 
                       var=features, 
                       obsm={"spatial": coordinates})

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

## Your code goes here
import cellcharter as cc
import anndata as ad
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import random
import torch

random.seed(seed)

if args.n_genes is not None:
    n_genes = args.n_genes
else:
    n_genes = config.get("n_genes", None)  # No filtering if no top genes are specified
    
# Load data
adata = get_anndata(args)
adata.var_names_make_unique()

# raw count input
import scvi

# Assume filtered or not?
# sc.pp.filter_genes(adata, min_counts=3)
#sc.pp.filter_cells(adata, min_counts=3)

#GPU possibility
use_cuda = torch.cuda.is_available()

# Follow https://github.com/CSOgroup/cellcharter_analyses/blob/main/src/benchmarking/CellCharter/individual.py
sc.pp.filter_genes(adata, min_counts=config["min_counts"])

# preprocessing for scvi
adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=n_genes,
    subset=True,
    layer="counts",
    flavor="seurat_v3",

)

# Set up scvi model for reduced_dimension embedding
scvi.settings.seed = seed
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts"
)
model = scvi.model.SCVI(adata, n_latent=config["n_latent"])
model.train(early_stopping=True, enable_progress_bar=False, progress_bar_refresh_rate=0)
adata.obsm['reduced_dimensions'] = model.get_latent_representation(adata).astype(np.float32)
# Find neighbors
coord = "grid" if technology in ["Visium", "ST"] else "generic"
delTri = technology not in ["Visium", "ST"] 

sq.gr.spatial_neighbors(adata, 
                        coord_type = coord, 
                        delaunay=delTri)

# Trim if only delaunay is used
if delTri:
    cc.gr.remove_long_links(adata)

# Aggregate neighborhoods
cc.gr.aggregate_neighbors(adata, 
                          n_layers=config["n_layers"], 
                          aggregations=config["aggregations"],
                          use_rep='reduced_dimensions', 
                          out_key='X_cellcharter', 
                          sample_key=None)

# clustering
model = cc.tl.Cluster(
    n_clusters=n_clusters, 
    random_state=seed,
    convergence_tolerance=config["convergence_tolerance"], 
    covariance_regularization=config["covariance_regularization"],
    # if running on gpu
    trainer_params=dict(accelerator='gpu', devices=1) if use_cuda else None
    )

model.fit(adata, use_rep='X_cellcharter')

# label prediction
label_df = pd.DataFrame({"label" : model.predict(adata, use_rep='X_cellcharter')}) 
label_df.index = adata.obs.index

embedding_df = pd.DataFrame(adata.obsm["X_cellcharter"])
embedding_df.index = adata.obs_names

## Write output
out_dir.mkdir(parents=True, exist_ok=True)

label_df.to_csv(label_file, sep="\t", index_label="")
embedding_df.to_csv(embedding_file, sep="\t", index_label="")
