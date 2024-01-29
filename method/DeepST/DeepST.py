#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: ENTER YOUR NAME AND CONTRIBUTION HERE

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

    if args.img is not None:
        adata.uns["image"] = np.array(Image.open(args.img))
    else:
        adata.uns["image"] = None

    return adata

adata = get_anndata(args)

#######

# TODO set the seed, if the method requires the seed elsewhere please pass it on
import random

random.seed(seed)
# np.random.seed(seed)
# torch.manual_seed(seed)

# First, set up DeepST environment
import tempfile
import os
from pathlib import Path
import scanpy as sc

# Work in a temprary folder
with tempfile.TemporaryDirectory() as tmpdir:
    gitdir = Path(tmpdir) / "DeepST"
        
    print(f"Created temporary directory at {tmpdir} to store DeepST repo")

    # Clone the repository to the specific commit
    os.system(
    f"""
    git clone https://github.com/JiangBioLab/DeepST.git {gitdir}
    cd {gitdir} 
    git reset --hard 1daa513
    """)
    
    # From DeepST import run
    import importlib.util
    spec=importlib.util.spec_from_file_location("run", gitdir/"deepst/DeepST.py")
    run = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(run)
    
    # from his_feat import image_feature, image_crop

    # Set up deepST object
    deepen = run(save_path=None,
                 task = "Identify_Domain", 
                 platform=technology,
                 pca_n_comps=200,
                 pre_epochs=800,  # change based on hardware
                 epochs=1000,  # change based on hardware
                 Conv_type="GCNConv",  # ["GCNConv", ]
                 )
    

    # TODO: get graph_dict structure from input
    graph_dict = deepen._get_graph(
        adata.obsm["spatial"], distType="BallTree", k=k)

    # TODO: Parameters changed by configurations
    adata = deepen._fit(adata, graph_dict, pretrain=True)

    # Get clustering results, set prior=True so we got n_clusters no of clusters
    adata = deepen._get_cluster_data(adata, n_domains=n_clusters, priori=True)
    
    # Output dataframes
    label_df = adata.obs[["DeepST_refine_domain"]]
    embedding_df = pd.DataFrame(adata.obsm['DeepST_embed'])
    embedding_df.index = adata.obs_names
    
    ## Write output
    out_dir.mkdir(parents=True, exist_ok=True)
    label_df.columns = ["label"]
    label_df.to_csv(label_file, sep="\t", index_label="")
    
    if embedding_df is not None:
        embedding_df.to_csv(embedding_file, sep="\t", index_label="")
