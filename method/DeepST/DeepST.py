#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Jieran Sun; Implemented method DeepST

import argparse

# description
parser = argparse.ArgumentParser(description="Method DeepST: identifying spatial domains in spatial transcriptomics by deep learning https://doi.org/10.1093/nar/gkac901")

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
    # Def config
    import json
    with open(args.config, "r") as f:
        config = json.load(f)

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

    if args.image is not None:
        adata.uns["image"] = np.array(Image.open(args.image))
    else:
        adata.uns["image"] = None

    return adata

adata = get_anndata(args)
adata.var_names_make_unique()

#######
# NOTE: Seed cannot be used here

# First, set up DeepST environment
import tempfile
import os, sys
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch

# Seeding for controlling randomization
import random
random.seed(seed)
torch.manual_seed(seed)
np.random.seed(seed)

# Work in a temprary folder
with tempfile.TemporaryDirectory() as tmpdir:
    gitdir = f"{str(tmpdir)}/DeepST"

    # Clone the repository to the specific commit
    os.system(
    f"""
    git clone https://github.com/JiangBioLab/DeepST.git {gitdir}
    cd {gitdir} 
    git reset --hard 1daa513
    """)

    # Set working directory as deepST directory
    sys.path.append(f"{gitdir}/deepst")

    # From DeepST import run
    import importlib.util
    # Import the main wrapper DeepST module
    spec=importlib.util.spec_from_file_location("deepST", f"{gitdir}/deepst/DeepST.py")
    deepST = importlib.util.module_from_spec(spec)
    sys.modules["deepST"] = deepST
    spec.loader.exec_module(deepST)

    use_gpu = True if torch.cuda.is_available() else False
    # Set up deepST object
    deepen = deepST.run(save_path = None,
                 task = "Identify_Domain", 
                 pre_epochs = 800,  # change based on hardware
                 epochs = 1000,  # change based on hardware
                 use_gpu = use_gpu,
                 )

    # Adopted from github tutorial
    use_morphological = args.image is not None
    if use_morphological:
        adata = deepen._get_image_crop(adata, data_name=technology) 

    # Data augmentation. spatial_type includes three kinds of "KDTree", "BallTree" and "LinearRegress", among which "LinearRegress" is only applicable to 10x visium and the remaining omics selects the other two. "use_morphological" defines whether to use morphological images.
    spatial_type = "BallTree" if technology != "Visium" and config["spatial_type"]=="LinearRegress" else config["spatial_type"]

    adata = deepen._get_augment(adata, spatial_type=spatial_type, use_morphological=use_morphological)

    graph_dict = deepen._get_graph(adata.obsm["spatial"], distType = "BallTree")

    data = deepen._data_process(adata, pca_n_comps = min(config["npcs"], adata.n_vars))

    """ Original SpaceHack implementation
    # adopted code from deepst/adj.py
    " Store original adjacency matrix (without diagonal entries) for later "
    adj_pre = adata.obsp["spatial_connectivities"]
    adj_pre = sp.coo_matrix(adj_pre)

    adj_pre = adj_pre - sp.dia_matrix((adj_pre.diagonal()[np.newaxis, :], [0]), shape=adj_pre.shape)
    adj_pre.eliminate_zeros()

    " Some preprocessing."
    # intiate a "graph" item, k is a placeholder doesn't matter
    graph = adj.graph(data = adata, k = 1)
    adj_norm = graph.pre_graph(adj_pre)
    adj_label = adj_pre + sp.eye(adj_pre.shape[0])
    adj_label = torch.FloatTensor(adj_label.toarray())
    norm = adj_pre.shape[0] * adj_pre.shape[0] / float((adj_pre.shape[0] * adj_pre.shape[0] - adj_pre.sum()) * 2)

    graph_dict = {
                 "adj_norm": adj_norm,
                 "adj_label": adj_label,
                 "norm_value": norm }

    data = adata.obsm["reduced_dimensions"].astype(np.float64)
    """

    deepst_embed = deepen._fit(data, graph_dict=graph_dict)

    adata.obsm["DeepST_embed"] = deepst_embed
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
