#!/usr/bin/env python

# Author_and_contribution: Peiying Cai; created script

import argparse

parser = argparse.ArgumentParser(
    description="""stLearn (https://www.biorxiv.org/content/10.1101/2020.05.31.125658v1)
Images can be used."""
)

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
    "--n_clusters", help="Number of clusters to return.", required=False, type=int
)

parser.add_argument(
    "--technology",
    help="The technology of the dataset (Visium, ST, ...).",
    required=False,
) # Required only when the path to the H&E image is defined.

parser.add_argument(
    "--seed", help="Seed to use for random operations.", required=True, type=int
)
parser.add_argument(
    "--config",
    help="Optional config file used to pass additional parameters.",
    required=False,
)

args = parser.parse_args()

from pathlib import Path

out_dir = Path(args.out_dir)

# Output files
label_file = out_dir / "domains.tsv"
embedding_file = out_dir / "embedding.tsv"
embedding_df = None

# If the value of n_clusters is undefned, use the stLearn default setting (n_clusters = 20).
if args.n_clusters is not None:
    n_clusters = args.n_clusters
else:
    n_clusters = 20

if args.image is not None:
    import os
    if args.technology is None:
        raise Exception(
            f"Invalid argument. The technology of the dataset should be defined when the path of H&E staining is defined."
        )
    if args.technology not in ["Visium", "Old_ST"]:
        raise Exception(
            f"Invalid {args.technology}. The technology of the dataset should be 'Visium' or 'Old_ST' not {args.technology}."
        )
    technology = args.technology
    image_dir = args.image
    positions_dir = os.path.join(image_dir, 'tissue_positions_list.csv')
    scale_dir = os.path.join(image_dir, 'scalefactors_json.json')
    image_lowers_dir = os.path.join(image_dir, 'tissue_lowres_image.png')


seed = args.seed


## Your code goes here
import json

with open(args.config, "r") as f:
    config = json.load(f)

import random
import numpy as np
import pandas as pd
import stlearn as st
from pathlib import Path
import torch


def get_anndata(args):
    import anndata as ad
    import scipy as sp
    # from PIL import Image

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

    adata.obs["imagecol"] = adata.obs["col"]
    adata.obs["imagerow"] = adata.obs["row"]

    if args.image is not None:
        if technology == 'Visium':
            adata.obs["array_row"] = adata.obsm["spatial_pixel"][:, 0]
            adata.obs["array_col"] = adata.obsm["spatial_pixel"][:, 1]
        # Add the image
        st.add.positions(adata=adata, position_filepath=positions_dir, 
                         scale_filepath=scale_dir, quality='low')
        st.add.image(adata=adata, imgpath=image_lowers_dir, 
                     library_id='sample_name', quality='lowres', visium=False)
        # adata.uns["spatial"] = np.array(Image.open(args.image))
        # Output of pre-processing image
        tile_dir = Path(out_dir / "tiles")
        tile_dir.mkdir(parents=True, exist_ok=True)
        # Pre-processing for spot image
        st.pp.tiling(adata, tile_dir)
        # This step uses deep learning model to extract high-level features from tile images
        # May need few minutes to be completed
        st.pp.extract_feature(adata)
    else:
        adata.uns["spatial"] = None

    return adata


adata = get_anndata(args)


# Set seed
random.seed(seed)
torch.manual_seed(seed)
np.random.seed(seed)

if adata.uns["spatial"] is not None:
    st.em.run_pca(adata, n_comps=config["n_comps"])
    # Apply stSME to normalise log transformed data
    # With weights from morphological Similarly and physcial distance
    # TODO whether put "weights": "weights_matrix_pd_md" in `config``
    st.spatial.SME.SME_normalize(adata, use_data="raw",
                                 weights="weights_matrix_pd_md",
                                 platform = technology
                                )
    adata.X = adata.obsm['raw_SME_normalized']
    embedding_df = adata.obsm["X_morphology"]

# TODO whether provide highly variable genes
# by default `st.em.run_pca` uses highly variable genes if they have been stored in `.var['highly_variable']` beforehand

# If a reduced dimensionality representation is given, use it; otherwise, calculate and use PCA.
if args.dim_red is not None:
    adata.obsm["dim_red"] = args.dim_red
    use_data = "dim_red"
else:
    use_data = "X_pca"

st.pp.scale(adata)
st.em.run_pca(adata, n_comps=config["n_comps"], random_state=seed)
# Run
if config["method"] == "kmeans":
    # K-means clustering on stSME normalised PCA
    st.tl.clustering.kmeans(adata,n_clusters=n_clusters, use_data=use_data, key_added="cluster", random_state=seed)
    st.pl.cluster_plot(adata, use_label="cluster")
elif config["method"] == "louvain":
    st.pp.neighbors(adata,n_neighbors=config["n_neighbors"],use_rep=use_data, random_state=seed)
    st.tl.clustering.louvain(adata, resolution=config["res"], key_added="cluster", random_state=seed)

label_df = adata.obs[["cluster"]]
## Write output
out_dir.mkdir(parents=True, exist_ok=True)

label_df.columns = ["label"]
label_df.to_csv(label_file, sep="\t", index_label="")

if embedding_df is not None:
    embedding_df.to_csv(embedding_file, sep="\t", index_label="")
