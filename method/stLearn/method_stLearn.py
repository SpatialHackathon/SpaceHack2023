#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Peiying Cai, implemented method. 
import argparse

parser = argparse.ArgumentParser(
    description="""stLearn (https://www.biorxiv.org/content/10.1101/2020.05.31.125658v1; https://www.nature.com/articles/s41467-023-43120-6)
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
parser.add_argument("--image", nargs='+', help="Path to H&E staining.", required=False)
parser.add_argument(
    "--n_clusters", help="Number of clusters to return.", required=True, type=int
)

parser.add_argument(
    "--technology",
    help="The technology of the dataset (Visium, ST, ...).",
    required=True,
) # Actually required only when the path to the H&E image is defined.

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

n_clusters = args.n_clusters
import os
import json
with open(args.config, "r") as f:
    config = json.load(f)

# Map technology spelling
def map_technology(input_technology):
    technology_mapping = {
        "Visium": "Visium",
        "ST": "Old_ST"
    }
    return technology_mapping.get(input_technology, None)
    
if args.image is not None and config["normalize"] is True:
    technology = map_technology(str(args.technology))
    if technology is None:
        raise Exception(
            f"Invalid technology. stLearn with stSME clustering only supports the spatial transcriptomics platform that includes H&E images, like 10X Visium and ST not {args.technology}. "
            )

# Return warning if no reduced dimensionality representation provided
if args.dim_red is None:
    raise Exception(
            f"No dimensionality reduction found. Provide path to representation (e.g., PCA). "
            )

seed = args.seed
## Your code goes here
import random
import numpy as np
import pandas as pd
import stlearn as st

def get_anndata(args):
    import anndata as ad
    import scipy as sp
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

    adata.obs["imagecol"] = adata.obsm["spatial_pixel"][:, 0]
    adata.obs["imagerow"] = adata.obsm["spatial_pixel"][:, 1]
    adata.obs["array_row"] = adata.obs["row"]
    adata.obs["array_col"] = adata.obs["col"]
    # Use morphology information
    # Modified from stLearn.adds.add_image.py
    if args.image is not None:
        image_tif = args.image[0]
        image_json = args.image[1] if len(args.image) == 2 else None
        from PIL import Image
        # Add image
        Image.MAX_IMAGE_PIXELS=None
        im = Image.open(image_tif)
        image = np.array(im)
        spatial_key = "spatial"
        library_id = "sample"
        quality = "hires"
        try:
            adata.uns[spatial_key] = {library_id: {}}
            adata.uns[spatial_key][library_id]["images"] = {}
            adata.uns[spatial_key][library_id]["images"][quality] = image
            adata.uns[spatial_key][library_id]["use_quality"] = quality

            if technology != 'Visium':
                adata.uns[spatial_key][library_id]["scalefactors"] = {}
                if image_json is None:
                    raise ValueError("Please provide the path to H&E staining (as json).")
                else:
                    with open(image_json) as json_file:
                        scale_factors = json.load(json_file)
                    scale = scale_factors.get('scale')
                    spot_diameter_fullres = scale_factors.get('spot_dimeter_fullres')
                    if scale is None:
                        raise ValueError("Please provide 'scale' in the path to H&E staining (as json).")
                    if spot_diameter_fullres is None:
                        raise ValueError("Please provide 'spot_diameter_fullres' in the path to H&E staining (as json).")

                    adata.uns[spatial_key][library_id]["scalefactors"]["tissue_hires_scalef"] = scale
                    adata.uns[spatial_key][library_id]["scalefactors"]["spot_diameter_fullres"] = spot_dimeter_fullres
                    adata.obsm[spatial_key] = adata.obs[["imagecol", "imagerow"]].values
                    adata.obs[["imagecol", "imagerow"]] = adata.obsm["spatial"] * scale

            print("Added tissue image to the object.")
        except:
            raise ValueError(
                f"""\
                {args.image!r} does not end on a valid extension.
                """
                )
        # Pre-processe image
        tile_dir = Path(out_dir / "tiles")
        tile_dir.mkdir(parents=True, exist_ok=True)
        # Pre-processing for spot image
        st.pp.tiling(adata, tile_dir)
        # This step uses deep learning model to extract high-level features from tile images
        # May need few minutes to be completed
        st.pp.extract_feature(adata)
        # Convert the NumPy array to Pandas DataFrame with row names
        embedding_df = pd.DataFrame(adata.obsm['X_morphology'], index=adata.obs_names)
    else:
        adata.uns["spatial"] = None

    return adata


adata = get_anndata(args)


# Set seed
random.seed(seed)
np.random.seed(seed)

# Use the provided reduced dimensionality representation.
adata.obsm["X_pca"] = pd.read_table(args.dim_red, index_col=0).loc[adata.obs_names].to_numpy()

use_data = "X_pca"
if config["image"] is True and adata.uns["spatial"] is not None:
    if config["normalize"]:
        # Apply stSME to normalise log transformed data
        # With weights from morphological Similarly and physcial distance
        st.spatial.SME.SME_normalize(adata, use_data="raw",
                                     weights=config["weights"],
                                     platform = technology
                                    )
        adata.X = adata.obsm['raw_SME_normalized']
        st.pp.scale(adata)
        st.em.run_pca(adata, n_comps=config["n_comps"])
    else:
        # Adjust the PCA/UMAP expression space based on identified neighbors and morphological distance by using the disk smoothing approach 
        st.spatial.smooth.disk(adata, use_data=use_data)
        use_data = use_data + "_disk"

def res_search_fixed_clus_louvain(adata, n_clusters, increment=0.01, random_seed=2023):
    for res in np.arange(0.2, 2, increment):
        st.tl.clustering.louvain(adata, resolution=res, key_added="res_search", random_state=seed)
        if len(adata.obs['res_search'].unique()) > n_clusters:
            break
    return res-increment
    
# Run
if config["method"] == "kmeans":
    st.tl.clustering.kmeans(adata,n_clusters=n_clusters, use_data=use_data, key_added="cluster", random_state=seed)
elif config["method"] == "louvain":
    st.pp.neighbors(adata,n_neighbors=config["n_neighbors"],use_rep=use_data, random_state=seed)
    res = res_search_fixed_clus_louvain(adata,n_clusters=n_clusters, increment=0.01, random_seed=seed)
    st.tl.clustering.louvain(adata, resolution=res, key_added="cluster", random_state=seed)

label_df = adata.obs[["cluster"]]
## Write output
out_dir.mkdir(parents=True, exist_ok=True)

label_df.columns = ["label"]
label_df.to_csv(label_file, sep="\t", index_label="")

if embedding_df is not None:
    embedding_df.to_csv(embedding_file, sep="\t", index_label="")
