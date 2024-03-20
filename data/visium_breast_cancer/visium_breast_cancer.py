#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Kirti Biharie; added dataset

import argparse

parser = argparse.ArgumentParser(description="Load data for Visium Breast Cancer")

parser.add_argument(
    "-o", "--out_dir", help="Output directory to write files to.", required=True
)

args = parser.parse_args()


from pathlib import Path

import pandas as pd

out_dir = Path(args.out_dir)

# The folder structure should look like the following
# out_dir
# |_______sample_1  (sample name can be chosen freely)
#         |_____coordinates.tsv
#         |_____features.tsv
#         |_____observations.tsv
#         |_____counts.mtx  (use scipy.io.mmwrite)
#         |_____labels.tsv  (optional)
#         |_____H_E.(tiff/png/...)  (optional)
#         |_____H_E.json  (optional, required if H_E is provided)
# |_______sample_2
#         | ...
# |_______samples.tsv
# |_______experiment.json
# if additional output files are required write it also to out_dir


## Your code goes here
import os
import tempfile
import urllib.request
import urllib.parse
import scanpy as sc
import tarfile
import json
import numpy as np
import shutil

sample_columns = ["patient","sample","position","replicate","n_clusters","directory"]
samples_df = pd.DataFrame(columns=sample_columns)

samples_info = [{
                "sample_id": "CytAssist_FFPE_Human_Breast_Cancer", 
                "filtered_matrix": "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_filtered_feature_bc_matrix.h5",
                "spatial": "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Human_Breast_Cancer/CytAssist_FFPE_Human_Breast_Cancer_spatial.tar.gz",
                "labels": "https://cdn.10xgenomics.com/raw/upload/v1695234604/Xenium%20Preview%20Data/Cell_Barcode_Type_Matrices.xlsx"
               }]

def download_link(link, filename, temp_dir):
    try:
        request = urllib.request.Request(url=link, headers={'User-Agent': 'Mozilla/6.0'})
        response = urllib.request.urlopen(request)
        # Extract filename from the URL
        filename = os.path.join(temp_dir, filename)
        with open(filename, 'wb') as file:
            file.write(response.read())
        print(f"Downloaded: {filename}")
    except Exception as e:
        print(f"Error downloading {link}: {e}")

# Example how a sample could be written
def write_sample(
    path,
    sample,
    coordinates_df,
    observations_df,
    features_df,
    counts,
    labels_df=None,
    img_path=None,
    scalefactor=None
):


    import scipy as sp

    sample_path = Path(path) / sample
    os.makedirs(sample_path, exist_ok=True)

    coordinates_df.to_csv(sample_path / "coordinates.tsv", sep="\t", index_label="")
    features_df.to_csv(sample_path / "features.tsv", sep="\t", index_label="")
    observations_df.to_csv(sample_path / "observations.tsv", sep="\t", index_label="")
    sp.io.mmwrite(sample_path / "counts.mtx", counts)

    if img_path is not None:
        shutil.copy(img_path, os.path.join(sample_path, "H_E.png"))

    if scalefactor is not None:
        with open(os.path.join(sample_path, "H_E.json"), "w") as f:
            json.dump(scalefactor, f)

    if labels_df is not None:
        labels_df.columns = ["label"]
        labels_df.to_csv(sample_path / "labels.tsv", sep="\t", index_label="")

with tempfile.TemporaryDirectory() as temp_dir: 
    for i, sample_info in enumerate(samples_info):
        sample_id = sample_info["sample_id"]
        sample_temp_dir = os.path.join(temp_dir, sample_id)
        os.makedirs(sample_temp_dir)

        download_link(sample_info["filtered_matrix"], "filtered_feature_bc_matrix.h5", sample_temp_dir)
        download_link(sample_info["spatial"], "spatial.tar.gz", sample_temp_dir)
        download_link(sample_info["labels"], "labels.xlsx", sample_temp_dir)

        TarFile = tarfile.open(os.path.join(sample_temp_dir, "spatial.tar.gz"),"r:gz")
        TarFile.extractall(sample_temp_dir)
        os.rename(os.path.join(sample_temp_dir, "spatial/tissue_positions.csv"), os.path.join(sample_temp_dir, "spatial/tissue_positions_list.csv"))
        
        adata = sc.read_visium(sample_temp_dir)

        all_labels_df = pd.read_excel(os.path.join(sample_temp_dir, "labels.xlsx"), sheet_name="Visium", index_col=0)
        all_labels_df = all_labels_df.loc[adata.obs_names]
        all_labels_df["label"] = all_labels_df["Cluster"].astype(str) + "-" + all_labels_df["Annotation"]
        labels_df = all_labels_df[["label"]]

        coordinates_df = pd.DataFrame(adata.obsm["spatial"], index=adata.obs_names, columns=["x","y"])

        counts = adata.X

        features_df = adata.var.copy()
        features_df["selected"] = True

        observations_df = adata.obs.copy()
        observations_df.columns = ["in_tissue", "row", "col"]
        observations_df["selected"] = True

        img_path = os.path.join(sample_temp_dir, "spatial", "tissue_hires_image.png")
        scalefactors_path = os.path.join(sample_temp_dir, "spatial", "scalefactors_json.json")
        
        with open(scalefactors_path, "r") as f:
            scalefactor = json.load(f)["tissue_hires_scalef"]
            scalefactor = {"scale": scalefactor}

        samples_df.loc[i] = [i, sample_id, np.nan, np.nan, labels_df["label"].nunique(), sample_id]
        
        write_sample(out_dir,
                    sample_id,
                    coordinates_df,
                    observations_df,
                    features_df,
                    counts,
                    labels_df,
                    img_path=img_path,
                    scalefactor=scalefactor)

technology = "Visium"

## Metadata files
samples_df.loc[
    :, ["patient", "sample", "position", "replicate", "directory", "n_clusters"]
].to_csv(out_dir / "samples.tsv", sep="\t", index_label="")

import json

with open(out_dir / "experiment.json", "w") as f:
    exp_info = {"technology": technology}
    json.dump(exp_info, f)
