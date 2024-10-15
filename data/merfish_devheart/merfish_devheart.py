#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Kirti Biharie; added dataset

import argparse

parser = argparse.ArgumentParser(description="Load data for Merfish Developing Heart")

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

links = [
         {"name": "ventricle_15pcw_merfish.h5ad", "link":"https://datadryad.org/stash/downloads/file_stream/2801011"}, 
         {"name": "overall_merfish.h5ad", "link": "https://datadryad.org/stash/downloads/file_stream/2801003"},
         {"name": "R123_N3S5_single_cell_raw_counts.csv", "link": "https://datadryad.org/stash/downloads/file_stream/2801014"},
         {"name": "R77_4C4_single_cell_raw_counts.csv", "link": "https://datadryad.org/stash/downloads/file_stream/2800996"},
         {"name": "R78_4C12_single_cell_raw_counts.csv", "link": "https://datadryad.org/stash/downloads/file_stream/2800997"},
         {"name": "R78_4C15_single_cell_raw_counts.csv", "link": "https://datadryad.org/stash/downloads/file_stream/2800999"},
        ]

def download_link(link, filename, temp_dir):
    try:
        request = urllib.request.Request(url=link, headers={'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_0) AppleWebKit/600.1.17 (KHTML, like Gecko) Version/8.0 Safari/600.1.17'})
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

    if type(counts) == pd.DataFrame:
        counts = counts.to_numpy()
    counts = sp.sparse.csr_matrix(counts)
    sp.io.mmwrite(sample_path / "counts.mtx", counts)

    if img_path is not None:
        shutil.copy(img_path, os.path.join(sample_path, "H_E.png"))

    if scalefactor is not None:
        with open(os.path.join(sample_path, "H_E.json"), "w") as f:
            json.dump(scalefactor, f)

    if labels_df is not None:
        # labels_df.columns = ["label"]
        labels_df.to_csv(sample_path / "labels.tsv", sep="\t", index_label="")

with tempfile.TemporaryDirectory() as temp_dir:
    for link in links:
        download_link(link["link"], link["name"], temp_dir)

    heart = sc.read_h5ad(os.path.join(temp_dir, "overall_merfish.h5ad"))
    ventricle = sc.read_h5ad(os.path.join(temp_dir, "ventricle_15pcw_merfish.h5ad"))

    # Process 3 13PCW heart samples
    for sample_name in ["R77_4C4", "R78_4C12", "R78_4C15"]:
        heart_sample = heart[heart.obs["sample_id"] == sample_name]
        
        counts = pd.read_csv(os.path.join(temp_dir, f"{sample_name}_single_cell_raw_counts.csv"), index_col=0)
        counts.index = counts.index.astype(str) + "-" + sample_name
        counts = counts.loc[heart_sample.obs_names]
        
        coordinates_df = pd.DataFrame(heart_sample.obsm["spatial"], index=heart_sample.obs_names, columns=["x", "y"])
        labels_df = heart_sample.obs[["communities", "populations"]]
        labels_df.columns = ["label", "populations"]
        features_df = pd.DataFrame(index=heart_sample.var_names)
        features_df["selected"] = True
        observations_df = pd.DataFrame(index=heart_sample.obs_names)
        observations_df["selected"] = True
        samples_df.loc[len(samples_df)] = [0, sample_name, np.nan, np.nan, len(labels_df["label"].cat.categories), sample_name]

        write_sample(out_dir,
                    sample_name,
                    coordinates_df,
                    observations_df,
                    features_df,
                    counts,
                    labels_df)

    # Process 15PCW ventricle sample
    sample_name = "R123_N3S5"

    counts = pd.read_csv(os.path.join(temp_dir, f"{sample_name}_single_cell_raw_counts.csv"), index_col=0)
    counts.index = counts.index.astype(str)
    counts = counts.loc[ventricle.obs_names]
    
    coordinates_df = pd.DataFrame(ventricle.obsm["X_spatial"], index=ventricle.obs_names, columns=["x", "y"])
    labels_df = ventricle.obs[["subpopulations"]]
    labels_df.columns = ["label"]
    features_df = pd.DataFrame(index=ventricle.var_names)
    features_df["selected"] = True
    observations_df = pd.DataFrame(index=ventricle.obs_names)
    observations_df["selected"] = True
    samples_df.loc[len(samples_df)] = [1, sample_name, np.nan, np.nan, len(labels_df["label"].cat.categories), sample_name]

    write_sample(out_dir,
                    sample_name,
                    coordinates_df,
                    observations_df,
                    features_df,
                    counts,
                    labels_df)

technology = "Merfish"

## Metadata files
samples_df.loc[
    :, ["patient", "sample", "position", "replicate", "directory", "n_clusters"]
].to_csv(out_dir / "samples.tsv", sep="\t", index_label="")

import json

with open(out_dir / "experiment.json", "w") as f:
    exp_info = {"technology": technology}
    json.dump(exp_info, f)
