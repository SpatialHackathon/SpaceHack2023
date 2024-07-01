#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Kirti Biharie; added dataset

import argparse

parser = argparse.ArgumentParser(description="Load data for cosmx lung dataset")

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
import json
import numpy as np
import shutil
import anndata as ad
import subprocess
import rpy2.robjects
import rpy2.robjects.pandas2ri
import tarfile

samples = ['Lung5_Rep1', 'Lung5_Rep2', 'Lung5_Rep3', 'Lung6', 'Lung9_Rep1', 'Lung9_Rep2', 'Lung12', 'Lung13']
patients = [0,0,0,1,2,2,3,4]
replicates = ["Rep1", "Rep2", "Rep3", np.nan, "Rep1", "Rep2", np.nan, np.nan]

sample_columns = ["patient","sample","position","replicate","n_clusters","directory"]
samples_df = pd.DataFrame(columns=sample_columns)

links = [
         {"name": "All+SMI+Flat+data.tar.gz", "link":"https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI-Compressed/All+SMI+Flat+data.tar.gz"}, 
         {"name": "All+SMI+Giotto+object.tar.gz", "link":"https://nanostring-public-share.s3.us-west-2.amazonaws.com/SMI-Compressed/All+SMI+Giotto+object.tar.gz"}, 
        ]

def download_link(link, filename, temp_dir):
    try:
        filename = os.path.join(temp_dir, filename)
        subprocess.run(["wget", "-O", filename,link])
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
        # labels_df.columns = ["label"]
        labels_df.to_csv(sample_path / "labels.tsv", sep="\t", index_label="")

with tempfile.TemporaryDirectory() as temp_dir:
    for link in links:
        download_link(link["link"], link["name"], temp_dir)

    with tarfile.open(os.path.join(temp_dir, "All+SMI+Flat+data.tar.gz")) as tfile: # 11:30min
        tfile.extractall(temp_dir, filter="data")
    print("Raw data extracted")
    
    with tarfile.open(os.path.join(temp_dir, "All+SMI+Giotto+object.tar.gz")) as tfile: # 11:30min
        tfile.extractall(temp_dir, filter="data")
    print("RData extracted")

    rpy2.robjects.r['load'](os.path.join(temp_dir, "Processed Data Giotto Object/SMI_Giotto_Object.RData"))
    print("RData loaded")

    rpy2.robjects.r('metadata=gem@cell_metadata$rna')
    metadata_r = rpy2.robjects.r("metadata")
    with (rpy2.robjects.default_converter + rpy2.robjects.pandas2ri.converter).context():
        metadata_r = rpy2.robjects.conversion.get_conversion().rpy2py(metadata_r)

    for sample_id, sample_name in enumerate(samples):
        print("Now processing", sample_name)

        counts = pd.read_csv(os.path.join(temp_dir, f"{sample_name}/{sample_name}-Flat_files_and_images/{sample_name}_exprMat_file.csv"))
        counts = counts.set_index((counts["fov"].astype(str) + "_" + counts["cell_ID"].astype(str)).values)

        metadata = pd.read_csv(os.path.join(temp_dir, f"{sample_name}/{sample_name}-Flat_files_and_images/{sample_name}_metadata_file.csv"))
        metadata = metadata.set_index((metadata["fov"].astype(str) + "_" + metadata["cell_ID"].astype(str)).values)

        metadata_r_sample = metadata_r[metadata_r["Run_Tissue_name"] == sample_name]
        metadata_r_sample = metadata_r_sample.set_index(metadata_r_sample["cell_ID"].str.rsplit("_",n=2).str[1:].str.join("_"))
    
        obs_names = metadata.index.intersection(counts.index).intersection(metadata_r_sample.index)
        metadata = metadata.loc[obs_names]
        counts = counts.loc[obs_names]
        var_names = counts.columns[2:]
        var_names = var_names[~var_names.str.startswith("NegPrb")]

        metadata["cell_type"] = metadata_r_sample["cell_type"]
        metadata["niche"] = metadata_r_sample["niche"]

        X = counts[var_names].values

        spatial_df = metadata[["CenterX_global_px", "CenterY_global_px"]]
        spatial_df.columns = ["x", "y"]

        observations_df = pd.DataFrame(index=obs_names)
        observations_df["selected"] = True
        features_df = pd.DataFrame(index=var_names)
        features_df["selected"] = True

        labels_df = metadata[["niche", "cell_type"]].astype("category")
        labels_df.columns = ["label", "cell_type"]

        samples_df.loc[len(samples_df)] = [patients[sample_id], sample_name, np.nan, replicates[sample_id], len(labels_df["label"].cat.categories), sample_name]

        write_sample(out_dir,
                    sample_name,
                    spatial_df,
                    observations_df,
                    features_df,
                    X,
                    labels_df)

technology = "cosmx"

## Metadata files
samples_df.loc[
    :, ["patient", "sample", "position", "replicate", "directory", "n_clusters"]
].to_csv(out_dir / "samples.tsv", sep="\t", index_label="")

with open(out_dir / "experiment.json", "w") as f:
    exp_info = {"technology": technology}
    json.dump(exp_info, f)
