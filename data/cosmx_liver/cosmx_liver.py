# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Kirti Biharie; added dataset

import argparse

parser = argparse.ArgumentParser(description="Load data for cosmx liver dataset")

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
import tiledbsc
import anndata as ad
import zipfile
import subprocess

sample_columns = ["patient","sample","position","replicate","n_clusters","directory"]
samples_df = pd.DataFrame(columns=sample_columns)

links = [
         {"name": "LiverDataReleaseTileDB.zip", "link":"https://smi-public.objects.liquidweb.services/LiverDataReleaseTileDB.zip"}, 
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

    with zipfile.ZipFile(os.path.join(temp_dir, "LiverDataReleaseTileDB.zip"), 'r') as zObject: 
        zObject.extractall(path=temp_dir)

    collection = tiledbsc.soma_collection.SOMACollection(os.path.join(temp_dir, "LiverDataRelease/"))
    print(collection)

    counts = collection['RNA'].X['counts'].csr()

    obs = collection['RNA'].obs.df()
    obs = obs.set_index("obs_id")
    var = collection['RNA'].var.df()
    var = var.set_index("var_id")
    coordinates = obs[["x_slide_mm", "y_slide_mm"]]

    adata = ad.AnnData(counts, obs=obs, obsm={"spatial": coordinates.values}, var=var)

    for sample_id, sample_name in enumerate(["NormalLiver", "CancerousLiver"]):
        sample = adata[adata.obs["Run_Tissue_name"] == sample_name]
        
        X = sample.X
        coordinates_df = pd.DataFrame(sample.obsm["spatial"], index=sample.obs_names, columns=["x", "y"])
        labels_df = sample.obs[["niche", "cellType"]].astype("category")
        labels_df.columns = ["label", "cellType"]
        features_df = pd.DataFrame(index=sample.var_names)
        features_df["selected"] = True
        observations_df = pd.DataFrame(index=sample.obs_names)
        observations_df["selected"] = True
        samples_df.loc[len(samples_df)] = [sample_id, sample_name, np.nan, np.nan, len(labels_df["label"].cat.categories), sample_name]

        write_sample(out_dir,
                    sample_name,
                    coordinates_df,
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
