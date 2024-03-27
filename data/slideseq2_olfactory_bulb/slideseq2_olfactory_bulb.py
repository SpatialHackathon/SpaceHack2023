#!/usr/bin/env python

# Made by Paul Kiessling pakiessling@ukaachen.de


import argparse
import json
import os
import tempfile
import urllib.request
from urllib.parse import urlparse

import anndata
import pandas as pd
import scipy

LINKS = [
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000172/stomics/STSP0005373/GSM5173925_OB1_Slide1_processed.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000172/stomics/STSP0005374/GSM5173926_OB1_Slide2_processed.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000172/stomics/STSP0005375/GSM5173927_OB1_Slide3_processed.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000172/stomics/STSP0005376/GSM5173928_OB1_Slide4_processed.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000172/stomics/STSP0005377/GSM5173929_OB1_Slide5_processed.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000172/stomics/STSP0005379/GSM5173931_OB1_Slide7_processed.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000172/stomics/STSP0005382/GSM5173934_OB1_Slide10_processed.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000172/stomics/STSP0005378/GSM5173930_OB1_Slide6_processed.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000172/stomics/STSP0005381/GSM5173933_OB1_Slide9_processed.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000172/stomics/STSP0005380/GSM5173932_OB1_Slide8_processed.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000172/stomics/STSP0005383/GSM5173935_OB1_Slide11_processed.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000172/stomics/STSP0005384/GSM5173936_OB1_Slide12_processed.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000172/stomics/STSP0005385/GSM5173937_OB1_Slide13_processed.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000172/stomics/STSP0005386/GSM5173938_OB1_Slide14_processed.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000172/stomics/STSP0005387/GSM5173939_OB1_Slide15_processed.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000172/stomics/STSP0005388/GSM5173940_OB1_Slide16_processed.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000172/stomics/STSP0005389/GSM5173941_OB1_Slide17_processed.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000172/stomics/STSP0005390/GSM5173942_OB1_Slide18_processed.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000172/stomics/STSP0005391/GSM5173943_OB1_Slide19_processed.h5ad"
    "https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000172/stomics/STSP0005392/GSM5173944_OB1_Slide20_processed.h5ad"
]


META_DICT = {"technology": "Slideseq2"}
LICENSE = """ """
SAMPLE_COLUMNS = [
    "patient",
    "sample",
    "position",
    "replicate",
    "n_clusters",
    "directory",
]


def download_links(links, temp_dir):
    for link in links:
        try:
            response = urllib.request.urlopen(link)
            # Extract filename from the URL
            filename = os.path.join(temp_dir, urlparse(link).path.split("/")[-1])
            with open(filename, "wb") as file:
                file.write(response.read())
            print(f"Downloaded: {filename}")
        except Exception as e:
            print(f"Error downloading {link}: {e}")


def process_adata(adata_path, output_folder, iteration, samples_list):
    adata = anndata.read_h5ad(adata_path)

    # extract sample name from path
    sample = os.path.basename(adata_path).split("_")[0]
    complete_path = os.path.join(output_folder, sample)
    print(f"Writing {sample}")
    os.makedirs(complete_path, exist_ok=True)

    # Observations
    obs = adata.obs.copy()
    obs["selected"] = "true"
    obs.loc[obs.layer == "unknown", "selected"] = "false"
    obs.to_csv(f"{complete_path}/observations.tsv", sep="\t", index_label="")

    # Features
    vars = adata.var.copy()
    vars["selected"] = "true"
    vars.to_csv(f"{complete_path}/features.tsv", sep="\t", index_label="")

    # Coordinates
    if "spatial" in adata.obsm:
        coords = pd.DataFrame(adata.obsm["spatial"], columns=["x", "y"])
    coords.index = adata.obs.index
    coords.to_csv(f"{complete_path}/coordinates.tsv", sep="\t", index_label="")

    # Matrix
    scipy.io.mmwrite(f"{complete_path}/counts.mtx", adata.X)

    # Write labels.tsv
    labels = adata.obs["layer"]
    labels.to_csv(f"{complete_path}/labels.tsv", sep="\t", index_label="")

    # add info for sample.tsv
    sample_data_basis = {
        "patient": "0",
        "sample": "0",
        "position": iteration,
        "replicate": "0",
        "n_clusters": adata.obs.layer.nunique() - 1,  # -1 for unknown which we exclude
        "directory": complete_path,
    }

    # Append to sample list to create df later
    samples_list.append(sample_data_basis)


def write_json(dict, output_path):
    with open(output_path, "w") as json_file:
        json.dump(dict, json_file)


def main():
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(
        description="Convert Slideseqv2 olfactory bulb to Spacehack format."
    )

    # Add arguments for output folder
    parser.add_argument(
        "-o", "--out_dir", help="Output directory to write files to.", required=True
    )

    # Parse the command-line arguments
    args = parser.parse_args()

    # Download and process
    with tempfile.TemporaryDirectory() as temp_dir:
        download_links(LINKS, temp_dir)
        os.makedirs(args.out_dir, exist_ok=True)

        samples = []
        anndatas = [
            os.path.join(
                temp_dir,
                file,
            )
            for file in os.listdir(temp_dir)
            if file.endswith(".h5ad")
        ]
        for iteration, adata in enumerate(anndatas):
            process_adata(adata, args.out_dir, iteration, samples)

        sample_df = pd.DataFrame(
            samples,
            index=range(len(samples)),
        )

        # write json
        write_json(META_DICT, f"{args.out_dir}/experiment.json")

        # write samples.tsv
        sample_df.to_csv(f"{args.out_dir}/samples.tsv", sep="\t", index_label=False)

        # write LICENSE
        with open(f"{args.out_dir}/LICENSE.md", "w") as file:
            file.write(LICENSE)


if __name__ == "__main__":
    main()
