#!/usr/bin/env python

# Made by Paul Kiessling pakiessling@ukaachen.de

import argparse
import json
import os
import tempfile
import urllib.request

import anndata
import pandas as pd
import scipy

LINKS = {
    "slice_1": "https://gene.ai.tencent.com/SpatialOmics/api/download_data/2492",
    "slice_2": "https://gene.ai.tencent.com/SpatialOmics/api/download_data/2494",
    "slice_3": "https://gene.ai.tencent.com/SpatialOmics/api/download_data/2493",
}

META_DICT = {"technology": "BARseq2"}
LICENSE = """CC BY 4.0"""
SAMPLE_COLUMNS = [
    "patient",
    "sample",
    "position",
    "replicate",
    "n_clusters",
    "directory",
]


def download_links(links, temp_dir):
    for name, link in links.items():
        try:
            response = urllib.request.urlopen(link)

            filename = os.path.join(temp_dir, name + ".h5ad")
            with open(filename, "wb") as file:
                file.write(response.read())
            print(f"Downloaded: {filename}")
        except Exception as e:
            print(f"Error downloading {link}: {e}")


def process_adata(adata_path, output_folder, iteration, sample_df):
    folder_name = os.path.splitext(os.path.basename(adata_path))[0]
    complete_path = os.path.join(output_folder, folder_name)
    os.makedirs(complete_path, exist_ok=True)
    adata = anndata.read_h5ad(adata_path)

    # Observations
    obs = adata.obs.copy()
    obs["selected"] = "true"
    obs.loc[obs.layer == "outside_VISp", "selected"] = "false"
    obs.to_csv(f"{complete_path}/observations.tsv", sep="\t", index_label="")

    # Features
    vars = adata.var.copy()
    vars["selected"] = "true"
    vars.to_csv(f"{complete_path}/features.tsv", sep="\t", index_label="")

    # Coordinates
    coords = pd.DataFrame(adata.obsm["spatial"], columns=["x", "y"])
    coords.index = adata.obs.index
    coords.to_csv(f"{complete_path}/coordinates.tsv", sep="\t", index_label="")

    # Matrix
    scipy.io.mmwrite(f"{complete_path}/counts.mtx", adata.X)

    # add info for sample.tsv
    # Your sample_data_basis dictionary
    sample_data_basis = {
        "patient": "0",
        "sample": iteration,
        "position": "0",
        "replicate": "0",
        "n_clusters": adata.obs.layer.nunique()
        - 1,  # -1 for outside_VISp which we exclude
        "directory": folder_name,
    }

    # Concatenating the new DataFrame to sample_df
    sample_df.iloc[iteration] = sample_data_basis
    # Write labels.tsv
    labels = adata.obs.loc[adata.obs.layer != "outside_VISp", "layer"]
    labels.to_csv(f"{complete_path}/labels.tsv", sep="\t", index_label="")


def write_json(dict, output_path):
    with open(output_path, "w") as json_file:
        json.dump(dict, json_file)


def main():
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(
        description="Convert Barseq2 Viscortex data to Spacehack format."
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

        sample_df = pd.DataFrame(columns=SAMPLE_COLUMNS, index=range(len(LINKS)))
        anndatas = [
            os.path.join(temp_dir, file)
            for file in os.listdir(temp_dir)
            if file.endswith(".h5ad")
        ]
        for iteration, adata in enumerate(anndatas):
            process_adata(adata, args.out_dir, iteration, sample_df)

        # write json
        write_json(META_DICT, f"{args.out_dir}/experiment.json")

        # write samples.tsv
        sample_df.to_csv(f"{args.out_dir}/samples.tsv", sep="\t", index_label=False)

        # write LICENSE
        with open(f"{args.out_dir}/LICENSE.md", "w") as file:
            file.write(LICENSE)


if __name__ == "__main__":
    main()
