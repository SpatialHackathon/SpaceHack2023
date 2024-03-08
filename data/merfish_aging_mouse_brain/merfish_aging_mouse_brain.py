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
    "aging": "https://datasets.cellxgene.cziscience.com/1892a5fe-c5f4-4fa7-9dac-54216566bdbf.h5ad",
    "inflammation": "https://datasets.cellxgene.cziscience.com/2bac088e-cc19-461f-a3dd-83466c44a122.h5ad",
}


META_DICT = {"technology": "Merfish"}
LICENSE = """CC-BY 4.0"""
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


def process_adata(adata_path, output_folder, iteration, samples_list):
    adata = anndata.read_h5ad(adata_path)
    for sample in adata.obs.donor_id.unique():
        print(f"Writing {sample}")
        # both files have donor 12 so we need to differentiate them
        complete_path = os.path.join(output_folder, f"Cohort_{iteration}_" + sample)
        os.makedirs(complete_path, exist_ok=True)
        subset = adata[adata.obs.donor_id == sample]

        # Observations
        obs = subset.obs.copy()
        obs["selected"] = "true"
        obs.to_csv(f"{complete_path}/observations.tsv", sep="\t", index_label="")

        # Features
        vars = subset.var.copy()
        vars["selected"] = "true"
        vars.to_csv(f"{complete_path}/features.tsv", sep="\t", index_label="")

        # Coordinates
        if "spatial" in subset.obsm:
            coords = pd.DataFrame(subset.obsm["spatial"], columns=["x", "y"])
        else:
            coords = pd.DataFrame(subset.obs[["center_x", "center_y"]])
        coords.index = subset.obs.index
        coords.to_csv(f"{complete_path}/coordinates.tsv", sep="\t", index_label="")

        # Matrix
        scipy.io.mmwrite(f"{complete_path}/counts.mtx", subset.X)

        # Write labels.tsv
        labels = subset.obs["tissue"]
        labels.to_csv(f"{complete_path}/labels.tsv", sep="\t", index_label="")

        # add info for sample.tsv
        # Your sample_data_basis dictionary
        sample_data_basis = {
            "sample": "0",
            "position": "0",
            "replicate": "0",
            "n_clusters": subset.obs.tissue.nunique(),
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
        description="Convert Merfish aging mouse brain to Spacehack format."
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
            os.path.join(temp_dir, file)
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
        sample_df.insert(0, "patient", range(len(sample_df)))
        sample_df.to_csv(f"{args.out_dir}/samples.tsv", sep="\t", index_label=False)

        # write LICENSE
        with open(f"{args.out_dir}/LICENSE.md", "w") as file:
            file.write(LICENSE)


if __name__ == "__main__":
    main()
