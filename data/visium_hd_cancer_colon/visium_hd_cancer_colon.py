#!/usr/bin/env python

# Made by Paul Kiessling pakiessling@ukaachen.de


import urllib.request
from urllib.parse import urlparse
import os
import argparse
import tempfile
import shutil
from spatialdata_io import visium_hd
import pandas as pd
import scipy
import json

LINKS = [
    "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_tissue_image.btf",
    "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_feature_slice.h5",
    "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_spatial.tar.gz",
    "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_square_002um_outputs.tar.gz",
    "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_square_008um_outputs.tar.gz",
    "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_square_016um_outputs.tar.gz",
    "https://zenodo.org/records/11402686/files/segmented_nuclei.zip",
    "https://zenodo.org/records/11402686/files/2um_squares_annotation.csv",
    "https://zenodo.org/records/11402686/files/8um_squares_annotation.csv",
    "https://zenodo.org/records/11402686/files/16um_squares_annotation.csv",
]

META_DICT = {"technology": "Visium HD"}
LICENSE = """

"""
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
            response = urllib.request.urlopen(
                urllib.request.Request(
                    link, headers={"User-Agent": "Mozilla"}
                )  # Add a User-Agent header to avoid HTTP 403 error
            )

            # Extract filename from the URL
            filename = os.path.join(temp_dir, urlparse(link).path.split("/")[-1])
            with open(filename, "wb") as file:
                file.write(response.read())
            print(f"Downloaded: {filename}")
        except Exception as e:
            print(f"Error downloading {link}: {e}")

    # Extract the tar.gz files
    for file in os.listdir(temp_dir):
        if file.endswith(".tar.gz") or file.endswith(".zip"):
            shutil.unpack_archive(os.path.join(temp_dir, file), temp_dir)

    print(os.listdir(temp_dir))


def process_adata(input_path, output_folder, sample_df):
    results = {"square_002um":"2um_squares_annotation.csv","square_008um" : "8um_squares_annotation.csv","square_016um":"16um_squares_annotation.csv"}
    for result in results.keys():
        os.makedirs(os.path.join(output_folder, f"visium_hd_cancer_colon_{result}"), exist_ok=True)
    sdata = visium_hd(input_path)

    for result,annotation_file in results.items():
        table = sdata.tables[result]
        table.var_names_make_unique()
        domain_annotation = pd.read_table(
            os.path.join(input_path, annotation_file),
            index_col=0,
            header=None,
            names=["annot_type"],
        )
        complete_path  = os.path.join(output_folder, f"visium_hd_cancer_colon_{result}")

        # Obs
        obs = table.obs.copy()
        obs["selected"] = "true"
        # A few bins are outside the image + disconnected, so we remove them
        obs.loc[domain_annotation.loc[:, "annot_type"] == "Outside", "selected"] = "false"

        obs.rename(columns={"array_row": "row", "array_col": "col"}, inplace=True)
        obs.to_csv(f"{complete_path}/observations.tsv", sep="\t", index_label="")

        # Features
        vars = table.var.copy()
        vars["selected"] = "true"
        vars.to_csv(f"{complete_path}/features.tsv", sep="\t", index_label="")

        # Coordinates
        coords = pd.DataFrame(table.obsm["spatial"], columns=["x", "y"])
        coords.index = table.obs.index
        coords.to_csv(f"{complete_path}/coordinates.tsv", sep="\t", index_label="")

        # Matrix
        scipy.io.mmwrite(f"{complete_path}/counts.mtx", table.X)

        # Write labels.tsv
        labels = domain_annotation.loc[:, "annot_type"]
        labels.index = table.obs.index
        labels = labels.rename("label")
        labels.to_csv(f"{complete_path}/labels.tsv", sep="\t", index_label="")

        # Move image files
        shutil.copy(
            os.path.join(input_path, "Visium_HD_Human_Colon_Cancer_tissue_image.btf"),
            os.path.join(complete_path, "H_E.tiff"),
    )


    # add info for sample.tsv
    # sample_data_basis = {
    #     "patient": "1",
    #     "sample": "1",
    #     "position": "0",
    #     "replicate": "1",
    #     "n_clusters": domain_annotation.loc[:, "annot_type"].nunique()
    #     - 1,  # -1 because we removed the "Outside" class
    #     "directory": "visium_hd_cancer_colon",
    # }

    # # Concatenating the new DataFrame to sample_df
    # sample_df.iloc[0] = sample_data_basis






def process_nuclei(input_path, output_folder):
    input = os.path.join(input_path, "output")
    output = os.path.join(output_folder, "visium_hd_cancer_colon_segmented_nuclei")
    shutil.move(input, output)


def write_json(dict, output_path):
    with open(output_path, "w") as json_file:
        json.dump(dict, json_file)


def main():
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(
        description="Convert Visium HD data to Spacehack format."
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
        sample_df = pd.DataFrame(columns=SAMPLE_COLUMNS, index=[0])

        # Segmented Nuclei are already in the correct format
        process_nuclei(temp_dir, args.out_dir)

        process_adata(temp_dir, args.out_dir, sample_df)

        # write json
        write_json(META_DICT, f"{args.out_dir}/experiment.json")

        # write samples.tsv
        sample_df.to_csv(f"{args.out_dir}/samples.tsv", sep="\t", index_label=False)

        # write LICENSE
        with open(f"{args.out_dir}/LICENSE.md", "w") as file:
            file.write(LICENSE)


if __name__ == "__main__":
    main()
