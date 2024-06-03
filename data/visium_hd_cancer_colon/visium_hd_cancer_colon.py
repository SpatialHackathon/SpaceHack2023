#!/usr/bin/env python

# Made by Paul Kiessling pakiessling@ukaachen.de


import os
import argparse
import tempfile
import shutil
from spatialdata_io import visium_hd
import pandas as pd
import scipy
import json
from pypdl import Downloader

LINKS = {
    "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_tissue_image.btf": "83bea5bc5761ccac5e54f7dc9b36d1d2",
    "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_feature_slice.h5": "ba5db688c27fa4b7203dc4f40c240453",
    "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_spatial.tar.gz": "18c123015ecad7dbb17e5862b427c21a",
    "https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Human_Colon_Cancer/Visium_HD_Human_Colon_Cancer_binned_outputs.tar.gz": "2a8a0df135d3d77ed77a465882f0bb2f",
    "https://zenodo.org/records/11402686/files/segmented_nuclei.zip": "952134e629e3861af44f70ccfc555a9d",
    "https://zenodo.org/records/11402686/files/2um_squares_annotation.csv": "f074d2a88b87bbd2d92d210bf9635849",
    "https://zenodo.org/records/11402686/files/8um_squares_annotation.csv": "e7eb740df3072cf7df505d61e8b62a9d",
    "https://zenodo.org/records/11402686/files/16um_squares_annotation.csv": "0c3644e6b6dcc9fa00d135bde6b8413b",
}

META_DICT = {"technology": "Visium HD"}
LICENSE = """

"""


def download_links(links, temp_dir):
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 6.1; Win64; x64; rv:47.0) Gecko/20100101 Firefox/47.0"
    }
    dl = Downloader(headers=headers)
    for link, checksum in links.items():
        print(f"Downloading {link}")
        file = dl.start(
            url=link,
            file_path=temp_dir,
            segments=10,
            display=True,
            multithread=True,
            block=True,
            retries=3,
        )
        if not file.validate_hash(checksum, "md5"):
            raise ValueError(f"File {file} is corrupted")

    # Extract the tar.gz files
    for file in os.listdir(temp_dir):
        if file.endswith(".tar.gz") or file.endswith(".zip"):
            shutil.unpack_archive(os.path.join(temp_dir, file), temp_dir)
    # Move binned folders to parent directory for visium_hd reader
    source_folder = os.path.join(temp_dir, "binned_outputs")
    parent_folder = os.path.dirname(source_folder)

    for folder in os.listdir(source_folder):
        item_path = os.path.join(source_folder, folder)
        if os.path.isdir(item_path):
            destination_path = os.path.join(parent_folder, folder)
            shutil.move(item_path, destination_path)


def process_adata(input_path, output_folder):
    results = {
        "square_002um": "2um_squares_annotation.csv",
        "square_008um": "8um_squares_annotation.csv",
        "square_016um": "16um_squares_annotation.csv",
    }
    for result in results.keys():
        os.makedirs(
            os.path.join(output_folder, f"visium_hd_cancer_colon_{result}"),
            exist_ok=True,
        )
    sdata = visium_hd(input_path)

    for result, annotation_file in results.items():
        table = sdata.tables[result]
        table.var_names_make_unique()
        domain_annotation = pd.read_table(
            os.path.join(input_path, annotation_file),
            index_col=0,
            header=None,
            names=["annot_type"],
        )
        complete_path = os.path.join(output_folder, f"visium_hd_cancer_colon_{result}")

        # Obs
        obs = table.obs.copy()
        obs["selected"] = "true"
        # A few bins are outside the image + disconnected, so we remove them
        obs.loc[
            domain_annotation.loc[:, "annot_type"] == "Outside", "selected"
        ] = "false"

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


def process_nuclei(input_path, output_folder):
    input = os.path.join(input_path, "output")
    output = os.path.join(output_folder, "visium_hd_cancer_colon_segmented_nuclei")
    os.makedirs(output, exist_ok=True)

    for filename in os.listdir(input):
        file_path = os.path.join(input, filename)
        if os.path.isfile(file_path):
            shutil.move(file_path, os.path.join(output, filename))
    shutil.copy(
        os.path.join(input_path, "Visium_HD_Human_Colon_Cancer_tissue_image.btf"),
        os.path.join(output, "H_E.tiff"),
    )


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
    with tempfile.TemporaryDirectory() as temp_dir:  #
        download_links(LINKS, temp_dir)
        os.makedirs(args.out_dir, exist_ok=True)

        sample_df = {
            "patient": [1, 1, 1, 1],
            "sample": [0, 0, 0, 0],
            "position": [1, 1, 1, 1],
            "replicate": [1, 1, 1, 1],
            "n_clusters": [6, 6, 6, 6],
            "directory": [
                "visium_hd_cancer_colon_segmented_nuclei",
                "visium_hd_cancer_colon_square_002um",
                "visium_hd_cancer_colon_square_008um",
                "visium_hd_cancer_colon_square_016um",
            ],
        }
        sample_df = pd.DataFrame(sample_df)

        # Segmented Nuclei are already in the correct format
        process_nuclei(temp_dir, args.out_dir)

        process_adata(temp_dir, args.out_dir)

        # write json
        write_json(META_DICT, f"{args.out_dir}/experiment.json")

        # write samples.tsv
        sample_df.to_csv(f"{args.out_dir}/samples.tsv", sep="\t", index_label=False)

        # write LICENSE
        with open(f"{args.out_dir}/LICENSE.md", "w") as file:
            file.write(LICENSE)


if __name__ == "__main__":
    main()
