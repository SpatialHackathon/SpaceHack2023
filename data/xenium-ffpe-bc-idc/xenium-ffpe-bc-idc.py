#!/usr/bin/env python

# Made by Paul Kiessling pakiessling@ukaachen.de

import argparse
import json
import os
import shutil
import tempfile

import pandas as pd
import scipy
from pypdl import Downloader
from spatialdata_io import xenium

LINKS = {
    "https://s3-us-west-2.amazonaws.com/10x.files/samples/xenium/1.0.2/Xenium_V1_FFPE_Human_Breast_IDC_With_Addon/Xenium_V1_FFPE_Human_Breast_IDC_With_Addon_outs.zip": "7d3374472092b320ee9b876cb56c520b",
    "https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FFPE_Human_Breast_ILC_With_Addon/Xenium_V1_FFPE_Human_Breast_ILC_With_Addon_outs.zip": "cf779754817893dc98ff4311df2db61e",
    "https://zenodo.org/records/15411357/files/idc.csv": "219392e7c41587efaeb315d172fba7b0",
    "https://zenodo.org/records/15411357/files/ilc.csv": "d1b150c706539f0f20d9df5496a08984",
}


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


def process_files(temp_folder, out_path, annotation_path, sample_name):
    print(f"Transferring {temp_folder} to {out_path}")
    sdata = xenium(
        temp_folder,
        cells_boundaries=False,
        nucleus_boundaries=False,
        cells_as_circles=False,
        cells_labels=False,
        nucleus_labels=False,
        transcripts=False,
        morphology_mip=False,
        morphology_focus=False,
    )
    sdata = sdata["table"]
    annotation = pd.read_csv(annotation_path)
    sdata.obs = sdata.obs.merge(
        annotation, how="left", left_on="cell_id", right_on="cell"
    )
    sdata = sdata[sdata.obs["domain"].notna()].copy()
    process_adata(sdata, out_path, sample_name)


def process_adata(adata, out_path, sample_name):
    complete_path = os.path.join(out_path, sample_name)
    os.makedirs(complete_path, exist_ok=True)

    # Observations
    obs = adata.obs.copy()
    obs["selected"] = "true"
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

    # Write labels.tsv
    labels = adata.obs["domain"]
    labels.to_csv(f"{complete_path}/labels.tsv", sep="\t", index_label="")


def write_json(out_path):
    experiment = {
        "technology": "Xenium",
        "species": "human",
        "is_3D": False,
    }
    with open(os.path.join(out_path, "experiment.json"), "w") as f:
        json.dump(experiment, f)


def write_table(out_path):
    data = {
        "patient": ["ILC", "IDC"],
        "sample": ["ILC", "IDC"],
        "position": [0, 0],
        "replicate": [0, 0],
        "directory": ["ILC", "IDC"],
        "n_clusters": [6, 8],
    }
    df = pd.DataFrame(data)
    df.to_csv(f"{out_path}/samples.tsv", sep="\t", index_label=False)


def main():
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(
        description="Convert Xenium data to Spacehack format."
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
        for file in os.listdir(temp_dir):
            if file.endswith(".zip"):
                sample_name = file.split(".")[0]
                sample_path = os.path.join(temp_dir, sample_name)

                shutil.unpack_archive(os.path.join(temp_dir, file), sample_path)

                if sample_name == "Xenium_V1_FFPE_Human_Breast_IDC_With_Addon_outs":
                    annotation_path = os.path.join(temp_dir, "idc.csv")
                    sample_name_short = "IDC"

                elif sample_name == "Xenium_V1_FFPE_Human_Breast_ILC_With_Addon_outs":
                    annotation_path = os.path.join(temp_dir, "ilc.csv")
                    sample_name_short = "ILC"
                else:
                    raise ValueError("Unknown sample found in downloaded files")

                output_path = os.path.join(args.out_dir, sample_name_short)
                process_files(
                    sample_path, output_path, annotation_path, sample_name_short
                )
        write_json(args.out_dir)
        write_table(args.out_dir)


if __name__ == "__main__":
    main()
