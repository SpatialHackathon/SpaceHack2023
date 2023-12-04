#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: ENTER YOUR NAME AND CONTRIBUTION HERE

import argparse

# TODO adjust description
parser = argparse.ArgumentParser(description="Load data for ...")

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
# TODO
# features_df = ...  # DataFrame with index (gene-id/name) and n columns (?)
# observations_df = ...  # DataFrame with index (cell-id/barcode) and n columns (?)
# coordinates_df = ...  # DataFrame with index (cell-id/barcode) and 2/3 columns (x, y, z?)
# counts = ...  # array with #observations rows x #features columns
# labels_df = None  # optional, DataFrame with index (cell-id/barcode) and 1 column (label)
# img = None  # optional
# technology = ...  # "Visium", "ST", "imaging-based"
# samples_df = ...  # DataFrame with information on samples. columns: (patient, sample, position, replicate, directory, n_clusters), columns can be NA

# Make sure to use consistent indexes for the DataFrames
# i.e. the index (not necessarily the order) of observations and coordinates should match
# But the order of observations and features must match counts (observations x features)


# Example how a sample could be written
def write_sample(
    path,
    sample,
    coordinates_df,
    observations_df,
    features_df,
    counts,
    labels_df=None,
    img=None,
):
    if img is not None:
        # TODO write to image_file
        # H_E.json must contain the scale
        pass

    import scipy as sp

    sample_path = Path(path) / sample

    coordinates_df.to_csv(sample_path / "coordinates.tsv", sep="\t", index_label="")
    features_df.to_csv(sample_path / "features.tsv", sep="\t", index_label="")
    observations_df.to_csv(sample_path / "observations.tsv", sep="\t", index_label="")
    sp.io.mmwrite(sample_path / "counts.mtx", counts)

    if labels_df is not None:
        labels_df.columns = ["label"]
        labels_df.to_csv(sample_path / "labels.tsv", sep="\t", index_label="")


## Metadata files
samples_df.loc[
    :, ["patient", "sample", "position", "replicate", "directory", "n_clusters"]
].to_csv(out_dir / "samples.tsv", sep="\t", index_label="")

import json

with open(out_dir / "experiment.json", "w") as f:
    exp_info = {"technology": technology}
    json.dump(exp_info, f)
