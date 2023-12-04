#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created script

import argparse

parser = argparse.ArgumentParser(
    description="Calculate Adjusted Rand Index (scikit-learn)"
)

parser.add_argument(
    "-l", "--labels", help="Labels from domain clustering.", required=True
)
parser.add_argument("-g", "--ground_truth", help="Groundtruth labels.", required=False)
parser.add_argument(
    "-e",
    "--embedding",
    help="Embedding of points in latent space. Potential usage for metrics without groundtruth.",
    required=False,
)
parser.add_argument(
    "-c",
    "--config",
    help="Optional config file used to pass additional parameters.",
    required=False,
)
parser.add_argument("-o", "--out_file", help="Output file.", required=True)

args = parser.parse_args()

# Use these filepaths as input
label_file = args.labels

if args.ground_truth is not None:
    groundtruth_file = args.ground_truth


## Your code goes here
if args.ground_truth is None:
    raise Exception("Groundtruth labels needed to calculate the Adjusted Rand Index")

import pandas as pd
from sklearn.metrics import adjusted_rand_score

domains = pd.read_table(label_file, index_col=0)["label"].astype("category").cat.codes
groundtruth = (
    pd.read_table(groundtruth_file, index_col=0)["label"].astype("category").cat.codes
)

common_index = domains.index.intersection(groundtruth.index)

metric = adjusted_rand_score(groundtruth.loc[common_index], domains.loc[common_index])


## Write output
with open(args.out_file, "w") as file:
    file.write(f"{metric:.5e}\n")
