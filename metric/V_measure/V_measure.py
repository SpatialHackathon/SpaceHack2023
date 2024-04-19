#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created script

import argparse

parser = argparse.ArgumentParser(description="Calculate V-measure (scikit-learn)")

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
if args.config is not None:
    config_file = args.config


## Your code goes here
if args.ground_truth is None:
    raise Exception("Groundtruth labels needed to calculate the Adjusted Rand Index.")

if args.config is None:
    raise Exception("Config file not provided.")

import json

import pandas as pd
import os
from sklearn.metrics import v_measure_score

with open(config_file, "r") as f:
    config = json.load(f)

domains = pd.read_table(label_file, index_col=0)["label"].astype("category").cat.codes
groundtruth = (
    pd.read_table(groundtruth_file, index_col=0)["label"].astype("category").cat.codes
)
common_index = domains.index.intersection(groundtruth.index)

metric = v_measure_score(
    groundtruth.loc[common_index], domains.loc[common_index], beta=config["beta"]
)

## Write output
out_file_path = args.out_file
out_dir = os.path.dirname(out_file_path)
os.makedirs(out_dir, exist_ok=True)

with open(out_file_path, "w") as file:
    file.write(f"{metric:.5e}\n")
