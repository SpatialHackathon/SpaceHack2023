#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Thomas Chartrand; created script

import argparse

parser = argparse.ArgumentParser(description="Calculate Matthew's correlation coefficient (scikit-learn)")

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
parser.add_argument(
    "--matched_labels",
    help="Flag indicating ground-truth and clustering labels have already been matched.",
    action='store_true',
)

args = parser.parse_args()

# Use these filepaths as input
label_file = args.labels

if args.ground_truth is not None:
    groundtruth_file = args.ground_truth
if args.config is not None:
    config_file = args.config


## Your code goes here
if args.ground_truth is None:
    raise Exception("Groundtruth labels needed.")

import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import matthews_corrcoef

domains = pd.read_table(label_file, index_col=0)["label"].astype("category").cat.codes
groundtruth = (
    pd.read_table(groundtruth_file, index_col=0)["label"].astype("category").cat.codes
)
common_index = domains.index.intersection(groundtruth.index)
groundtruth = groundtruth.loc[common_index]
domains = domains.loc[common_index]

if not args.matched_labels:
    contingency_table = pd.crosstab(domains, groundtruth)
    row_ind, col_ind = linear_sum_assignment(contingency_table, maximize=True)
    domains = domains.map(dict(zip(row_ind, col_ind)))

metric = matthews_corrcoef(groundtruth, domains)

## Write output
from pathlib import Path

Path(args.out_file).parent.mkdir(parents=True, exist_ok=True)

with open(args.out_file, "w") as file:
    file.write(f"{metric:.5e}\n")
