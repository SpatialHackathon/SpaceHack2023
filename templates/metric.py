#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: ENTER YOUR NAME AND CONTRIBUTION HERE

import argparse

# TODO adjust description
parser = argparse.ArgumentParser(description="Calculate metric ...")

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
    help="Optional config file (json) used to pass additional parameters.",
    required=False,
)  # format should be json
parser.add_argument("-o", "--out_file", help="Output file.", required=True)

args = parser.parse_args()

# Use these filepaths as input
label_file = args.labels

if args.ground_truth is not None:
    groundtruth_file = args.ground_truth
if args.embedding is not None:
    embedding_file = args.embedding
if args.config is not None:
    config_file = args.config


## Your code goes here
# TODO
# metric: float = ... output of the metric as float


## Write output
from pathlib import Path

Path(args.out_file).parent.mkdir(parents=True, exist_ok=True)

with open(args.out_file, "w") as file:
    file.write(f"{metric:.5e}\n")
