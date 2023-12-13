#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Kirti Biharie; implemented Calinski-Harabasz score

import argparse

parser = argparse.ArgumentParser(description="Calculate Calinski-Harabasz Score (scikit-learn)")

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
if args.embedding is None:
    raise Exception("Embeddings needed to calculate the Calinski-Harabasz Score")

import pandas as pd
import sklearn.metrics

embeddings = pd.read_table(embedding_file, index_col=0)
labels = pd.read_table(label_file, index_col=0)["label"].astype("category").cat.codes

common_index = labels.index.intersection(embeddings.index)
embeddings = embeddings.loc[common_index]
labels = labels.loc[common_index]

metric = sklearn.metrics.calinski_harabasz_score(embeddings, labels)

## Write output
from pathlib import Path

Path(args.out_file).parent.mkdir(parents=True, exist_ok=True)

with open(args.out_file, "w") as file:
    file.write(f"{metric:.5e}\n")