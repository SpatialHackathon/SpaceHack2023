#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created script

import argparse

# TODO adjust description
parser = argparse.ArgumentParser(description="PCA (with standard-scaling)")

parser.add_argument(
    "-c", "--coordinates", help="Path to coordinates (as tsv).", required=True
)
parser.add_argument(
    "-m", "--matrix", help="Path to (transformed) counts (as mtx).", required=True
)
parser.add_argument(
    "-f", "--features", help="Path to features (as tsv).", required=True
)
parser.add_argument(
    "-o", "--observations", help="Path to observations (as tsv).", required=True
)
parser.add_argument(
    "-n",
    "--n_components",
    help="Number of components/factors to generate.",
    required=False,
    type=int,
)
parser.add_argument("-d", "--out_dir", help="Output directory.", required=True)
parser.add_argument(
    "--config",
    help="Optional config file (json) used to pass additional parameters.",
    required=False,
)

args = parser.parse_args()

# Output files
from pathlib import Path

out_dir = Path(args.out_dir)

# Output files
dim_red_file = out_dir / "dimensionality_reduction.tsv"
# if additional output files are required write it also to out_dir

# Use these filepaths as input ...
coord_file = args.coordinates
matrix_file = args.matrix
feature_file = args.features
observation_file = args.observations

if args.n_components is not None:
    n_components = args.n_components
if args.config is not None:
    config_file = args.config

## Your code goes here
import pandas as pd
import scipy as sp
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

matrix = sp.io.mmread(matrix_file)
if sp.sparse.issparse(matrix):
    matrix = matrix.tocsr()

observations = pd.read_table(observation_file, index_col=0)
features = pd.read_table(feature_file, index_col=0)

# Filter features and observations
if "selected" in observations.columns:
    matrix = matrix[observations["selected"].to_numpy().nonzero()[0], :]
    observations = observations.loc[lambda df: df["selected"]].index
else:
    observations = observations.index
if "selected" in features.columns:
    matrix = matrix[:, features["selected"].to_numpy().nonzero()[0]]
    features = features.loc[lambda df: df["selected"]].index
else:
    features = features.index

matrix = matrix.toarray() if sp.sparse.issparse(matrix) else matrix
matrix = pd.DataFrame(matrix, columns=features, index=observations)

scaler = StandardScaler().set_output(transform="pandas")
matrix = scaler.fit_transform(matrix)

pca = PCA(n_components=n_components).set_output(transform="pandas")
dim_red_df = pca.fit_transform(matrix)


## Write output
out_dir.mkdir(parents=True, exist_ok=True)
dim_red_df.to_csv(dim_red_file, sep="\t", index_label="", float_format="%g")
