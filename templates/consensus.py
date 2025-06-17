#!/usr/bin/env python

# Author_and_contribution: Jieran Sun & Mark Robinson; implmented method
# Author_and_contribution: Peiying Cai; created template
# Author_and_contribution: ENTER YOUR NAME AND CONTRIBUTION HERE

import argparse

# TODO adjust description
parser = argparse.ArgumentParser(description="Calculate consensus ... for selected BCs")

parser.add_argument(
    "-i", "--input_file", help="Input containing the aggregated labels.", required=True)
parser.add_argument(
    "--seed", type=int, default=None, help="Seed for random number generator")
parser.add_argument(
    "-b", "--base_clusterings", help="Path to base-clustering ranking file", required=True)
parser.add_argument(
    "-o", "--output_file", help="Desired output file", required=True)
# TODO adjust default numbers in `n_bcs` and `n_clusters`
# make sure that `n_clusters` exists among the column names of the base-clustering ranking file
parser.add_argument(
    "--n_clusters", type=int, default=7, help="Desired number of clusters in the consensus output")
parser.add_argument(
    "--n_bcs", type=int, default=8, help="Desired number of base clustering results fed into the algorithm")

args = parser.parse_args()
from pathlib import Path
import pandas as pd
import sys
import warnings

seed = args.seed
output_file = args.output_file
output_path = Path(args.output_file)

# Read input label data
label_df = pd.read_csv(args.input_file, sep="\t", index_col=0)
    
# Read base clustering rankings
bc_df = pd.read_csv(args.base_clusterings, sep="\t", index_col=0)

n_clust_str = str(args.n_clusters)
if n_clust_str not in bc_df.columns:
    sys.exit(f"Error: n_clusters={args.n_clusters} not found in base clustering file columns.")

# bc_list stores all 'method_config_n_clust_label' entries matching n_clust
bc_list = bc_df[n_clust_str].dropna().tolist()

if len(bc_list) < args.n_bcs:
    warnings.warn(f"Not enough ({args.n_bcs}) base clusterings (BCs) are available, use {len(bc_list)} BCs instead.")
bc_list = bc_list[:min(args.n_bcs, len(bc_list))]

# Subset the label data to keep only the selected base clusterings
label_selected = label_df[bc_list]

# Make sure clusters are ranked 1 to n without jumps (SOTIP)
def rank_labels(u):
    unique_labels = sorted(u.dropna().unique())
    if unique_labels == list(range(1, len(unique_labels) + 1)):
        # Already consecutive starting from 1
        return u.astype('Int64')
        
    freq = u.dropna().value_counts()
    # Rank by descending frequency, ties.method='first' equivalent
    rank_map = {label: rank+1 for rank, label in enumerate(freq.index)}
    return u.map(rank_map).astype('Int64')

label_selected = label_selected.apply(rank_labels, axis=0)


# TODO set the seed, if the algorithm requires the seed elsewhere please pass it on
import random

random.seed(seed)
# np.random.seed(seed)
# torch.manual_seed(seed)

## Your code goes here
# TODO
# Input: label_selected (DataFrame) with samples as rows and base clusterings as columns

# output_df = ...

## Write output
output_path.parent.mkdir(parents=True, exist_ok=True)
output_df.to_csv(output_file, sep="\t", index=True)

