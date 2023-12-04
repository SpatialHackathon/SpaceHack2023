#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created script

import argparse

# TODO adjust description
parser = argparse.ArgumentParser(description="HVG selection using scanpy")

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
# parser.add_argument(
#     "-n", "--n_top_genes", help="Number of genes to keep.", required=False, type=int
# )
parser.add_argument("-d", "--out_dir", help="Output directory.", required=True)
parser.add_argument(
    "--config",
    help="Optional config file (json) used to pass additional parameters.",
    required=False,
)

args = parser.parse_args()

from pathlib import Path

out_dir = Path(args.out_dir)

# Output files
feature_selection_file = out_dir / "features.tsv"
# if additional output files are required write it also to out_dir

# Use these filepaths and inputs ...
coord_file = args.coordinates
matrix_file = args.matrix
feature_file = args.features
observation_file = args.observations

# if args.n_top_genes is not None:
#     n_top_genes = args.n_top_genes
if args.config is not None:
    config_file = args.config


# ... or AnnData if you want
def get_anndata(args):
    # Untested template
    import anndata as ad
    import pandas as pd
    import scipy as sp

    X = sp.io.mmread(args.matrix)
    if sp.sparse.issparse(X):
        X = X.tocsr()
    observations = pd.read_table(args.observations, index_col=0)
    features = pd.read_table(args.features, index_col=0)
    coordinates = (
        pd.read_table(args.coordinates, index_col=0)
        .loc[observations.index, :]
        .to_numpy()
    )

    adata = ad.AnnData(
        X=X, obs=observations, var=features, obsm={"spatial": coordinates}
    )

    return adata


adata = get_anndata(args)

## Your code goes here
import scanpy as sc

features_df = adata.var.copy()

sc.pp.highly_variable_genes(adata, flavor="seurat")

features_df["selected"] = adata.var["highly_variable"]


## Write output
out_dir.mkdir(parents=True, exist_ok=True)
features_df.to_csv(feature_selection_file, sep="\t", index_label="")
