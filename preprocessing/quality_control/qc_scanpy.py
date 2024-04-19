#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Jieran Sun; Implemented quality control code

import argparse

# TODO adjust description
parser = argparse.ArgumentParser(description="quality control (gene/cell filtering) using scanpy")

parser.add_argument(
    "-c", "--coordinates", help="Path to coordinates (as tsv).", required=True
)
parser.add_argument("-m", "--matrix", help="Path to counts (as mtx).", required=True)
parser.add_argument(
    "-f", "--features", help="Path to features (as tsv).", required=True
)
parser.add_argument(
    "-o", "--observations", help="Path to observations (as tsv).", required=True
)
parser.add_argument("-d", "--out_dir", help="Output directory.", required=True)
parser.add_argument(
    "--config",
    help="Optional config file (json) used to pass additional parameters.",
    required=False,
)
parser.add_argument(
    "--min_genes",
    help="Minimum number of genes expressed required for a cell to pass filtering",
    required=False,
)
parser.add_argument(
    "--min_cells",
    help="Minimum number of cell expressed required for a gene to pass filtering",
    required=False,
)
parser.add_argument(
    "--min_counts",
    help="Minimum number of counts required for a cell to pass filtering",
    required=False,
)

args = parser.parse_args()

from pathlib import Path

out_dir = Path(args.out_dir)

# Output files
filter_counts_file = out_dir / "counts.mtx"
filter_obserations_file = out_dir / "observations.tsv"
filter_features_file = out_dir / "features.tsv"
filter_coordinates_file = out_dir / "coordinates.tsv"
# if additional output files are required write it also to out_dir

# Use these filepaths as input ...
coord_file = args.coordinates
matrix_file = args.matrix
feature_file = args.features
observation_file = args.observations

if args.config is not None:
    config_file = args.config

# Set default to 0 (no qc) if no info is provided
min_cells = int(args.min_cells) if args.min_cells is not None else 0
min_genes = int(args.min_genes) if args.min_genes is not None else 0
min_counts = int(args.min_counts) if args.min_counts is not None else 0

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
adata.var_names_make_unique()

## Your code goes here
import scanpy as sc
import pandas as pd
import numpy as np

sc.pp.filter_cells(adata, min_genes=min_genes)
sc.pp.filter_cells(adata, min_counts=min_counts)
sc.pp.filter_genes(adata, min_cells=min_cells)

# Find out mitochondria percentage 
mito_detect_list = [adata.var[col].astype(str).str.startswith("MT-") 
                    for col in adata.var.columns
                   ] + [adata.var_names.astype(str).str.startswith("MT-")]

mito_sum = [np.sum(mito_detect) for mito_detect in mito_detect_list]

if not all(x == 0 for x in mito_sum):
    adata.var["mt"] = mito_detect_list[np.argmax(mito_sum)]

    # Calculate QC metric
    sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False
    )

    # mitochondira needs to be less than 30%
    adata = adata[adata.obs["pct_counts_mt"] <= 30, ]

# Filter if tissue is on spot
if "in_tissue" in adata.obs.columns:
    adata.obs["in_tissue"] = adata.obs["in_tissue"].astype("bool")
    adata = adata[adata.obs["in_tissue"], ]

# each spots needs to contain more than 10 cells
if "cell_count" in adata.obs.columns:
    adata.obs["cell_count"] = adata.obs["cell_count"].astype("int64")
    adata = adata[adata.obs["cell_count"] > 10, ]


filter_counts = adata.X
observation_df = adata.obs
feature_df = adata.var
coord_cols = pd.read_table(args.coordinates, index_col=0).columns
coordinates_df = pd.DataFrame(adata.obsm["spatial"], index = adata.obs_names, columns = coord_cols)


## Write output
import scipy as sp

out_dir.mkdir(parents=True, exist_ok=True)
sp.io.mmwrite(filter_counts_file, filter_counts, precision=5)
observation_df.to_csv(filter_obserations_file, sep="\t", index_label="")
feature_df.to_csv(filter_features_file, sep="\t", index_label="")
coordinates_df.to_csv(filter_coordinates_file, sep="\t", index_label="")