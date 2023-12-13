#!/usr/bin/env python

# Author_and_contribution: Qirong Mao; created script

import argparse

# TODO adjust description
parser = argparse.ArgumentParser(description="SVG selection using Squidpy (Moran's I)")

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
     "-n", "--n_top_genes", help="Number of genes to keep (Default:3000).", required=False, type=int, default=3000)

parser.add_argument("-d", "--out_dir", help="Output directory.", required=True)

parser.add_argument(
    "--config",
    help="Optional config file (json) used to pass additional parameters.",
    required=False,
)

args = parser.parse_args()

from pathlib import Path
import pandas as pd

out_dir = Path(args.out_dir)

# Output files
feature_selection_file = out_dir / "features_MoranI.tsv"
# if additional output files are required write it also to out_dir

# Use these filepaths and inputs ...
coord_file = args.coordinates
matrix_file = args.matrix
feature_file = args.features
observation_file = args.observations


if args.n_top_genes is not None:
    n_top_genes = int(args.n_top_genes)
    
    
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
import squidpy as sq

## Warning if the input dataset has less features than the number of genes want to keep

n_input_features = len(adata.var)

if n_input_features < n_top_genes:
    raise ValueError('Input data has less features than the number of genes you want to keep, please clarify the numbers of genes you need to keep (--n_top_genes)')

features_df = adata.var.copy()

## Compute the graph based from Delaunay triangulation
sq.gr.spatial_neighbors(adata,coord_type="generic",delaunay= True)
## Calculate Moran's I score for each gene
sq.gr.spatial_autocorr(adata, mode="moran", genes=adata.var_names,show_progress_bar=True)


## Selecting top spatially variable genes based on Moran's I score
SVG = (
    adata.uns["moranI"]["I"].sort_values(ascending=False).head(n_top_genes).index.tolist()
)

### Create boolean indication to specify spatially variable genes      
adata.var['spatially_variable'] = adata.var.index.isin(SVG)

features_df["spatially_variable"] = adata.var["spatially_variable"]


## Write output
out_dir.mkdir(parents=True, exist_ok=True)
features_df.to_csv(feature_selection_file, sep="\t", index_label="")
