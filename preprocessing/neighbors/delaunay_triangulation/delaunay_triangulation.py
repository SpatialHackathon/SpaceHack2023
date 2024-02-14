#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created script
# Author_and_contribution: Qirong Mao; modifying output file format

import argparse

# TODO adjust description
parser = argparse.ArgumentParser(
    description="Neighbor definition using Delaunay triangulation"
)

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

spatial_connectivities_file = out_dir / "spatial_connectivities.mtx"
##spatial_distances_file = out_dir / "spatial_distances.mtx"


# Use these filepaths and inputs ...
coord_file = args.coordinates
matrix_file = args.matrix
feature_file = args.features
observation_file = args.observations


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
import squidpy as sq

sq.gr.spatial_neighbors(adata, delaunay=True, coord_type="generic")

neighbors = adata.obsp["spatial_connectivities"].astype(int)
##distance = adata.obsp["spatial_distances"].astype(float)

## Write output
import scipy as sp

out_dir.mkdir(parents=True, exist_ok=True)

sp.io.mmwrite(spatial_connectivities_file, neighbors)
##sp.io.mmwrite(spatial_distances_file, distance)            
