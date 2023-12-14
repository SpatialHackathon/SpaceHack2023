#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created script, 
# Author_and_contribution: Qirong Mao; implemented method

import argparse

# TODO adjust description
parser = argparse.ArgumentParser(
    description="Constructing spatial neighborhood with custom parameters"
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

parser.add_argument(
    "--coord_type", default=None, choices=['grid', 'generic', None],help="Type of coordinate system.", required=False
)

parser.add_argument(
    "--n_neighs", type=int, default=6,help="Number of neighboring tiles (grid data) or neighborhoods (non-grid data)", required=False
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
spatial_distances_file = out_dir / "spatial_distances.mtx"

# Use these filepaths and inputs ...
coord_file = args.coordinates
matrix_file = args.matrix
feature_file = args.features
observation_file = args.observations

## Custom parameters
coord_type = args.coord_type
n_neighs= args.n_neighs

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
import squidpy as sq

sq.gr.spatial_neighbors(adata, coord_type=coord_type, n_neighs=n_neighs)

neighbors = adata.obsp["spatial_connectivities"].astype(int)
distance = adata.obsp["spatial_distances"].astype(int)

## Write output
import scipy as sp

out_dir.mkdir(parents=True, exist_ok=True)

sp.io.mmwrite(spatial_connectivities_file, neighbors)
sp.io.mmwrite(spatial_distances_file, distance)
             
