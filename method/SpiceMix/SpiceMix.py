#!/usr/bin/env python

# Author_and_contribution: Shahul Alam; created script

import argparse

parser = argparse.ArgumentParser(
    description="""SpiceMix (https://www.nature.com/articles/s41592-021-01255-8)
Requirements: Any spatially resolved transcriptomics."""
)

parser.add_argument(
    "-c", "--coordinates", help="Path to coordinates (as tsv).", required=True
)
parser.add_argument(
    "-m", "--matrix", help="Path to (transformed) counts (as mtx).", required=True
)
# parser.add_argument(
#     "-f", "--features", help="Path to features (as tsv).", required=True
# )
parser.add_argument(
    "-o", "--observations", help="Path to observations (as tsv).", required=True
)

parser.add_argument(
    "-n",
    "--neighbors",
    help="Path to neighbor definitions. Square matrix (not necessarily symmetric) where each row contains the neighbors of this observation (as mtx).",
    required=False,
)

parser.add_argument("-d", "--out_dir", help="Output directory.", required=True)

parser.add_argument(
    "--n_clusters", help="Number of clusters to return.", required=True, type=int
)
parser.add_argument(
    "--seed", help="Seed to use for random operations.", required=True, type=int
)
parser.add_argument(
    "--config",
    help="Config file used to pass additional parameters.",
    required=True,
)

args = parser.parse_args()

from pathlib import Path

out_dir = Path(args.out_dir)
out_dir.mkdir(parents=True, exist_ok=True)

label_file = out_dir / "domains.tsv"

## Your code goes here
import json
import tempfile

with open(args.config, "r") as f:
    config = json.load(f)

import numpy as np
import pandas as pd

import anndata as ad
import scanpy as sc

from scipy.io import mmread
from scipy.sparse import issparse

import torch

from popari.components import PopariDataset
from popari.io import save_anndata
from popari.model import SpiceMix
from popari import tl

preprocess_parameters = config.pop("preprocess", None)
def get_anndata(args):
    """Convert data input into SpiceMix anndata format."""
    
    counts_matrix = mmread(args.matrix)

    if issparse(counts_matrix):
        counts_matrix = counts_matrix.tocsr()

    observations = pd.read_csv(args.observations, delimiter="\t", index_col=0)
    coordinates = pd.read_csv(args.coordinates, delimiter="\t", index_col=0)

    adata = ad.AnnData(
        X=counts_matrix, obs=observations, obsm={"spatial": coordinates[["x", "y"]].values}
    )

    if preprocess_parameters:
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=preprocess_parameters["hvgs"])
        
        adata = PopariDataset(adata[:, adata.var["highly_variable"]], "processed") 
        adata.compute_spatial_neighbors()
        
    else:
        adjacency_matrix = mmread(args.neighbors)
        if issparse(adjacency_matrix):
            adjacency_matrix = adjacency_matrix.tocsr()
    
        # Symmetrize matrix
        transpose_condition = adjacency_matrix.T > adjacency_matrix
        adjacency_matrix = adjacency_matrix + adjacency_matrix.T.multiply(transpose_condition) - adjacency_matrix.multiply(transpose_condition)
        
        num_cells = adjacency_matrix.shape[0]
        adjacency_list = [[] for _ in range(num_cells)]                                         
        for x, y in zip(*adjacency_matrix.nonzero()):                              
            adjacency_list[x].append(y)                                                         

        adata.obsp["adjacency_matrix"] = adjacency_matrix                                                                                         
        adata.obs["adjacency_list"] = adjacency_list 
    
        adata = PopariDataset(adata, "raw")

    return adata

adata = get_anndata(args)

num_preiterations = config.pop("num_preiterations", 5)
num_iterations = config.pop("num_iterations", 200)

device = config.pop("device", "cpu")
dtype = config.pop("dtype", "float32")

if dtype == "float32":
    dtype = torch.float32
elif dtype == "float64":
    dtype = torch.float64

torch_context = {
    "device": device,
    "dtype": dtype
}




temp_dir = tempfile.TemporaryDirectory()
temp_path = Path(temp_dir.name)
save_anndata(temp_path / "input.h5ad", [adata])

model = SpiceMix(
    dataset_path=temp_path / "input.h5ad",
    random_state=args.seed,
    initial_context=torch_context,
    torch_context=torch_context,
    **config
)

# Run
model.train(num_preiterations=num_preiterations, num_iterations=num_iterations)

tl.preprocess_embeddings(model, normalized_key="normalized_X")
tl.leiden(model, joint=True, target_clusters=args.n_clusters, use_rep="normalized_X")

# TODO: add optional smoothing step
temp_dir.cleanup()
label_df = model.datasets[0].obs[["leiden"]]

## Write output


label_df.columns = ["label"]
label_df.to_csv(label_file, sep="\t", index_label="")