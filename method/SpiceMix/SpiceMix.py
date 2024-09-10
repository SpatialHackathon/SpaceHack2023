#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template, minor revisions
# Author_and_contribution: Shahul Alam; implemented SpiceMix

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
parser.add_argument(
    "-f", "--features", help="Path to features (as tsv).", required=True
)
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
    "--dim_red",
    help="Reduced dimensionality representation (e.g. PCA).",
    required=False,
)
parser.add_argument("--image", help="Path to H&E staining.", required=False)
parser.add_argument(
    "--n_clusters", help="Number of clusters to return.", required=True, type=int
)
parser.add_argument(
    "--n_genes", help="Number of genes to use.", required=False, type=int
)
parser.add_argument(
    "--technology",
    help="The technology of the dataset (Visium, ST, imaging-based).",
    required=True,
)
parser.add_argument(
    "--seed", help="Seed to use for random operations.", required=True, type=int
)
parser.add_argument(
    "--config",
    help="Optional config file (json) used to pass additional parameters.",
    required=False,
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

if args.n_genes is not None:
    n_genes = args.n_genes
elif preprocess_parameters:
    n_genes = preprocess_parameters["n_genes"]
    
def get_anndata(args):
    """Convert data input into SpiceMix anndata format."""

    counts_matrix = mmread(args.matrix)
    if issparse(counts_matrix):
        counts_matrix = counts_matrix.tocsr()

    observations = pd.read_csv(args.observations, delimiter="\t", index_col=0)
    features = pd.read_table(args.features, index_col=0)

    coordinates = (
        pd.read_table(args.coordinates, index_col=0)
        .loc[observations.index, :]
        .iloc[:, :2]
        .to_numpy()
    )

    # adjacency_matrix = mmread(args.neighbors).T.tocsr()

    adata = ad.AnnData(
        X=counts_matrix,
        obs=observations,
        obsm={"spatial": coordinates} #,
        # obsp={"spatial_connectivities": adjacency_matrix},
    )

    # Filter by selected samples
    if "selected" in adata.obs.columns:
        adata = adata[observations["selected"].astype(bool), :]

    if preprocess_parameters:
        # del adata.obsp["spatial_connectivities"]
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=n_genes)

        adata = PopariDataset(adata[:, adata.var["highly_variable"]], "processed")
        adata.compute_spatial_neighbors()

    else:
        if "selected" in adata.var.columns:
            adata = adata[:, features["selected"].astype(bool)]

        adjacency_matrix = adata.obsp["spatial_connectivities"]

        # Symmetrize matrix
        transpose_condition = adjacency_matrix.T > adjacency_matrix
        adjacency_matrix = (
            adjacency_matrix
            + adjacency_matrix.T.multiply(transpose_condition)
            - adjacency_matrix.multiply(transpose_condition)
        )

        num_cells = adjacency_matrix.shape[0]
        adjacency_list = [[] for _ in range(num_cells)]
        for x, y in zip(*adjacency_matrix.nonzero()):
            adjacency_list[x].append(y)

        adata.obsp["adjacency_matrix"] = adjacency_matrix
        adata.obs["adjacency_list"] = adjacency_list

        del adata.obsp["spatial_connectivities"]

        adata = PopariDataset(adata, "processed")

    return adata


adata = get_anndata(args)

num_preiterations = config.pop("num_preiterations", 5)
num_iterations = config.pop("num_iterations", 200)

device = config.pop("device", "cpu") if torch.cuda.is_available() else 'cpu'
dtype = config.pop("dtype", "float32")

if dtype == "float32":
    dtype = torch.float32
elif dtype == "float64":
    dtype = torch.float64

torch_context = dict(
    device=device,
    dtype=dtype,
)

# ignore "source" in configs
keys_to_ignore = ["source", "device"]
filtered_config = {key: value for key, value in config.items() if key not in keys_to_ignore}

with tempfile.TemporaryDirectory() as temp_dir:
    temp_path = Path(temp_dir)
    save_anndata(temp_path / "input.h5ad", [adata])

    model = SpiceMix(
        dataset_path=temp_path / "input.h5ad",
        random_state=args.seed,
        initial_context=torch_context,
        torch_context=torch_context,
        **filtered_config,
    )

    # Run
    model.train(num_preiterations=num_preiterations, num_iterations=num_iterations)

    tl.preprocess_embeddings(model, normalized_key="normalized_X")
    tl.leiden(
        model, joint=True, target_clusters=args.n_clusters, use_rep="normalized_X"
    )
    # TODO: add optional smoothing step
    label_df = model.datasets[0].obs[["leiden"]]

## Write output
label_df.columns = ["label"]
label_df.to_csv(label_file, sep="\t", index_label="")
