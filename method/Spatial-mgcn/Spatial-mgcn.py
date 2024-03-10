#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created script
# Author_and_contribution: Kirti Biharie; Added Spatial-MGCN

import argparse

CALC_ARI = False

parser = argparse.ArgumentParser(
    description="""Spatial-MGCN (https://doi.org/10.1093/bib/bbad262)"""
)

parser.add_argument(
    "-c", "--coordinates", help="Path to coordinates (as tsv).", required=True
)
parser.add_argument(
    "-m", "--matrix", help="Path to (transformed) counts (as mtx).", required=False
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
    "--technology",
    help="The technology of the dataset (Visium, ST, ...).",
    required=True,
)
parser.add_argument(
    "--seed", help="Seed to use for random operations.", required=True, type=int
)
parser.add_argument(
    "--config",
    help="Optional config file used to pass additional parameters.",
    required=False,
)

# Uncomment to calculate ARI every epoch as in original implementation
if CALC_ARI:
    parser.add_argument(
        "-g", "--groundtruth",
        help="Groundtruth.",
        required=False,
    )

args = parser.parse_args()

from pathlib import Path

out_dir = Path(args.out_dir)

# Output files
label_file = out_dir / "domains.tsv"
embedding_file = out_dir / "embedding.tsv"

n_clusters = args.n_clusters
technology = args.technology
seed = args.seed

## Your code goes here
import json

with open(args.config, "r") as f:
    config = json.load(f)

import random
import numpy as np
import pandas as pd
import torch
import tempfile
import os
import sys
import scanpy as sc
import scipy as sp
import torch.optim
import sklearn.cluster
import sklearn.metrics
import tqdm

def get_anndata(args):
    import anndata as ad
    
    from PIL import Image
    import scipy.io

    X = sp.io.mmread(args.matrix)
    if sp.sparse.issparse(X):
        X = X.tocsr()

    observations = pd.read_table(args.observations, index_col=0)
    features = pd.read_table(args.features, index_col=0)

    # Filter
    if "selected" in observations.columns:
        X = X[observations["selected"].to_numpy().nonzero()[0], :]
        observations = observations.loc[lambda df: df["selected"]]
    if "selected" in features.columns:
        X = X[:, features["selected"].to_numpy().nonzero()[0]]
        features = features.loc[lambda df: df["selected"]]

    coordinates = (
        pd.read_table(args.coordinates, index_col=0)
        .loc[observations.index, :]
        .to_numpy()
    )

    adata = ad.AnnData( # Rename spatial_pixel to spatial for Spatial-mgcn
        X=X, obs=observations, var=features, obsm={"spatial": coordinates}
    )

    if args.image is not None:
        adata.uns["image"] = np.array(Image.open(args.image))
    else:
        adata.uns["image"] = None

    return adata


adata = get_anndata(args)

if CALC_ARI:
    labels = pd.read_table(args.groundtruth, index_col=0)
    adata = adata[adata.obs_names.isin(labels.index)]
    labels = labels.loc[adata.obs_names]

# Set seed
random.seed(seed)
torch.manual_seed(seed)
np.random.seed(seed)
torch.cuda.manual_seed(seed)

# Work in a temprary folder
with tempfile.TemporaryDirectory() as tmpdir:
    gitdir = f"{str(tmpdir)}/Spatial-MGCN"

    # Clone the repository to the specific commit
    os.system(
    f"""
    git clone https://github.com/cs-wangbo/Spatial-MGCN.git {gitdir}
    cd {gitdir} 
    git reset --hard cf4412d
    """)

    # Set working directory as Spatial-MGCN directory
    sys.path.append(f"{gitdir}/Spatial-MGCN")
    import utils
    import models

    # Normalize data: min_cells, calculate HVG and scale
    # sc.pp.filter_genes(adata, min_cells=100)
    # sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=config["fdim"])
    # adata = adata[:, adata.var['highly_variable']].copy()
    # adata.X = adata.X / np.sum(adata.X, axis=1).reshape(-1, 1) * 10000
    # adata.X = sp.sparse.csr_matrix(adata.X)
    # sc.pp.scale(adata, zero_center=False, max_value=10)

    config["fdim"] = len(adata.var_names)

    # Calculate graphs
    fadj = utils.features_construct_graph(adata.X, k=config["k"]) 
    sadj, graph_nei, graph_neg = utils.spatial_construct_graph1(adata, radius=config["radius"])

    adata.obsm["fadj"] = fadj
    adata.obsm["sadj"] = sadj
    adata.obsm["graph_nei"] = graph_nei.numpy()
    adata.obsm["graph_neg"] = graph_neg.numpy()

    features = torch.FloatTensor(adata.X.todense())
    
    nfadj = utils.normalize_sparse_matrix(fadj + sp.eye(fadj.shape[0]))
    nfadj = sp.sparse.csr_matrix(nfadj)
    nfadj = utils.sparse_mx_to_torch_sparse_tensor(nfadj)
    
    nsadj = utils.normalize_sparse_matrix(sadj + sp.eye(sadj.shape[0]))
    nsadj = sp.sparse.csr_matrix(nsadj)
    nsadj = utils.sparse_mx_to_torch_sparse_tensor(nsadj)
    
    graph_nei = torch.LongTensor(adata.obsm['graph_nei'])
    graph_neg = torch.LongTensor(adata.obsm['graph_neg'])

    # Create model
    cuda = torch.cuda.is_available()

    if cuda:
            features = features.cuda()
            nsadj = nsadj.cuda()
            nfadj = nfadj.cuda()
            graph_nei = graph_nei.cuda()
            graph_neg = graph_neg.cuda()
    
    model = models.Spatial_MGCN(nfeat=config["fdim"],
                             nhid1=config["nhid1"],
                             nhid2=config["nhid2"],
                             dropout=config["dropout"])
    
    if cuda:
        model.cuda()
    
    optimizer = torch.optim.Adam(model.parameters(), lr=config["lr"], weight_decay=config["weight_decay"])

    # Train model
    epoch_max = 0
    ari_max = 0
    idx_max = []
    mean_max = []
    emb_max = []

    for epoch in tqdm.tqdm(range(config["epochs"])):
        model.train()
        optimizer.zero_grad()
        com1, com2, emb, pi, disp, mean = model(features, nsadj, nfadj)
        zinb_loss = utils.ZINB(pi, theta=disp, ridge_lambda=0).loss(features, mean, mean=True)
        reg_loss = utils.regularization_loss(emb, graph_nei, graph_neg)
        con_loss = utils.consistency_loss(com1, com2)
        total_loss = config["alpha"] * zinb_loss + config["beta"] * con_loss + config["gamma"] * reg_loss
        emb = pd.DataFrame(emb.cpu().detach().numpy()).fillna(0).values
        total_loss.backward()
        optimizer.step()
        
        kmeans = sklearn.cluster.KMeans(n_clusters=args.n_clusters, n_init=10).fit(emb)
        idx = kmeans.labels_

        if CALC_ARI:
            ari_res = sklearn.metrics.adjusted_rand_score(labels.to_numpy()[:,0], idx)
            if ari_res > ari_max:
                ari_max = ari_res
                idx_max = idx
                emb_max = emb
        else:
            idx_max = idx
            emb_max = emb
    
    # Write output
    emb_df = pd.DataFrame(emb_max, index=adata.obs_names)
    label_df = pd.DataFrame(idx_max, index=adata.obs_names, columns=["label"])

    emb_df.to_csv(embedding_file, sep="\t", index_label="")
    label_df.to_csv(label_file, sep="\t", index_label="")
