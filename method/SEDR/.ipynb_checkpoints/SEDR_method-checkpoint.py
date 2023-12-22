#!/usr/bin/env python

# Author_and_contribution: 
# Niklas Mueller-Boetticher; created template
# Jieran Sun; Implement SEDR method

import argparse

# TODO adjust description
parser = argparse.ArgumentParser(
    description="SEDR – Unsupervised spatially embedded deep representation of spatial transcriptomics. See https://doi.org/10.1101/2021.06.15.448542 for details"
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

## Session for code
args = parser.parse_args()

# anndata input
def get_anndata(args):
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy as sp

    observations = pd.read_table(args.observations, index_col=0)
    features = pd.read_table(args.features, index_col=0)

    coordinates = (
        pd.read_table(args.coordinates, index_col=0)
        .loc[observations.index, :]
        .to_numpy()
    )

    adata = ad.AnnData(obs=observations, var=features, obsm={"spatial": coordinates})

    if args.matrix is not None:
        X = sp.io.mmread(args.matrix)
        if sp.sparse.issparse(X):
            X = X.tocsr()
        adata.X = X

    if args.neighbors is not None:
        adata.obsp["spatial_connectivities"] = sp.io.mmread(args.neighbors).T.tocsr()

    # Filter by selected samples/features
    if "selected" in adata.obs.columns:
        adata = adata[observations["selected"].astype(bool), :]
    if "selected" in adata.var.columns:
        adata = adata[:, features["selected"].astype(bool)]

    if args.dim_red is not None:
        adata.obsm["reduced_dimensions"] = (
            pd.read_table(args.dim_red, index_col=0).loc[adata.obs_names].to_numpy()
        )

    return adata

## Cluster code is SEDR had some bugs. Modified below
def res_search_fixed_clus_leiden(adata, n_clusters, increment=0.01, random_seed=2023):
    import scanpy as sc
    import pandas as pd
    import numpy as np
    
    for res in np.arange(0.2, 2, increment):
        sc.tl.leiden(adata, random_state=random_seed, resolution=res)
        if len(adata.obs['leiden'].unique()) > n_clusters:
            break
    return res-increment


def leiden(adata, n_clusters, use_rep='SEDR', key_added='SEDR', random_seed=2023):
    import scanpy as sc
    import pandas as pd
    import numpy as np
    
    sc.pp.neighbors(adata, use_rep=use_rep)
    res = res_search_fixed_clus_leiden(adata, n_clusters, increment=0.01, random_seed=random_seed)
    sc.tl.leiden(adata, random_state=random_seed, resolution=res)

    adata.obs[key_added] = adata.obs['leiden']
    adata.obs[key_added] = adata.obs[key_added].astype('int')
    adata.obs[key_added] = adata.obs[key_added].astype('category')

    return adata


def res_search_fixed_clus_louvain(adata, n_clusters, increment=0.01, random_seed=2023):
    import scanpy as sc
    import pandas as pd
    import numpy as np
    
    for res in np.arange(0.2, 2, increment):
        sc.tl.louvain(adata, random_state=random_seed, resolution=res)
        if len(adata.obs['louvain'].unique()) > n_clusters:
            break
    return res-increment

def louvain(adata, n_clusters, use_rep='SEDR', key_added='SEDR', random_seed=2023):
    import scanpy as sc
    import pandas as pd
    import numpy as np
    
    sc.pp.neighbors(adata, use_rep=use_rep)
    res = res_search_fixed_clus_louvain(adata, n_clusters, increment=0.01, random_seed=random_seed)
    sc.tl.louvain(adata, random_state=random_seed, resolution=res)

    adata.obs[key_added] = adata.obs['louvain']
    adata.obs[key_added] = adata.obs[key_added].astype('int')
    adata.obs[key_added] = adata.obs[key_added].astype('category')

    return adata

# modified SDER mclust function cuz in their they define their own R home pathway
def mclust_R(adata, n_clusters, use_rep='SEDR', key_added='SEDR', modelNames= "EEE", random_seed=2023):
    """\
    Clustering using the mclust algorithm.
    The parameters are the same as those in the R package mclust.
    """
    import numpy as np

    np.random.seed(random_seed)
    import rpy2.robjects as robjects
    robjects.r.library("mclust")

    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()
    r_random_seed = robjects.r['set.seed']
    r_random_seed(random_seed)
    rmclust = robjects.r['Mclust']

    res = rmclust(rpy2.robjects.numpy2ri.numpy2rpy(adata.obsm[use_rep]), n_clusters, modelNames)
    mclust_res = np.array(res[-2])

    adata.obs[key_added] = mclust_res
    adata.obs[key_added] = adata.obs[key_added].astype('int')
    adata.obs[key_added] = adata.obs[key_added].astype('category')

    return adata
    
# Import packages for SEDR
import scanpy as sc
import numpy as np
import pandas as pd
import scipy as sp
import torch
import SEDR

import os
import warnings
warnings.filterwarnings('ignore')

# Def config
import json
with open(args.config, "r") as f:
    config = json.load(f)

# Set up output files
from pathlib import Path
out_dir = Path(args.out_dir)
label_file = out_dir / "domains.tsv"
embedding_file = out_dir / "embedding.tsv"

# Def vars
n_clusters = args.n_clusters
technology = args.technology
seed = args.seed

# Set SEED for SEDR
import random
random.seed(seed)
SEDR.fix_seed(seed)

# Load data
adata = get_anndata(args)
adata.var_names_make_unique()

# if dim-red is not provided, use default PCA dimRed
if "reduced_dimensions" not in adata.obsm_keys():
    from sklearn.decomposition import PCA 
    adata_X = PCA(n_components=200, 
                  random_state=seed).fit_transform(adata.X)
    adata.obsm['reduced_dimensions'] = adata_X

# Constructing neighborhood graphs if neighbors not provided
if "spatial_connectivities" in adata.obsp.keys():
    # import intermediate functions
    from SEDR.graph_func import preprocess_graph
    
    # Copy from source code in order for customization
    adj_m1 = adata.obsp["spatial_connectivities"]
    adj_m1 = sp.sparse.coo_matrix(adj_m1)

    # Store original adjacency matrix (without diagonal entries) for later
    adj_m1 = adj_m1 - sp.sparse.dia_matrix((adj_m1.diagonal()[np.newaxis, :], [0]), shape=adj_m1.shape)
    adj_m1.eliminate_zeros()

    # Some preprocessing
    adj_norm_m1 = preprocess_graph(adj_m1)
    adj_m1 = adj_m1 + sp.sparse.eye(adj_m1.shape[0])

    adj_m1 = adj_m1.tocoo()
    shape = adj_m1.shape
    values = adj_m1.data
    indices = np.stack([adj_m1.row, adj_m1.col])
    adj_label_m1 = torch.sparse_coo_tensor(indices, values, shape)

    norm_m1 = adj_m1.shape[0] * adj_m1.shape[0] / float((adj_m1.shape[0] * adj_m1.shape[0] - adj_m1.sum()) * 2)

    graph_dict = {
        "adj_norm": adj_norm_m1,
        "adj_label": adj_label_m1.coalesce(),
        "norm_value": norm_m1
    }
    
else: 
    graph_dict = SEDR.graph_construction(adata, 
                                         n=config['graph_dmax'], 
                                         dmax=config['graph_n'], 
                                         mode=config['graph_mode'])
    

# Training SEDR
# device: using cpu or gpu (if avaliable)
# using_dec: boolean, whether to use the unsupervised deep embedded clustering (DEC) method to improve clustering results 
sedr_net = SEDR.Sedr(adata.obsm['reduced_dimensions'], 
                     graph_dict, 
                     mode='clustering', 
                     device=config["device"])

if config["using_dec"]:
    sedr_net.train_with_dec(N=1)
else:
    sedr_net.train_without_dec(N=1)
sedr_feat, _, _, _ = sedr_net.process()
# latent embedding
adata.obsm['SEDR'] = sedr_feat

# Clustering 
match config['cluster_method']:
    case "mclust":
        adata = mclust_R(adata, 
                              n_clusters = n_clusters, 
                              use_rep='SEDR', 
                              key_added='SEDR', 
                              random_seed=seed
                             )
    case "louvain":
         adata = louvain(adata, 
                              n_clusters = n_clusters, 
                              use_rep='SEDR', 
                              key_added='SEDR', 
                              random_seed=seed
                             )
    case "leiden":
         adata = leiden(adata, 
                              n_clusters = n_clusters, 
                              use_rep='SEDR', 
                              key_added='SEDR', 
                              random_seed=seed
                             )
        
# Output dataframes
label_df = adata.obs[["SEDR"]]
embedding_df = pd.DataFrame(adata.obsm['SEDR'])
embedding_df.index = adata.obs_names

## Write output
out_dir.mkdir(parents=True, exist_ok=True)
label_df.columns = ["label"]
label_df.to_csv(label_file, sep="\t", index_label="")

if embedding_df is not None:
    embedding_df.to_csv(embedding_file, sep="\t", index_label="")
