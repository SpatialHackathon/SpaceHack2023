#!/usr/bin/env python

# Author_and_contribution: 
# Niklas Mueller-Boetticher; created template
# Jieran Sun; Implement SEDR method

import argparse

# TODO adjust description
parser = argparse.ArgumentParser(
    description="SEDR : Unsupervised spatially embedded deep representation of spatial transcriptomics. See https://doi.org/10.1101/2021.06.15.448542 for details"
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
    from PIL import Image

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

    if args.image is not None:
        adata.uns["image"] = np.array(Image.open(args.img))
    else:
        adata.uns["image"] = None

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
import sys
from pathlib import Path

import os
import warnings
warnings.filterwarnings('ignore')

# Add the parent directory of the current file to sys.path
method_dir = Path(__file__).resolve().parent.parent  # Navigate two levels up
sys.path.append(str(method_dir))

from search_res import binary_search

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

device = 'cuda:0' if torch.cuda.is_available() else 'cpu'

# Load data
adata = get_anndata(args)
adata.var_names_make_unique()

# Add preprocessing steps of the method
if not isinstance(adata.X, np.ndarray):
    adata.layers['count'] = adata.X.toarray()
else:
    adata.layers['count'] = adata.X

sc.pp.normalize_total(adata, target_sum=1e6)

if adata.n_vars > 2000 and config["HVG"] is True:
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", layer='count', n_top_genes=2000)
    adata = adata[:, adata.var['highly_variable'] == True]

sc.pp.scale(adata)

from sklearn.decomposition import PCA  # sklearn PCA is used because PCA in scanpy is not stable.
adata_X = PCA(n_components=200, random_state=seed).fit_transform(adata.X)
adata.obsm['X_pca'] = adata_X

graph_dict = SEDR.graph_construction(adata, config["n"])

# Training SEDR
# device: using cpu or gpu (if avaliable)
# using_dec: boolean, whether to use the unsupervised deep embedded clustering (DEC) method to improve clustering results 
sedr_net = SEDR.Sedr(adata.obsm['X_pca'],
                     graph_dict,
                     mode='clustering',
                     device=device)

if config["using_dec"]:
    sedr_net.train_with_dec(N=1)
else:
    sedr_net.train_without_dec(N=1)
sedr_feat, _, _, _ = sedr_net.process()
# latent embedding
adata.obsm['SEDR'] = sedr_feat

# Clustering 
if config['cluster_method'] == "mclust":
    adata = mclust_R(adata, 
                            n_clusters = n_clusters, 
                            use_rep='SEDR', 
                            key_added='SEDR', 
                            random_seed=seed
                            )
    label_df = adata.obs[["SEDR"]]
    
else:
    sc.pp.neighbors(adata, use_rep="SEDR")
    label_df = binary_search(adata, n_clust_target=n_clusters, method=config['cluster_method'], seed = seed)

# Output dataframes

embedding_df = pd.DataFrame(adata.obsm['SEDR'])
embedding_df.index = adata.obs_names

## Write output
out_dir.mkdir(parents=True, exist_ok=True)
label_df.columns = ["label"]
label_df.to_csv(label_file, sep="\t", index_label="")

if embedding_df is not None:
    embedding_df.to_csv(embedding_file, sep="\t", index_label="")