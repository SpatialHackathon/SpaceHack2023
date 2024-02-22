#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Florian Heyl (heylf); created code

import argparse
import os
import json

###### Template arguments ##################################################
parser = argparse.ArgumentParser(description="conST: an interpretable multi-modal contrastive learning \
                                 framework for spatial transcriptomics")

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

args = parser.parse_args()

from pathlib import Path

out_dir = Path(args.out_dir)

# Output files
label_file = out_dir / "domains.tsv"
embedding_file = out_dir / "embedding.tsv"

coord_file = args.coordinates
feature_file = args.features
observation_file = args.observations

if args.neighbors is not None:
    neighbors_file = args.neighbors
if args.matrix is not None:
    matrix_file = args.matrix
if args.dim_red is not None:
    dimred_file = args.dim_red
if args.image is not None:
    image_file = args.image
if args.config is not None:
    with open(args.config, "r") as f:
        config = json.load(f)

n_clusters = args.n_clusters
technology = args.technology
seed = args.seed

# ... or AnnData if you want
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


adata = get_anndata(args)
adata.var_names_make_unique()


##############
### conST ####
########################################################################################################################

if not os.path.exists(f'{out_dir}'):
    os.makedirs(f'{out_dir}')

params = args

# Add config to params
for key, value in config.items():
    setattr(params, key, value)

# Tool imports
import random
import torch
import numpy as np
import pandas as pd
import tempfile
import anndata
import matplotlib.pyplot as plt
import scanpy as sc
import warnings
from sklearn import metrics
warnings.filterwarnings('ignore')

random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
print('Using device: ' + device)
params.device = device

os.environ['PYTHONHASHSEED'] = str(seed)
torch.backends.cudnn.benchmark = False
torch.backends.cudnn.deterministic = True

# TODO this is what conST uses for data preprocessing
def adata_preprocess(i_adata, min_cells=3, pca_n_comps=300):
    print('===== Preprocessing Data ')
    sc.pp.filter_genes(i_adata, min_cells=min_cells)
    adata_X = sc.pp.normalize_total(i_adata, target_sum=1, exclude_highly_expressed=True, inplace=False)['X']
    adata_X = sc.pp.scale(adata_X)
    adata_X = sc.pp.pca(adata_X, n_comps=pca_n_comps)
    return adata_X

# TODO should data processing be changed to the processing of scanpy?
# scanpy starts here
# if not args.dim_red:
#     sc.pp.normalize_total(adata)
#     sc.pp.log1p(adata)
#     sc.tl.pca(adata)

# Work in a temprary folder
with tempfile.TemporaryDirectory() as tmpdir:
    gitdir = f"{str(tmpdir)}/conST"

    # Clone the repository to the specific commit
    os.system(
    f"""
    git clone https://github.com/ys-zong/conST {gitdir}
    cd {gitdir} 
    git reset --hard a32d747
    """)

    # Further tool imports
    from src.graph_func import graph_construction
    from src.utils_func import res_search_fixed_clus, plot_clustering
    from src.training import conST_training

    # Preprocessing
    adata_X = adata_preprocess(adata, min_cells=5, pca_n_comps=params.cell_feat_dim)

    # Graph construction
    graph_dict = graph_construction(adata.obsm['spatial'], adata.shape[0], params)

    params.cell_num = adata.shape[0]

    # Use image data if provided
    if params.use_img:
        # TODO please check if this is the correct way to use the image data
        img_transformed = adata.uns["image"]
        img_transformed = (img_transformed - img_transformed.mean()) / img_transformed.std() * adata_X.std() + adata_X.mean()
        conST_net = conST_training(adata_X, graph_dict, params, n_clusters, img_transformed)
    else:
        conST_net = conST_training(adata_X, graph_dict, params, n_clusters)

    # Training phase
    conST_net.pretraining()
    conST_net.major_training()

    # Get embeddings
    conST_embedding = conST_net.get_embedding()

    # Clustering
    adata_conST = anndata.AnnData(conST_embedding)
    adata_conST.obsm['spatial'] = adata.obsm['spatial']

    sc.pp.neighbors(adata_conST, n_neighbors=params.eval_graph_n)

    eval_resolution = res_search_fixed_clus(adata_conST, n_clusters)
    print(eval_resolution)

    sc.tl.leiden(adata_conST, key_added="conST_leiden", resolution=eval_resolution)

    # Plot leiden clsuters without refinement
    plot_clustering(adata_conST, "conST_leiden", savepath = f'{out_dir}/conST_leiden_plot.jpg')

    # Visium or ST data refinement
    index = np.arange(start=0, stop=adata_X.shape[0]).tolist()
    index = [str(x) for x in index]

    def refine(sample_id, pred, dis, shape="hexagon"):
        refined_pred=[]
        pred=pd.DataFrame({"pred": pred}, index=sample_id)
        dis_df=pd.DataFrame(dis, index=sample_id, columns=sample_id)
        if shape=="hexagon":
            num_nbs=6 
        elif shape=="square":
            num_nbs=4
        else:
            print("Shape not recongized, shape='hexagon' for Visium data, 'square' for ST data.")
        for i in range(len(sample_id)):
            index=sample_id[i]
            dis_tmp=dis_df.loc[index, :].sort_values(ascending=False)
            nbs=dis_tmp[0:num_nbs+1]
            nbs_pred=pred.loc[nbs.index, "pred"]
            self_pred=pred.loc[index, "pred"]
            v_c=nbs_pred.value_counts()
            if (v_c.loc[self_pred]<num_nbs/2) and (np.max(v_c)>num_nbs/2):
                refined_pred.append(v_c.idxmax())
            else:           
                refined_pred.append(self_pred)
        return refined_pred

    dis = graph_dict['adj_norm'].to_dense().numpy() + np.eye(graph_dict['adj_norm'].shape[0])
    refine = refine(sample_id = index, pred = adata_conST.obs['leiden'].tolist(), dis=dis, shape=params.shape)
    adata_conST.obs['refine'] = refine

    # Plot leiiden clusters with refinement
    plot_clustering(adata_conST, 'refine', savepath = f'{out_dir}/conST_leiden_plot_refined.jpg')

    # Select labels with or without refinement
    cluster_key = 'conST_leiden'
    if (params.refinement):
        cluster_key = 'refine'

    ## Write output
    out_dir.mkdir(parents=True, exist_ok=True)

    label_df = pd.DataFrame({'labels': adata_conST.obs[cluster_key].tolist()}, index=adata.obs.index.tolist())
    label_df.to_csv(label_file, sep="\t", index_label="")

    conST_embedding = pd.DataFrame(conST_embedding)
    conST_embedding.index = adata.obs.index.tolist()

    if conST_embedding is not None:
        conST_embedding.to_csv(embedding_file, sep="\t", index_label="")
