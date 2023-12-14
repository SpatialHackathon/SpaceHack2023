#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Louis B. Kuemmerle; wrapped cellcharter's purity metric
# Author_and_contribution: Marco Varrone; gave input on setting hyperparameters and heuristics that generalise well and helped with debugging

# At the time of development (12th Dec 2023) cellcharter doesn't have release tags. The code used for development can
# be accessed via (commit 2b15549968559d9023d1fc92f9798a8c38dffcac):
# 
# git clone https://github.com/CSOgroup/cellcharter.git
# cd cellcharter
# git checkout -b [new-branch-name] 2b15549968559d9023d1fc92f9798a8c38dffcac

import argparse
import json
import numpy as np
import pandas as pd
import anndata as ad
import squidpy as sq
import cellcharter as cc

# TODO: template arguments are strict. Need to readjust here.
parser = argparse.ArgumentParser(description="Calculate cellcharter's purity metric that measures if cells assigned to other connected components are located within the boundaries of given components.")
parser.add_argument("--coordinates", help="Path to coordinates (as tsv).", required=True)
parser.add_argument("-l", "--labels", help="Labels from domain clustering.", required=True)
parser.add_argument(
    "-c", 
    "--config",
    help="Config file (json) used to pass additional parameters.",
    required=True,
)
parser.add_argument("-o", "--out_file", help="Output file.", required=True)
parser.add_argument("--fig_path", help="File path for plot to check if components look reasonable.", required=False, default=False) #TODO: Delete

SUPPORTED_TECHNOLOGIES = ["CODEX", "CosMx", "MERFISH", "IMC", "Visium", "RNA + ATAC"]
# for each technology we need to specify how to compute the graph connectivity and how to restrict connected components


###### helper functions ######

from anndata import AnnData
from scipy.sparse import csr_matrix
from squidpy._constants._constants import CoordType, Transform
from squidpy._constants._pkg_constants import Key
from squidpy.gr._utils import _assert_connectivity_key

def remove_long_links(
    adata: AnnData,
    distance_percentile: float | None = 99.0,
    mean_distance_factor: float | None = 4.0,
    connectivity_key: str | None = None,
    distances_key: str | None = None,
    neighs_key: str | None = None,
    copy: bool = False,
) -> tuple[csr_matrix, csr_matrix] | None:
    """
    Remove links between cells at a distance bigger than a certain percentile or mean distance of all positive distances.

    It is designed for data with generic coordinates.
    
    Note: This function is copied from cellcharter.graph.remove_long_links and extended to allow for a mean distance 
    factor.

    Parameters
    ----------
    %(adata)s

    distance_percentile
        Percentile of the distances between cells over which links are trimmed after the network is built.
    mean_distance_factor
        Factor to multiply the mean distance between cells by to get the threshold for trimming links. If
        ``distance_percentile`` is not ``None``, this parameter is ignored.
    %(conn_key)s

    distances_key
        Key in :attr:`anndata.AnnData.obsp` where spatial distances are stored.
        Default is: :attr:`anndata.AnnData.obsp` ``['{{Key.obsp.spatial_dist()}}']``.
    neighs_key
        Key in :attr:`anndata.AnnData.uns` where the parameters from gr.spatial_neighbors are stored.
        Default is: :attr:`anndata.AnnData.uns` ``['{{Key.uns.spatial_neighs()}}']``.

    %(copy)s

    Returns
    -------
    If ``copy = True``, returns a :class:`tuple` with the new spatial connectivities and distances matrices.

    Otherwise, modifies the ``adata`` with the following keys:
        - :attr:`anndata.AnnData.obsp` ``['{{connectivity_key}}']`` - the new spatial connectivities.
        - :attr:`anndata.AnnData.obsp` ``['{{distances_key}}']`` - the new spatial distances.
        - :attr:`anndata.AnnData.uns`  ``['{{neighs_key}}']`` - :class:`dict` containing parameters.
    """
    connectivity_key = Key.obsp.spatial_conn(connectivity_key)
    distances_key = Key.obsp.spatial_dist(distances_key)
    neighs_key = Key.uns.spatial_neighs(neighs_key)
    _assert_connectivity_key(adata, connectivity_key)
    _assert_connectivity_key(adata, distances_key)

    conns, dists = adata.obsp[connectivity_key], adata.obsp[distances_key]

    if copy:
        conns, dists = conns.copy(), dists.copy()

    if distance_percentile is not None:
        threshold = np.percentile(np.array(dists[dists != 0]).squeeze(), distance_percentile)
    else:
        threshold = np.mean(np.array(dists[dists != 0]).squeeze()) * mean_distance_factor
    conns[dists > threshold] = 0
    dists[dists > threshold] = 0

    conns.eliminate_zeros()
    dists.eliminate_zeros()

    if copy:
        return conns, dists
    else:
        adata.uns[neighs_key]["params"]["radius"] = threshold

def get_alpha_start(adata, factor = 0.05):
    """ Compute alpha start for boundary computation.
    
    Alpha start is computed based on a heuristic. We take the minimum length of the diagonals of the bounding boxes
    over connected components and multiply it with a factor.
    
    The alpha start shouldn't be too large, otherwise alpha shapes get too close to the convex hull. Alpha start should
    also not be too small to reduce computation time when searching for the optimal alpha.
    
    """
    distances = []
    for component in adata.obs["component"].dropna().unique():
        coords = adata[adata.obs["component"]==component].obsm["spatial"]
        x1 = np.array([np.min(coords[:,0]), np.min(coords[:,1])])
        x2 = np.array([np.max(coords[:,0]), np.max(coords[:,1])])
        distances.append(np.linalg.norm(x1-x2))
        
    return np.min(distances) * factor


###### functions for debug plots ######
# It's important that reasonable connected components are identified. This needs to be checked on different 
# technologies, tissues, and spatial domains. Potentially the adjacency graph parameters need to be adjusted for better
# generalisation.

def get_n_colors(n):
    from random import shuffle
    
    # Generate a slightly larger list of colors to account for the removal of black and white
    cmap = plt.cm.get_cmap('nipy_spectral', n + 2)  # Adjusting for potential removal

    # Generate colors and convert to hex
    colors = cmap(np.linspace(0, 1, n + 2))
    hex_colors = ['#%02x%02x%02x' % (int(r*255), int(g*255), int(b*255)) for r, g, b, _ in colors]

    # Remove black and white from the list
    hex_colors = [color for color in hex_colors if color not in ['#000000', '#ffffff']]

    # Ensure the list has exactly n colors by taking the first n colors
    shuffle(hex_colors)
    return hex_colors[:n]

import matplotlib.pyplot as plt

def plot_components(adata, save=False):
    """ Plot components and their boundaries.
    
    """    
    n = len(adata.obs["spatial_cluster"].dropna().unique())
    colors = get_n_colors(n)
    component_to_spatial_cluster = dict(zip(adata.obs["component"], adata.obs["spatial_cluster"]))
    
    # color_code strings 
    df = adata.obs
    color_code = df.loc[~df['spatial_cluster'].isnull(),'spatial_cluster'].unique()
    label_to_color = dict(zip(color_code,colors))
    df['color'] = 'lightgray'
    df.loc[~df['spatial_cluster'].isnull(),'color'] = df.loc[~df['spatial_cluster'].isnull(),'spatial_cluster'].map(label_to_color)
    
    fig = plt.figure(figsize=(10,10))
    for i, polygon1 in adata.uns[f"shape_component"]["boundary"].items():
        plt.plot(*polygon1.exterior.xy, color=label_to_color[component_to_spatial_cluster[i]])
        
    plt.scatter(adata.obsm["spatial"][:,0],adata.obsm["spatial"][:,1],c=adata.obs['color'], s=1)

    if save:
        fig.savefig(save)


###### main function ######

def script():
    
    args = parser.parse_args()
    
    # Set paths and variables
    coord_file = args.coordinates
    label_file = args.labels
    fig_path = args.fig_path #TODO: Delete?
    
    if args.config is not None:
        config_file = args.config
        with open(config_file) as f:
            config = json.load(f)
        sample_file = config["sample_file"] if "sample_file" in config else None
        technology = config["technology"] if "technology" in config else None
        if technology not in SUPPORTED_TECHNOLOGIES:
            raise ValueError(f"Technology {technology} not supported. Supported technologies are {SUPPORTED_TECHNOLOGIES}.")
    
    # Init adata
    obs = pd.DataFrame(
        data={
            "spatial_cluster":pd.read_table(label_file,index_col=0).iloc[:,0].values,
            "sample": "sample1" if sample_file is None else pd.read_table(sample_file,index_col=0).iloc[:,0].values,
            }
    )
    obs.index = obs.index.astype(str)
    obsm = {"spatial": pd.read_table(coord_file,index_col=0).values[:,:2]}
    adata = ad.AnnData(X=None, obs=obs, obsm=obsm) 
    
    # Compute graph connectivity
    """
    From cellcharter paper:
    Spatial network construction
    Depending on the technology used, we employed different approaches implemented in the 
    Squidpy library (v.1.3.0) to construct the network. For the Visium and RNA + ATAC data, 
    which exhibit a regular structure, we assigned the six and four closest surrounding spots, 
    respectively, as neighbors for each spot. For the CODEX, CosMx, MERFISH and IMC data, we 
    constructed the network using Delaunay triangulation. However, Delaunay triangulation 
    can result in long edges between cells, especially at the slide borders. Therefore, we 
    eliminated edges between nodes if the distance between them exceeded the 99th percentile 
    of all distances between connected nodes.
    """
    if technology in ["CODEX", "CosMx", "MERFISH", "IMC"]:
        adata_domains = adata[adata.obs["spatial_cluster"].notna()].copy()
        
        # Initial version from cellcharter
        #sq.gr.spatial_neighbors(adata, delaunay=True, coord_type='generic')
        #cc.gr.remove_long_links(adata, distance_percentile = 99.0)
        
        # Alternative cellcharter like version that should be more robust than percentiles
        #sq.gr.spatial_neighbors(adata_domains, delaunay=True, coord_type='generic')
        #remove_long_links(adata_domains, distance_percentile=None, mean_distance_factor=2.0)
        
        # Simple n neighbors version (danger: long range connections are not removed)
        sq.gr.spatial_neighbors(adata_domains, delaunay=False, n_neighs=20, coord_type='generic')

    elif technology in ["Visium"]:
        sq.gr.spatial_neighbors(adata, n_neighs=6, coord_type="grid")
    elif technology in ["RNA + ATAC"]: 
        sq.gr.spatial_neighbors(adata, n_neighs=4, coord_type="grid")
    else:
        raise NotImplementedError(f"Graph connectivity calculation for the given technology {technology} is not implemented.")
    
    # Compute connected components
    if technology in ["CODEX", "CosMx", "MERFISH", "IMC"]:
        min_cells = 100
        cc.gr.connected_components(adata_domains, cluster_key='spatial_cluster',min_cells = min_cells)
        adata.obs["component"] = adata_domains.obs["component"]
    elif technology in ["Visium", "RNA + ATAC"]:
        min_cells = 3
        cc.gr.connected_components(adata, cluster_key='spatial_cluster',min_cells = min_cells)
    else:
        raise NotImplementedError(f"Min Nr of cells for connected component identification for the given technology {technology} is not implemented.")
        
    # Compute boundaries of domain clusters (saved in adata.uns[f"shape_component"]["boundary"]}
    a_start = get_alpha_start(adata)    
    cc.tl.boundaries(adata, cluster_key = "component", min_hole_area_ratio = 0.1, alpha_start = a_start, copy = False)
    
    # Compute purities for each connected component
    purities = cc.tl.purity(adata, cluster_key="component", library_key="sample", exterior=False, copy=True)
    
    # Compute purity metric (mean over purities of all connected components)
    metric: float = np.mean([p for _,p in purities.items()]) 
        
    # Write files for debug purposes
    if fig_path:
        plot_components(adata, save=fig_path)
    # If adata is saved for debug purposes writing fails when the polygons are still in adata
    # del adata.uns[f"shape_component"]
    # adata.write(args.out_file.replace(".tsv",".h5ad")) #TODO: DELETE
        
    ## Write output
    from pathlib import Path
    
    Path(args.out_file).parent.mkdir(parents=True, exist_ok=True)
    
    with open(args.out_file, "w") as file:
        file.write(f"{metric:.5e}\n")
    
# NOTE: __name__ == "__main__" is required, otherwise the multiprocessing might fail
if __name__ == "__main__":
    script()