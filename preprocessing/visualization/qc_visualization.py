#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Jieran Sun; Implemented visualization

import argparse

# TODO adjust description
parser = argparse.ArgumentParser(description="Visualization of preprocessing-QC")

parser.add_argument(
    "-c", "--coordinates", help="Path to coordinates (as tsv).", required=True
)
parser.add_argument(
    "-m", "--matrix", help="Path to counts (as mtx).", required=True
)
parser.add_argument(
    "-f", "--features", help="Path to features (as tsv).", required=True
)
parser.add_argument(
    "-o", "--observations", help="Path to observations (as tsv).", required=True
)
parser.add_argument(
    "-d", "--out_dir", help="Output directory.", required=True
)
parser.add_argument(
    "--opt", help="optarg.json file for the dataset quality control", required=False
)
parser.add_argument(
    "--feature_qc", help="Path to features (as tsv) after quality control",
    required=True,
)
parser.add_argument(
    "--observation_qc", help="Path to observation (as tsv) after quality control",
    required=True,
)
parser.add_argument(
    "--technology", help="which technology is used for this dataset",
    required=True,
)

args = parser.parse_args()

from pathlib import Path

out_dir = Path(args.out_dir)

# Output files
vis_pdf = out_dir / "pp_report_sample.pdf"

# Use these filepaths as input ...
coord_file = args.coordinates
matrix_file = args.matrix
feature_file = args.features
observation_file = args.observations
technology = args.technology
spotsize_dic = {
    "Visium": 135,
    "MERFISH": 50,
    "Xenium": 25,
}

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
adata.var_names_make_unique()

## Your code goes here
import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

############### sample preparation ###############
obs_qc = pd.read_table(args.observation_qc, index_col=0).index.tolist()
fea_qc = pd.read_table(args.feature_qc, index_col=0).index.tolist()

adata.var['QC'] = adata.var.index.isin(fea_qc)
adata.obs['QC'] = adata.obs.index.isin(obs_qc)
adata.var['QC'] = adata.var['QC'].astype('category')
adata.obs['QC'] = adata.obs['QC'].astype('category')

# Find out mitochondria percentage
mito_detect_list = [adata.var[col].astype(str).str.startswith("MT-") 
                    for col in adata.var.columns
                   ] + [adata.var_names.astype(str).str.startswith("MT-")]
mito_sum = [np.sum(mito_detect) for mito_detect in mito_detect_list]

qc_var = ()
if not np.all(mito_sum == 0):
    adata.var["mt"] = mito_detect_list[np.argmax(mito_sum)]
    qc_var = ["mt"]

# Calculate QC metric for plotting
sc.pp.calculate_qc_metrics(
adata, qc_vars=qc_var, inplace=True, percent_top=None, log1p=True
)
#adata.write_h5ad("adata.h5ad")

# Find the QC threshold:
min_counts = adata.obs[adata.obs["QC"]]["total_counts"].min()
min_genes  = adata.obs[adata.obs["QC"]]["n_genes_by_counts"].min()
min_cells  = adata.var[adata.var["QC"]]["n_cells_by_counts"].min()

min_counts_log1p = adata.obs[adata.obs["QC"]]["log1p_total_counts"].min()
min_genes_log1p  = adata.obs[adata.obs["QC"]]["log1p_n_genes_by_counts"].min()

############### Plotting visualization ###############
# TODO: adding line for current QC
fig, axs = plt.subplots(2, 4, figsize=(25, 13))
# plt.rcParams["figure.figsize"] = (8, 8)
spot_size = spotsize_dic.pop(technology, 80)
axf = axs.flatten()

# Total expression counts per cell
sns.histplot(adata.obs["total_counts"], kde=False, ax=axf[0])
axf[0].axvline(x=min_counts, color='r', linestyle='--')
axf[0].set_title("Cell: Total counts")

# no. of genes expressed in the cell, for all cells
sns.histplot(adata.obs["n_genes_by_counts"], kde=False, bins=60, ax=axf[1])
axf[1].axvline(x=min_genes, color='r', linestyle='--')
axf[1].set_title("Cell: # gene expressed")

# Total expression counts per gene
sns.scatterplot(data=adata.obs, x="log1p_total_counts", y = "n_genes_by_counts", color="k", alpha = 0.8, ax= axf[2], linewidth=0)
sns.kdeplot(data=adata.obs, x="log1p_total_counts", y = "n_genes_by_counts", ax= axf[2])
axf[2].axvline(x=min_counts_log1p, color='r', linestyle='--')
axf[2].axhline(y=min_genes_log1p, color='r', linestyle='--')
axf[2].set_title("Cell: total vs #gene per cell")

sns.scatterplot(data=adata.var, x="log1p_total_counts", y = "n_cells_by_counts", color="k", alpha = 0.8, ax= axf[3], linewidth=0)
sns.kdeplot(data=adata.var, x="log1p_total_counts", y = "n_cells_by_counts", ax= axf[3])
# sns.kdeplot(data=adata.var, x="total_counts", y = "n_cells_by_counts", ax= axf[3])
axf[3].axhline(y=min_cells, color='r', linestyle='--')
axf[3].set_title("Gene: total vs #cells per gene")

sc.pl.spatial(adata, color="total_counts", spot_size = spot_size, show = False,
               title = "Cell: Total expression", ax= axf[4])

sc.pl.spatial(adata, color="n_genes_by_counts", spot_size = spot_size, show = False,
               title = "Cell: #gene expressed", ax= axf[5])

sc.pl.spatial(adata, color="QC", spot_size = spot_size, show = False,
               title = "QC results", ax= axf[6])

# For the last plot
if "mt" in adata.var.keys() and sum(adata.var["mt"]>0)>0:
    sns.scatterplot(x=adata.obs["log1p_total_counts"], 
                    y =np.sqrt(adata.obs["pct_counts_mt"]), 
                    ax= axf[7], color  = "black", alpha = 0.7, linewidth=0)
    sns.kdeplot(x=adata.obs["log1p_total_counts"], 
                    y =np.sqrt(adata.obs["pct_counts_mt"]), 
                    ax= axf[7])
    max_mt = np.sqrt(adata.obs[adata.obs["QC"]]["pct_counts_mt"].max())
    axf[7].axhline(y=max_mt, color='r', linestyle='--')
    axf[7].set_title("mt: percentage")
    # sc.pl.violin(adata, "pct_counts_mt",show = False, ax= axf[7], jitter=False)

elif "cell_count" in adata.obs.keys():
    sc.pl.violin(adata, "cell_count",show = False, title = "gene count distribution", ax= axf[7])

elif "in_tissue" in adata.obs.keys():
    adata.obs["in_tissue"] = adata.obs["in_tissue"].astype('category')
    sc.pl.spatial(adata, color="in_tissue", spot_size = spot_size, show = False,
               title = "in_tissue results", ax= axf[7])

elif "cell_area" in adata.obs.keys():
    sns.scatterplot(x=adata.obs["log1p_total_counts"],
                    y =adata.obs["cell_area"],
                    ax= axf[7], color  = "black", alpha = 0.7, linewidth=0)
    sns.kdeplot(x=adata.obs["log1p_total_counts"],
                    y =adata.obs["cell_area"],
                    ax= axf[7])
    max_ca = adata.obs[adata.obs["QC"]]["cell_area"].max()
    min_ca = adata.obs[adata.obs["QC"]]["cell_area"].min()
    axf[7].axhline(y=max_ca, color='r', linestyle='--')
    axf[7].axhline(y=min_ca, color='r', linestyle='--')
    axf[7].set_title("cell area")

elif "control_probe_counts" in adata.obs.keys():
    sns.histplot(adata.obs["control_probe_counts"], kde=False, bins=60, ax=axf[7])
    axf[7].axvline(x=adata.obs[adata.obs["QC"]]["control_probe_counts"].max(), color='r', linestyle='--')
    axf[7].set_title("Control probe numbers")

else:
    # no. of positively expressed cells of this gene, for all genes
    sns.histplot(adata.var["n_cells_by_counts"], kde=False, bins=60, ax=axf[7])
    axf[7].axvline(x=min_cells, color='r', linestyle='--')  # Add vertical line
    axf[7].set_title("Gene: # cell expressed")

fig.suptitle(f"{out_dir.parent.parent.name}: {out_dir.parent.name}")
plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Add some space at the top for the title

# Save the figure to a PDF file
out_dir.mkdir(parents=True, exist_ok=True)
pdf_file = out_dir / "pp_report_sample.pdf"
with PdfPages(pdf_file) as pdf:
    pdf.savefig(fig)