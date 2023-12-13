import os
import argparse
import shutil
import scipy
import anndata 
import numpy as np
import pandas as pd
import json


BAD_GENES = ["eGFP","mCherry2","tdToma"]
META_DICT = {"technology":"Merfish"}
SAMPLE_INFO = {"patient":"1","sample":"1","position":"0","replicate":"1","n_clusters":"0","directory":f"os.path.basename(args.output)"}
LICENSE = """
This dataset was created by AG Kuppe at the University Hospital Aachen, Germany.

It may only be used in the context of the Spacehack 2023 event.

In case of any questions feel free to contact Paul Kiessling, pakiessling@ukaachen.de.
"""


def copy_images(input_folder, output_folder):
    # Ensure the output folder exists, create if not
    os.makedirs(output_folder, exist_ok=True)

    # Get a list of files in the input folder
    files = os.listdir(input_folder)
    files = [file for file in files if file.endswith(".tif")]
    # Copy image
    for file in files:
        input_path = os.path.join(input_folder, file)
        output_path = os.path.join(output_folder, file)
        shutil.copy2(input_path, output_path)
        print(f"Copied: {input_path} to {output_path}")


def load_into_anndata(input_folder):
    data = pd.read_csv(input_folder + "/cell_by_gene.csv", index_col=0, dtype={"cell": str})
    obs = pd.read_csv(input_folder + "/cell_metadata.csv", index_col=0, dtype={"EntityID": str})
    is_gene = ~data.columns.str.lower().str.contains("blank")
    adata = anndata.AnnData(data.loc[:, is_gene], dtype=data.values.dtype, obs=obs)
    adata.obsm["blank"] = data.loc[:, ~is_gene]
    adata = adata[:,~adata.var_names.isin(BAD_GENES)]
    adata.obsm["spatial"] = adata.obs[["center_x", "center_y"]].values
    adata.obs["EntityID"] = adata.obs.index
    return adata

def convert_data(input_folder, output_folder,ct_file):
    os.makedirs(output_folder, exist_ok=True)
    adata = load_into_anndata(input_folder)

    # Observations
    obs = adata.obs.copy()
    obs["selected"] = "true"
    if ct_file != None:
        print("adding ct")
        ct = pd.read_table(ct_file, index_col=0)
        ct.index = ct.index.astype("str")
        obs["cell_type"] = ct["cell_type"]
        obs['cell_type'].fillna('filtered', inplace=True)

        adata.obs["cell_type"] = ct["cell_type"]
        adata.obs["cell_type"].fillna('filtered', inplace=True)
        
    obs.to_csv(f"{output_folder}/observations.tsv",sep="\t",index_label="")

    # Features
    vars = adata.var.copy()
    vars["selected"] = "true"

    vars.to_csv(f"{output_folder}/features.tsv",sep="\t",index_label="")

    # Coordinates
    coords = pd.DataFrame(adata.obsm["spatial"],columns=["x","y"])
    coords.index = adata.obs.index
    coords.to_csv(f"{output_folder}/coordinates.tsv",sep="\t",index_label="")

    # Matrix
    scipy.io.mmwrite(f"{output_folder}/counts.mtx",adata.X)

    # Anndata
    adata.write_h5ad(f"{output_folder}/anndata.h5ad")


def write_json(dict,output_path):
    with open(output_path, 'w') as json_file:
        json.dump(dict, json_file)


 

def main():
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(description="Convert Vizgen Merfish Data to Spacehack format.")
    
    # Add arguments for input and output folders
    parser.add_argument("--input", help="Path to the input folder containing Vizgen Merscope output",required=True)
    parser.add_argument("--output", help="Path to the output folder",required=True)
    parser.add_argument("--ct", help="Path to tsv containing cell-barcode and ct,columname should be 'cell'",required=False)
    

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the function to copy files
    convert_data(args.input, args.output,args.ct)
    #copy_images(os.path.join(args.input_folder, "images"), args.output_folder)

    # write json
    write_json(META_DICT,os.path.join(os.path.dirname(args.output), "experiment.json"))

    # write samples.tsv
    sample_df = pd.DataFrame.from_dict(SAMPLE_INFO, orient='index').T
    output_directory = os.path.dirname(args.output)
    sample_df.to_csv(f"{output_directory}/samples.tsv", sep="\t", index_label=False)

    # write LICENSE
    with open(f"{os.path.dirname(args.output)}/LICENSE.md", 'w') as file:
        file.write(LICENSE)

if __name__ == "__main__":
    main()
