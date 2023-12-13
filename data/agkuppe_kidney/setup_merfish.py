import os
import argparse
import shutil
import scipy
import squidpy as sq # Needs to be installed from github not pypi
import numpy as np
import pandas as pd

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


def convert_data(input_folder, output_folder,ct_file):
    os.makedirs(output_folder, exist_ok=True)
    adata = sq.read.vizgen(input_folder,counts_file="cell_by_gene.csv",meta_file="cell_metadata.csv")

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
    bad_genes = ["eGFP","mCherry2","tdToma"]
    vars.loc[vars.index.isin(bad_genes), "selected"] = "false"
    vars.to_csv(f"{output_folder}/features.tsv",sep="\t",index_label="")

    # Coordinates
    coords = pd.DataFrame(adata.obsm["spatial"],columns=["x","y"])
    coords.index = adata.obs.index
    coords.to_csv(f"{output_folder}/coordinates.tsv",sep="\t",index_label="")

    # Matrix
    scipy.io.mmwrite(f"{output_folder}/counts.mtx",adata.X)

    # Anndata
    adata.write_h5ad(f"{output_folder}/anndata.h5ad")
 

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
    

if __name__ == "__main__":
    main()
