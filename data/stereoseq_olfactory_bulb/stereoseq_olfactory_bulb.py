#!/usr/bin/env python

import urllib.request
from urllib.parse import urlparse
import os
import anndata
import argparse
import shutil
import pandas as pd
import scipy
import json
import tempfile

# 6 available but only 2 contain region label and coordinates
LINKS = [
"https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/Mouse_olfa_S1.h5ad",
"https://ftp.cngb.org/pub/SciRAID/stomics/STDS0000058/stomics/Mouse_olfa_S2.h5ad"]


META_DICT = {"technology":"Stereo-seq"}

SAMPLE_COLUMNS = ["patient","sample","position","replicate","n_clusters","directory"]


def download_links(links, temp_dir):
    for link in links:
        try:
            response = urllib.request.urlopen(link)
            # Extract filename from the URL
            filename = os.path.join(temp_dir, urlparse(link).path.split("/")[-1])
            with open(filename, 'wb') as file:
                file.write(response.read())
            print(f"Downloaded: {filename}")
        except Exception as e:
            print(f"Error downloading {link}: {e}")

def process_adata(adata_path,output_folder,iteration,sample_df):
    folder_name = os.path.splitext(os.path.basename(adata_path))[0]
    complete_path = os.path.join(output_folder,folder_name)
    os.makedirs(complete_path, exist_ok=True)
    adata = anndata.read_h5ad(adata_path)

    # Observations
    obs = adata.obs.copy()
    obs["selected"] = "true"
    obs.to_csv(f"{complete_path}/observations.tsv",sep="\t",index_label="")

    # Features
    vars = adata.var.copy()
    vars["selected"] = "true"
    vars.to_csv(f"{complete_path}/features.tsv",sep="\t",index_label="")

    # Coordinates
    coords = pd.DataFrame(adata.obsm["spatial"],columns=["x","y"])
    coords.index = adata.obs.index
    coords.to_csv(f"{complete_path}/coordinates.tsv",sep="\t",index_label="")

    # Matrix

    # Check if "count" key exists in adata.layers
    if "count" in adata.layers:
        matrix_to_write = adata.layers["count"]
        file_path = f"{complete_path}/counts.mtx"
        scipy.io.mmwrite(file_path, matrix_to_write)
    elif "counts" in adata.layers:
        matrix_to_write = adata.layers["counts"]
        file_path = f"{complete_path}/counts.mtx"
        scipy.io.mmwrite(file_path, matrix_to_write)
        print(f"Matrix written to {file_path}")
    else:
        print("Neither 'count' nor 'counts' key found in adata.layers.")


    # add info for sample.tsv
    # Your sample_data_basis dictionary
    sample_data_basis = {"patient": iteration, "sample": "1", "position": "0", "replicate": "2", "n_clusters": adata.obs.annotation.nunique(), "directory": folder_name}
    
    # Creating a DataFrame from the dictionary
    sample_data = pd.DataFrame([sample_data_basis])
    
    # Concatenating the new DataFrame to sample_df
    sample_df.iloc[iteration] = sample_data_basis

    # Write labels.tsv
    if "annotation" in adata.obs.columns:
        labels = adata.obs["annotation"]
        labels.to_csv(f"{complete_path}/labels.tsv",sep="\t",index_label="")
        
    

def write_json(dict,output_path):
    with open(output_path, 'w') as json_file:
        json.dump(dict, json_file)
    

def main():
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(description="Convert Stereo-seq Olfactory bulb data to Spacehack format.")
    
    # Add arguments for output folder
    parser.add_argument('-o','--out_dir', help="Output directory to write files to.",required=True)
    

    # Parse the command-line arguments
    args = parser.parse_args()

    # Download and process
    with tempfile.TemporaryDirectory() as temp_dir: 
        download_links(LINKS,temp_dir)
        os.makedirs(args.out_dir, exist_ok=True)
    
        
        sample_df = pd.DataFrame(columns=SAMPLE_COLUMNS,index=range(len(LINKS)))
        anndatas = [os.path.join(temp_dir, file) for file in os.listdir(temp_dir) if file.endswith(".h5ad")]
        for iteration, adata in enumerate(anndatas):
            process_adata(adata, args.out_dir,iteration,sample_df)
        
    
        # write json 
        write_json(META_DICT,f"{args.out_dir}/experiment.json")
    
        # write samples.tsv
        sample_df.to_csv(f"{args.out_dir}/samples.tsv", sep="\t", index_label="")
    

    


if __name__ == "__main__":
    main()



    