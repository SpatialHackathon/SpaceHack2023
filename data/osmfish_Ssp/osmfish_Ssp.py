#!/usr/bin/env python

# Made by Paul Kiessling pakiessling@ukaachen.de 


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


LINKS = [
"https://linnarssonlab.org/osmFISH/osmFISH_SScortex_mouse_all_cells.loom"]


META_DICT = {"technology":"osmFISH"}
LICENSE = """

"""
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
    adata = anndata.read_loom(adata_path)
    # Observations
    obs = adata.obs.copy()
    obs["selected"] = "true"
    obs.loc[obs["Valid"] == 0, "selected"] = "false"
    obs.to_csv(f"{complete_path}/observations.tsv",sep="\t",index_label="")

    # Features
    vars = adata.var.copy()
    vars["selected"] = "true"
    vars.to_csv(f"{complete_path}/features.tsv",sep="\t",index_label="")

    # Coordinates
    coords = pd.DataFrame(adata.obs,columns=["X","Y"])
    coords.index = adata.obs.index
    coords.to_csv(f"{complete_path}/coordinates.tsv",sep="\t",index_label="")

    # Matrix
    scipy.io.mmwrite(f"{complete_path}/counts.mtx",adata.X)

    # add info for sample.tsv
    # Your sample_data_basis dictionary
    sample_data_basis = {"patient": iteration, "sample": "1", "position": "0", "replicate": "1", "n_clusters": adata.obs.Region.nunique(), "directory": folder_name}
    
    # Creating a DataFrame from the dictionary
    sample_data = pd.DataFrame([sample_data_basis])
    
    # Concatenating the new DataFrame to sample_df
    sample_df.iloc[iteration] = sample_data_basis

    # Write labels.tsv
    if "Region" in adata.obs.columns:
        labels = adata.obs["Region"]
        labels.to_csv(f"{complete_path}/labels.tsv",sep="\t",index_label="")
        
    

def write_json(dict,output_path):
    with open(output_path, 'w') as json_file:
        json.dump(dict, json_file)
    

def main():
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(description="Convert osmFISH data to Spacehack format.")
    
    # Add arguments for output folder
    parser.add_argument('-o','--out_dir', help="Output directory to write files to.",required=True)
    

    # Parse the command-line arguments
    args = parser.parse_args()

    # Download and process
    with tempfile.TemporaryDirectory() as temp_dir: 
        download_links(LINKS,temp_dir)
        os.makedirs(args.out_dir, exist_ok=True)
    
        
        sample_df = pd.DataFrame(columns=SAMPLE_COLUMNS,index=range(len(LINKS)))
        anndatas = [os.path.join(temp_dir, file) for file in os.listdir(temp_dir) if file.endswith(".loom")]
        for iteration, adata in enumerate(anndatas):
            process_adata(adata, args.out_dir,iteration,sample_df)
        
    
        # write json 
        write_json(META_DICT,f"{args.out_dir}/experiment.json")
    
        # write samples.tsv
        sample_df.to_csv(f"{args.out_dir}/samples.tsv", sep="\t", index_label=False)
    
        # write LICENSE
        with open(f"{args.out_dir}/LICENSE.md", 'w') as file:
            file.write(LICENSE)
    


if __name__ == "__main__":
    main()



    