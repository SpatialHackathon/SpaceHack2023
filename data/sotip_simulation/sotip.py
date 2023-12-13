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


# 6 available but only 2 contain region label and coordinates
LINKS = [
"https://github.com/TencentAILabHealthcare/SOTIP/raw/master/SOTIP_analysis/simulation/simulation_data/simulation2.h5ad",
"https://github.com/TencentAILabHealthcare/SOTIP/raw/master/SOTIP_analysis/simulation/simulation_data/simulation5.h5ad"]

WORKDIR = "./.workdir"
META_DICT = {"technology":"Visium"}
LICENSE = """
Copyright 2022 Zhiyuan Yuan

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""
SAMPLE_COLUMNS = ["patient","sample","position","replicate","n_clusters","directory"]


def download_links(links):
    # Create a directory if it doesn't exist
    os.makedirs(WORKDIR, exist_ok=True)

    for link in links:
        try:
            response = urllib.request.urlopen(link)
            # Extract filename from the URL
            filename = os.path.join(WORKDIR, urlparse(link).path.split("/")[-1])
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
    scipy.io.mmwrite(f"{complete_path}/counts.mtx",adata.X)

    # Anndata
    adata.write_h5ad(f"{complete_path}/anndata.h5ad")

    # add info for sample.tsv
    # Your sample_data_basis dictionary
    sample_data_basis = {"patient": iteration, "sample": "1", "position": "0", "replicate": "1", "n_clusters": adata.obs.region_id.nunique(), "directory": folder_name}
    
    # Creating a DataFrame from the dictionary
    sample_data = pd.DataFrame([sample_data_basis])
    
    # Concatenating the new DataFrame to sample_df
    sample_df.iloc[iteration] = sample_data_basis
    

def write_json(dict,output_path):
    with open(output_path, 'w') as json_file:
        json.dump(dict, json_file)
    

def main():
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(description="Convert Vizgen Merfish Data to Spacehack format.")
    
    # Add arguments for output folder
    parser.add_argument("--output", help="Path to the output folder",required=True)
    

    # Parse the command-line arguments
    args = parser.parse_args()

    # Download and process
    download_links(LINKS)
    os.makedirs(args.output, exist_ok=True)

    
    sample_df = pd.DataFrame(columns=SAMPLE_COLUMNS,index=range(len(LINKS)))
    anndatas = [os.path.join(WORKDIR, file) for file in os.listdir(WORKDIR) if file.endswith(".h5ad")]
    for iteration, adata in enumerate(anndatas):
        process_adata(adata, args.output,iteration,sample_df)
    

    # write json 
    write_json(META_DICT,f"{args.output}/experiment.json")

    # write samples.tsv
    sample_df.to_csv(f"{args.output}/samples.tsv", sep="\t", index_label=False)

    # write LICENSE
    with open(f"{args.output}/LICENSE.md", 'w') as file:
        file.write(LICENSE)

    # Delete working directory
    shutil.rmtree(WORKDIR)

if __name__ == "__main__":
    main()



    