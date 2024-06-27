#!/usr/bin/env python

# Made by Paul Kiessling pakiessling@ukaachen.de 


import urllib.request
from urllib.parse import urlparse
import os
import scanpy as sc
import argparse
import shutil
import pandas as pd
import scipy
import json
import tempfile


LINKS = [
"https://raw.githubusercontent.com/JinmiaoChenLab/SEDR_analyses/master/data/BRCA1/metadata.tsv",
"https://github.com/JinmiaoChenLab/SEDR_analyses/raw/master/data/BRCA1/V1_Human_Breast_Cancer_Block_A_Section_1/filtered_feature_bc_matrix.h5",
"https://raw.githubusercontent.com/JinmiaoChenLab/SEDR_analyses/master/data/BRCA1/V1_Human_Breast_Cancer_Block_A_Section_1/spatial/scalefactors_json.json",
    "https://raw.githubusercontent.com/JinmiaoChenLab/SEDR_analyses/master/data/BRCA1/V1_Human_Breast_Cancer_Block_A_Section_1/spatial/tissue_hires_image.png",
    "https://raw.githubusercontent.com/JinmiaoChenLab/SEDR_analyses/master/data/BRCA1/V1_Human_Breast_Cancer_Block_A_Section_1/spatial/tissue_lowres_image.png",
    "https://raw.githubusercontent.com/JinmiaoChenLab/SEDR_analyses/master/data/BRCA1/V1_Human_Breast_Cancer_Block_A_Section_1/spatial/tissue_positions_list.csv"
]


META_DICT = {"technology":"Visium"}
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

def process_adata(input_path,output_folder,sample_df):
    complete_path = os.path.join(output_folder,"visium_breast_cancer_SEDR")
    # Bring data in correct structure for scanpy
    matrix_path = os.path.join(  input_path,"spatial/filtered_feature_bc_matrix.h5")
    parent_directory = os.path.dirname(matrix_path)
    shutil.move(matrix_path, input_path)


    
    os.makedirs(complete_path, exist_ok=True)
    adata = sc.read_visium(input_path, count_file="filtered_feature_bc_matrix.h5",source_image_path=os.path.join(complete_path,"H_E.tiff"))
    adata.var_names_make_unique()
    obs = adata.obs.copy()
    obs["selected"] = "true"
    obs.rename(columns={'array_row': 'row', 'array_col': 'col'}, inplace=True)
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

    # add info for sample.tsv
    domain_annotation = pd.read_table(os.path.join(input_path,"spatial/metadata.tsv"))

    # Write labels.tsv
    labels = domain_annotation
    labels.index = adata.obs.index
    labels = labels.rename(columns={'fine_annot_type':'label'})
    labels.to_csv(f"{complete_path}/labels.tsv",sep="\t",index_label="")

    sample_data_basis = {"patient": "1", "sample": "1", "position": "0", "replicate": "1", "n_clusters": labels.label.nunique(), "directory": "visium_breast_cancer_SEDR"}
    
    # Concatenating the new DataFrame to sample_df
    sample_df.iloc[0] = sample_data_basis
        
    # Move image files
    shutil.move(os.path.join(input_path,"spatial/tissue_hires_image.png"),os.path.join(complete_path,"H_E.tiff"))
    shutil.move(os.path.join(input_path,"spatial/scalefactors_json.json"),os.path.join(complete_path,"H_E.json"))
    
    
def write_json(dict,output_path):
    with open(output_path, 'w') as json_file:
        json.dump(dict, json_file)
    

def main():
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(description="Convert Visium data to Spacehack format.")
    
    # Add arguments for output folder
    parser.add_argument('-o','--out_dir', help="Output directory to write files to.",required=True)
    

    # Parse the command-line arguments
    args = parser.parse_args()

    # Download and process
    with tempfile.TemporaryDirectory() as temp_dir:
        os.makedirs(os.path.join(temp_dir,"spatial"), exist_ok=True)
        download_links(LINKS,os.path.join(temp_dir,"spatial"))
        os.makedirs(args.out_dir, exist_ok=True)

        
        
        sample_df = pd.DataFrame(columns=SAMPLE_COLUMNS,index=[0])

        process_adata(temp_dir, args.out_dir,sample_df)
        
    
        # write json 
        write_json(META_DICT,f"{args.out_dir}/experiment.json")
    
        # write samples.tsv
        sample_df.to_csv(f"{args.out_dir}/samples.tsv", sep="\t", index_label=False)
    
        # write LICENSE
        with open(f"{args.out_dir}/LICENSE.md", 'w') as file:
            file.write(LICENSE)
    


if __name__ == "__main__":
    main()



    