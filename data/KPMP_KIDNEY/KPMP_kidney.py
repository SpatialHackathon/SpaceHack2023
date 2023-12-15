#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Fadhl Alakwaa

import argparse

# TODO adjust description
parser = argparse.ArgumentParser(description="Load visium data for kidney tissue from the Kidney Precsion Medicine Project (KPMP) ")

parser.add_argument(
    "-o", "--out_dir", help="'/home/jovyan/scratch/SpaceHack2/data/KPMP_KIDNEY'", required=True
)

args = parser.parse_args()


from pathlib import Path

import pandas as pd

out_dir = Path(args.out_dir)

# The folder structure should look like the following
# out_dir
# |_______sample_1  (sample name can be chosen freely)
#         |_____coordinates.tsv
#         |_____features.tsv
#         |_____observations.tsv
#         |_____counts.mtx  (use scipy.io.mmwrite)
#         |_____labels.tsv  (optional)
#         |_____H_E.(tiff/png/...)  (optional)
#         |_____H_E.json  (optional, required if H_E is provided)
# |_______sample_2
#         | ...
# |_______samples.tsv
# |_______experiment.json
# if additional output files are required write it also to out_dir


## Your code goes here
import pandas as pd
import os
import anndata
from scipy.io import mmwrite
import shutil
import requests
from io import BytesIO
import json

technology='Visium'
out_dir = '/home/jovyan/scratch/SpaceHack2/data/KPMP_KIDNEY'

# Read the CSV file
csv_path = "KPMP_hachathon_tiff_h5ad.csv"
df = pd.read_csv(csv_path)

# Iterate over rows in the DataFrame
for i in range(0, len(df), 2):
    os.chdir(out_dir)
    
    # Create a folder with the "Participant ID"
    folder_name = str(df.loc[i, 'Participant_ID'])
    os.makedirs(folder_name, exist_ok=True)

    # Construct the paths for H5AD and TIFF files
    h5_file_path = os.path.join(folder_name, f"{df.loc[i, 'File_Name']}")
    tif_file_path = os.path.join(folder_name, f"{df.loc[i + 1, 'File_Name']}")
    
    # Read the H5AD file from the URL into an AnnData object
    h5_url = f"https://atlas.kpmp.org/api/v1/file/download/{df.loc[i, 'Package_ID']}/{df.loc[i, 'File_Name']}"
    response = requests.get(h5_url)
    adata = anndata.read_h5ad(BytesIO(response.content))
    # Read the H5AD file into an AnnData object
    #adata = anndata.read_h5ad(h5_file_path)
    # Download the TIFF file using wget
    tif_url = f"https://atlas.kpmp.org/api/v1/file/download/{df.loc[i + 1, 'Package_ID']}/{df.loc[i + 1, 'File_Name']}"
    urlretrieve(tif_url, tif_file_path)

    # Save coordinates.tsv
    adata.obsm.to_df().rename(columns={'spatial1': 'x', 'spatial2': 'y'}).to_csv(os.path.join(folder_name, 'coordinates.tsv'), sep='\t', index=True, header=True)

    # Save counts.mtx in Matrix Market format
    adata.X = adata.X.astype(int)  # Ensure counts are integers
    mmwrite(os.path.join(folder_name, 'counts.mtx'), adata.X)

    # Save features.tsv
    adata.var.to_csv(os.path.join(folder_name, 'features.tsv'), sep='\t', index=True, header=True)

    # Save observations.tsv
    adata.obs = adata.obs.iloc[:, 1:]
    adata.obs.rename(columns={'array_row': 'row', 'array_col': 'col'}).to_csv(os.path.join(folder_name, 'observations.tsv'), sep='\t', index=True, header=True)

    # Copy TIFF file to the folder

donors=[]
samples =[]
positions = []
replicates = []
directories = []
n_clusters=[]
csv_path = "KPMP_hachathon_tiff_h5ad.csv"
df = pd.read_csv(csv_path)

# Iterate over rows in the DataFrame
for i in range(0, len(df), 2):
    os.chdir(parent_directory)
    donor=df.loc[i, 'Participant_ID']
    sample=df.loc[i, 'Participant_ID']
    donors.append(donor)
    samples.append(sample)
    #positions.append(np.nan)
    #replicates.append(replicate_counter)
    #replicate_counter+=1
    #directories.append(donor_section)
    n_clusters.append(5)
samples_df = pd.DataFrame(data={'patient': donors, 
                                    'sample':samples, 
                                    'n_clusters':n_clusters
                                   })    
                                   
samples_df.to_csv(os.path.join(out_dir, "samples.tsv"), sep='\t', index_label='')
with open(os.path.join(out_dir, "experiment.json"), 'w') as f:
      exp_info = {'technology': technology,
                    'species': 'Homo sapiens'}
      json.dump(exp_info, f)