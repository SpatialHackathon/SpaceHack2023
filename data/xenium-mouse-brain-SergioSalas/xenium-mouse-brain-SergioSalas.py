#!/usr/bin/env python

import argparse
import os
import sys
from pathlib import Path
import anndata
import json
import gdown
import pandas as pd
import scipy
import shutil
import tempfile

def convert_data(out_dir):
    
    # TODO: Add downloader as soon as data is published on zenodo. 
    # THIS IS A PRELIMINARY LINK TO SEBASTIAN's PERSONAL GDRIVE....
    file_download_link = "https://drive.google.com/drive/folders/13dT5PgHypGpp_tsHxoYho8z3JwnvX34r?usp=sharing"
    # destination_path = "./"

    # temp_dir = Path(out_dir)/"temp"# Path(__file__).parent.resolve()    
    # os.makedirs(out_dir,)

    with tempfile.TemporaryDirectory() as temp_dir:

        # download the folder 'temporary_file_host_spacehack'
        print(f'Downloading annotated data from gdrive into {temp_dir}...')
        gdown.download_folder(url=file_download_link, output=temp_dir, quiet=False)
        print('...done')

        temp_dir=Path(temp_dir)
    
        in_file = temp_dir / 'adata_multisection_nuclei_r1_with_annotations.h5ad'
        adata = anndata.read_h5ad(in_file)
        adata.obs_names = "cell_" + adata.obs_names.astype(str)

        out_dir = Path(out_dir)
        project_root = out_dir 

        os.makedirs(project_root, exist_ok=True)
        
        # create experiment.json
        experiment_json_dict = {"technology":"Xenium"}
        with open(project_root / 'experiment.json', 'w') as f:
            json.dump(experiment_json_dict, f)

        # create samples.tsv
        df_sample_tsv = pd.DataFrame({"sample": "region_main",
                                     "directory": "region_main",
                                     "n_clusters": adata.obs.region_level2.unique().shape[0]},
                                    index = ["region_main"])
        df_sample_tsv.to_csv(project_root/'samples.tsv',sep="\t",index_label="")
    
        
        for i,sample in enumerate(df_sample_tsv.index):
    
            print(f'Processing sample {i}: "{sample}"')
    
            sample_output_folder = project_root / sample
    
            os.makedirs(sample_output_folder,exist_ok=True)
    
            # Labels:
            df_labels_tsv = adata.obs[['region_level2']]
            df_labels_tsv.columns=[['label']]
            df_labels_tsv['label_confidence']=(~adata.obs.region_level2.isna()).astype(int)
            print(f"Writing labels: {df_labels_tsv.head()}")
            df_labels_tsv.to_csv(sample_output_folder/'labels.tsv',sep="\t",index_label="")
    
            # Observations
            obs = adata.obs.copy()
            obs["selected"] = "true"
            df_observations_tsv = adata.obs #[['Class','cell_area','n_counts']]
            # df_observations_tsv.columns = [['cell_type','cell_area','transcript_counts']]
            print(f"Writing observations: {df_observations_tsv.head()}")           
            df_observations_tsv.to_csv(sample_output_folder/'observations.tsv',sep="\t",index_label="")
        
            # Features
            vars = adata.var.copy()
            vars["selected"] = "true"
            df_features_tsv = adata.var[[]]
            print(f"Writing features: {df_features_tsv.head()}")
            df_features_tsv.to_csv(sample_output_folder/'features.tsv',sep="\t",index_label="")
        
            # Coordinates
            df_coordinates_tsv = adata.obs[['x_centroid','y_centroid']]
            df_coordinates_tsv.columns=[['x','y']]
            df_coordinates_tsv.index = adata.obs.index
            print(f"Coordinates: {df_coordinates_tsv.head()}")
            df_coordinates_tsv.to_csv(sample_output_folder/'coordinates.tsv',sep="\t",index_label="")
        
            # Matrix
            scipy.io.mmwrite(f"{sample_output_folder}/counts.mtx",scipy.sparse.coo_matrix(adata.X))  # MR: use sparse version
            #scipy.io.mmwrite(f"{sample_output_folder}/counts.mtx",adata.X)
    
            print(f'{sample} ready.')
    
    
    print(f"{'xenium-brain-SergioSalas'} ready.")

def main():
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(
        description="""Xenium mouse brain data, annotated by Sergio Marco Salas."""
    )
        
    # Add arguments for input and output folders
    parser.add_argument('-o','--out_dir', help="Output directory to write files to.",required=True)
    
    # Parse the command-line arguments
    args = parser.parse_args()
    print(args)
    print(args.out_dir)
    convert_data(args.out_dir)

if __name__ == '__main__':
    main()
