#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Brian Long; wrote code in file

import argparse
from pathlib import Path
import json
import numpy as np
import pandas as pd
import anndata as ad
from scipy.io import mmwrite

import tempfile
import boto3
from botocore import UNSIGNED
from botocore.client import Config





# download URL: https://sea-ad-spatial-transcriptomics.s3.amazonaws.com/middle-temporal-gyrus/all_donors-h5ad/MTG_Version1.h5ad

def write_sample(path, sample, coordinates_df, observations_df, features_df,
                 counts, labels_df=None, img=None):
    
    sample_path = Path(path) / sample
    sample_path.mkdir(parents=True, exist_ok=True) 

    coordinates_df.to_csv(sample_path / 'coordinates.tsv', sep='\t', index_label='')
    features_df.to_csv(sample_path / 'features.tsv', sep='\t', index_label='')
    observations_df.to_csv(sample_path / 'observations.tsv', sep='\t', index_label='')
    mmwrite(sample_path / 'counts.mtx', counts)

    if labels_df is not None:
        labels_df.to_csv(sample_path / 'labels.tsv', sep='\t', index_label='')

    if img is not None:
        # TODO write to image_file
        # H_E.json must contain the scale
        pass

'''
--------------------------------------------------------------------------------
DOWNLOAD FUNCTIONS ---------------------------------------------------
--------------------------------------------------------------------------------
Function for downloading the SEA-AD MERSCOPE dataset from the AWS S3 
bucket.
'''
# AWS S3 bucket: https://sea-ad-spatial-transcriptomics.s3.amazonaws.com/middle-temporal-gyrus/all_donors-h5ad/MTG_Version1.h5ad
#https://sea-ad-single-cell-profiling.s3.amazonaws.com/MTG/MERFISH/SEAAD_MTG_MERFISH_all-nuclei.2023-05-08.h5ad
S3_BUCKET_NAME = 'sea-ad-single-cell-profiling'
S3_REGION_NAME = 'us-west-2'
S3_FILE_NAME = 'MTG/MERFISH/SEAAD_MTG_MERFISH_all-nuclei.2023-05-08.h5ad'
LOCAL_NAME = 'SEA_AD_MTG_Version1.h5ad'
def get_SEA_AD_h5ad():
    '''Load gene expression counts matrix as AnnData obj directly into 
    memory using tempfile.'''
    with tempfile.TemporaryDirectory() as temp_dir:
        # Initialize S3 client with unsigned credentials for public bucket
        s3_client = boto3.client('s3', region_name=S3_REGION_NAME,
                                 config=Config(signature_version=UNSIGNED))
        
        # Load one set of gene expr counts, raw or log2 (takes ~2min + ~3GB)
        counts_temp_path = str(Path(temp_dir).joinpath(LOCAL_NAME))
        s3_client.download_file(S3_BUCKET_NAME, S3_FILE_NAME, 
                                counts_temp_path)
        adata = ad.read_h5ad(counts_temp_path) #, backed='r')
        adata.obs["x"] =  adata.obsm["X_selected_cell_spatial_tiled"][:,0]
        adata.obs["y"] =  adata.obsm["X_selected_cell_spatial_tiled"][:,1]

    return adata

def subset_anndata(anndata_obj, obs_filter_key_values, combining_operation):
    '''
    Subset an AnnData object based on a dictionary of key-value pairs combined with a combining operation.
    returns a dict with the filter mask used and a view of AnnData object
    '''
    # Filter AnnData obj
    to_mask = [anndata_obj.obs[key].isin(obs_filter_key_values[key]).values for key in obs_filter_key_values]
    

    selection_mask = to_mask[0]

    for mask in to_mask:
        selection_mask = combining_operation(selection_mask, mask)
    
    return {"filter_mask" : selection_mask,
            "filtered_anndata" : anndata_obj[selection_mask,:]}








if __name__=='__main__':
    # parse output directory from user
    parser = argparse.ArgumentParser(
        description="Load data for SEA-AD middle temporal gyrus dataset from AWS S3 bucket and write to file."
    )
    parser.add_argument(
        "-o", "--out_dir", help="Output directory to write files to.", required=True
    )
    args = parser.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    
    ANNOTATED_LAYERS = ['L2/3', 'L5', 'L4', 'L6', 'L1']     

    
    # Load SEA_AD dataset
    sea_ad_anndata = get_SEA_AD_h5ad() 
    # select cells used in analysis:
    # inside selected pia-to-white matter rectangle and with layer annotations
    
    selection_dict = subset_anndata(sea_ad_anndata,
                      {"Used in analysis":[True], 
                       "Layer annotation":ANNOTATED_LAYERS},
                      np.logical_and)
    group_anndata = selection_dict["filtered_anndata"]

    COLUMNS_TO_SAVE = ['Section', 'Class', 'Subclass',
           'Supertype','Genes detected', 'Number of spots', 'Depth from pia',
           'Normalized depth from pia',
           'Layer annotation',"x","y"]
    SPATIAL_COLUMNS = ["x", "y"]
    technology = 'MERSCOPE'

    donors=[]
    samples =[]
    positions = []
    replicates = []
    directories = []
    n_clusters=[]


    for donor in list(group_anndata.obs["Donor ID"].unique()):
        print(donor)
        replicate_counter = 0
        for donor_section in list(
            group_anndata.obs.loc[group_anndata.obs["Donor ID"]==donor,:]["Section"].unique()):
            print(donor_section)
            #subset the input anndata object to just this section
            donor_section_results= subset_anndata(group_anndata,
                          {"Used in analysis":[True], 
                           "Layer annotation":ANNOTATED_LAYERS,
                          "Section":[donor_section]},
                          np.logical_and)                     
            adata = donor_section_results["filtered_anndata"]    
    
            features_df = adata.var.loc[:,:]
            observations_df = adata.obs.loc[
                                :, COLUMNS_TO_SAVE
                              ]
            coordinates_df = adata.obs.loc[:,SPATIAL_COLUMNS]

            counts = adata.X  # CSR sparse matrix
            labels_df_all = adata.obs.loc[:, ['Layer annotation']]
            labels_df_all.rename(columns={'Layer annotation':'label'}, 
                                 inplace=True)



            # Write dataset to file


            
            donors.append(donor)
            samples.append(donor_section)
            positions.append(np.nan)
            replicates.append(replicate_counter)
            replicate_counter+=1
            directories.append(donor_section)
            n_clusters.append(len(ANNOTATED_LAYERS))



            write_sample(out_dir, donor_section, coordinates_df, observations_df, features_df,
                     counts, labels_df=labels_df_all, img=None)


    samples_df = pd.DataFrame(data={'patient': donors, 
                                    'sample':samples, 
                                    'position':positions, 
                                    'replicate':replicates,
                                    'directory':directories, 
                                    'n_clusters':n_clusters
                                   })    


    ## Write metadata files
    samples_df.to_csv(out_dir / "samples.tsv", sep='\t', index_label='')

    with open(out_dir / 'experiment.json', 'w') as f:
        exp_info = {'technology': technology,
                    'species': 'Homo sapiens'}
        json.dump(exp_info, f)
        
        
       # Write license file -------------------------------------------------------
    license_text = (
        'Seattle Alzheimers Disease'+ 
        '(https://portal.brain-map.org/explore/seattle-alzheimers-disease)'+'\n'+ 
        'MERSCOPE v1 MTG Dataset is licensed under the Allen Institute terms of use'+'\n'+ 
        'See https://alleninstitute.org/citation-policy/ for the Allen Institute Citation'+'\n'+ 
        'Policy and https://alleninstitute.org/terms-of-use/ for the Allen Institute'+'\n'+ 
        'Terms of Use.'+'\n'
        )
    with open(out_dir / 'LICENSE.txt', "w") as f:
        f.write(license_text)