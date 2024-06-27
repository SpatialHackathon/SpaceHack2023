#!/usr/bin/env python

import argparse
import os
import re
import shutil
import scanpy as sc
import anndata as an
import pandas as pd
import scipy
import tarfile
import git
import json
from urllib import request
from PIL import Image, ImageSequence


def data_retrieval(out):
    """
    general function with all steps
    """
    Image.MAX_IMAGE_PIXELS = 933120000
    
    path = os.path.join(out, 'tmp')
    os.makedirs(path) 

    print('loading data')
    
    # download hires images 
    request.urlretrieve('https://figshare.com/ndownloader/files/42985297', os.path.join(path,'D4.tif'))
    request.urlretrieve('https://figshare.com/ndownloader/files/42985300', os.path.join(path,'D7.tif'))
    request.urlretrieve('https://figshare.com/ndownloader/files/42985303', os.path.join(path,'D10.tif'))
    request.urlretrieve('https://figshare.com/ndownloader/files/42985294', os.path.join(path,'D14.tif'))
    
    #download and untar archive from GEO
    
    request.urlretrieve('https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149457/suppl/GSE149457_RAW.tar', os.path.join(path,'GSE149457_RAW.tar'))
    with tarfile.open(os.path.join(path,'GSE149457_RAW.tar')) as f:
        f.extractall(path)
        
    #download metadata from github
    
    git.Git(path).clone('https://github.com/madhavmantri/chicken_heart.git')
    
    # temporal folders rearrangements
    print('download is complete')
    
    samples = ['D4', 'D7', 'D10', 'D14']
    for sample in samples:
        os.makedirs(os.path.join(path,sample, "spatial")) 
    
    h5_files = [f for f in os.listdir(path) if re.match(r'.*spatial_RNAseq.*\.h5', f)]
    for sample in samples:
        name = next(obj for obj in h5_files if sample in obj)
        os.rename(os.path.join(path,name), os.path.join(path, sample, 'filtered_feature_bc_matrix.h5'))

    processed_path = os.path.join(path,'chicken_heart/data/chicken_heart_spatial_RNAseq_processed')
    lowres_files = [f for f in os.listdir(processed_path) if re.match(r'.*\.png', f)]
    for sample in samples:
        name = next(obj for obj in lowres_files if sample in obj)
        os.rename(os.path.join(processed_path,name), os.path.join(path, sample, 'spatial','tissue_lowres_image.png'))
    
    json_files = [f for f in os.listdir(processed_path) if re.match(r'.*\.json', f)]
    for sample in samples:
        name = next(obj for obj in json_files if sample in obj)
        os.rename(os.path.join(processed_path,name), os.path.join(path, sample, 'spatial','scalefactors_json.json'))
    
    csv_files = [f for f in os.listdir(processed_path) if re.match(r'.*\.csv', f)]
    for sample in samples:
        name = next(obj for obj in csv_files if sample in obj)
        os.rename(os.path.join(processed_path,name), os.path.join(path, sample, 'spatial', 'tissue_positions_list.csv'))
        
        
    # transform into png

    print('transforming tif into png, this will take some time')
    tif_files = [f for f in os.listdir(path) if re.match(r'.*\.tif', f)]
    
    for i in tif_files:
        im = Image.open(os.path.join(path,i))
        im.thumbnail(im.size)
        folder = i.split('.')[0]
        im.save(os.path.join(path, folder, 'spatial', 'tissue_hires_image.png'), "png")
    
    # anndata creation

    print('output writing')
    
    for sample in samples:
        res_path = os.path.join(out, f"sample_{sample}")
        os.makedirs(res_path) 
        adata = sc.read_visium(os.path.join(path, sample))
        adata.var_names_make_unique()
        
        # saving tables
        scipy.io.mmwrite(os.path.join(res_path,"counts.mtx"),adata.X)
        
        features = adata.var.copy()
        features.to_csv(os.path.join(res_path,"features.tsv"), sep="\t", index_label="")
        
        coord = pd.DataFrame(adata.obsm["spatial"], index=adata.obs_names, columns=["x","y"])
        coord.to_csv(os.path.join(res_path,"coordinates.tsv"), sep="\t", index_label="")
        
        
        meta = pd.read_csv(os.path.join(path,'chicken_heart/data/spatialRNAseq_metadata.csv'))
        meta['Unnamed: 0'] = meta['Unnamed: 0'].str.split('_').str[1]
        meta = meta.set_index('Unnamed: 0')
        meta = meta[meta['orig.ident']== sample]
        meta_label = meta[['region']]
        meta_label = meta_label.rename(columns={'region':'label'})
        meta_label.to_csv(os.path.join(res_path,"labels.tsv"), sep="\t", index_label="")

        obs = adata.obs.copy()
        obs.columns = ["in_tissue", "row", "col"]
        obs = obs.join(meta, how="left")
        obs.to_csv(os.path.join(res_path,"observations.tsv"), sep="\t", index_label="")

        shutil.copy2(os.path.join(path, sample, 'spatial','tissue_lowres_image.png'), os.path.join(res_path,'H_E_lowres.png'))
        shutil.copy2(os.path.join(path,f"{sample}.tif"), os.path.join(res_path,'H_E.tiff'))
        shutil.copy2(os.path.join(path, sample, 'spatial','scalefactors_json.json'), os.path.join(res_path,'images.json'))
    
    # common files
    with open(os.path.join(out,"experiment.json"), "w") as f:
        pl_info = {"technology": "Visium"}
        json.dump(pl_info, f)
    
    ## creation of metadata
    
    d = {"sample": ["D4", "D7", "D10", "D14"], "tissue": ["embryonic chicken heart", "embryonic chicken heart", "embryonic chicken heart", "embryonic chicken heart"], "time": ["day_4", "day_7", "day_10", "day_14"], 
         "capture_area": ["A1", "B1", "C1", "D1"], "n_sections": [5,4,2,1], "directory": ["sample_D4", "sample_D7", "sample_D10", "sample_D14"], "n_clusters": [6, 10, 11, 10]}
    df = pd.DataFrame(data=d)
    
    df.to_csv(os.path.join(out,"samples.tsv"), sep="\t", index_label="")
    
    # removing extra files
    print('cleaning')
    shutil.rmtree(path)
    print('DONE')

def main():
    
    parser = argparse.ArgumentParser(
        description="Load data for embryonic chicken heart dataset"
    )
        
    parser.add_argument('-o','--out_dir', help="Output directory to write files to.",required=True)
    
    args = parser.parse_args()
    data_retrieval(args.out_dir)

if __name__ == '__main__':
    main()
