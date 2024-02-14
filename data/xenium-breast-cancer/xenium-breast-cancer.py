#!/usr/bin/env python

# this was implemented by heylf

import argparse
import os
import sys
import tempfile
import requests
import shutil
import gzip
import pandas as pd

from zipfile import ZipFile

# pip install openpyxl

def download_data(url, destination_folder, boolzip):
    print(f'Downloading annotated data from {url} and put it into {destination_folder}...') 
    
    # Create the destination folder if it doesn't exist
    if not os.path.exists(destination_folder):
        os.makedirs(destination_folder)

    # Get the file name from the URL
    file_name = os.path.join(destination_folder, url.split("/")[-1])

    # Download the zip file
    response = requests.get(url)
    with open(file_name, 'wb') as file:
        file.write(response.content)

    if ( boolzip ):
        # Extract the contents of the zip file
        with ZipFile(file_name, 'r') as zip_ref:
            zip_ref.extractall(destination_folder)

    print('...done')

def gunzip_file(gzipped_file_path, output_file_path):
    with gzip.open(gzipped_file_path, 'rb') as f_in, open(output_file_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

def get_data(out_dir):
    with tempfile.TemporaryDirectory() as tmpdir:

        # download the folder 'temporary_file_host_spacehack'   
        download_data('https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip', 
                      f'{tmpdir}/replicate1', True)
        download_data('https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image.tif', 
                      f'{tmpdir}/replicate1', False)
        download_data('https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep2/Xenium_FFPE_Human_Breast_Cancer_Rep2_outs.zip', 
                      f'{tmpdir}/replicate2', True)
        download_data('https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep2/Xenium_FFPE_Human_Breast_Cancer_Rep2_he_image.tif', 
                      f'{tmpdir}/replicate2', False)
        download_data('https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep2/Xenium_FFPE_Human_Breast_Cancer_Rep2_gene_panel.json', 
                      f'{tmpdir}', False)
        download_data('https://cdn.10xgenomics.com/raw/upload/v1695234604/Xenium%20Preview%20Data/Cell_Barcode_Type_Matrices.xlsx', 
                      f'{tmpdir}', False)
    
        shutil.move(f'{tmpdir}/Xenium_FFPE_Human_Breast_Cancer_Rep2_gene_panel.json', f'{out_dir}/experiment.json')
    
        for replicate in ['replicate1', 'replicate2']:
            print(f'Extract data for {replicate}...')
            
            tif = ''
            if (replicate == 'replicate1'):
                tif = 'Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image'
            else:
                tif = 'Xenium_FFPE_Human_Breast_Cancer_Rep2_he_image'
    
            os.makedirs(f'{out_dir}/{replicate}')
            shutil.move(f'{tmpdir}/{replicate}/{tif}.tif', f'{out_dir}/{replicate}/{tif}.tif')
            shutil.move(f'{tmpdir}/{replicate}/outs/gene_panel.json', f'{out_dir}/{replicate}/experiment.json')    
            gunzip_file(f'{tmpdir}/{replicate}/outs/cell_feature_matrix/features.tsv.gz', f'{out_dir}/{replicate}/features.tsv')
            gunzip_file(f'{tmpdir}/{replicate}/outs/cell_feature_matrix/barcodes.tsv.gz', f'{out_dir}/{replicate}/observations.tsv')
            gunzip_file(f'{tmpdir}/{replicate}/outs/cell_feature_matrix/matrix.mtx.gz', f'{out_dir}/{replicate}/counts.mtx')
            gunzip_file(f'{tmpdir}/{replicate}/outs/cells.csv.gz', f'{out_dir}/{replicate}/coordinates.tsv')
            
            # Read the Excel file into a pandas DataFrame
            sheet = ''
            if (replicate == 'replicate1'):
                sheet = 'Xenium R1 Fig1-5 (supervised)'
            else:
                sheet = 'Xenium R2 Fig1-5 (supervised)'
    
            excel_file = f'{tmpdir}/Cell_Barcode_Type_Matrices.xlsx'
            
            df_labels = pd.read_excel(excel_file, sheet_name=sheet)  # Change 'Sheet1' to the sheet name you want to export
            df_labels.index = df_labels['Barcode']
            df_labels = df_labels.drop(columns='Barcode')
            df_labels.to_csv(f'{out_dir}/{replicate}/labels.tsv.csv', sep="\t", index_label="")
    
            print('...done')
    
def main():
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(
        description="""Xenium breast cancer data, annotated by 10x."""
    )

    # Add arguments for input and output folders
    parser.add_argument('-o','--out_dir', help="Output directory to write files to.",required=True)

    # Parse the command-line arguments
    args = parser.parse_args()
    print(args)
    print(args.out_dir)
    get_data(args.out_dir)

if __name__ == '__main__':
    main()