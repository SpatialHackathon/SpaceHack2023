#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Florian Heyl (heylf); created code

# In[]
import argparse
import os
import tempfile
import requests
import pandas as pd
import scipy
import json
import pandas as pd
import numpy as np

from scipy.io import mmwrite
from pathlib import Path

# In[]
# TODO adjust description
parser = argparse.ArgumentParser(description="Load data for ...")

parser.add_argument(
    "-o", "--out_dir", help="Output directory to write files to.", required=True
)

args = parser.parse_args()

out_dir = Path(args.out_dir)

# In[]
def download_data(url, destination_folder, file_name):
    print(f'[INFO] Downloading annotated data from {url} and put it into {destination_folder}...') 
    
    # Create the destination folder if it doesn't exist
    if not os.path.exists(destination_folder):
        os.makedirs(destination_folder)

    # Get the file name from the URL
    file_name = os.path.join(destination_folder, file_name)

    # Download the file
    response = requests.get(url)
    with open(file_name, 'wb') as file:
        file.write(response.content)

    print('...done')

# In[]

#def get_data(out_dir):
def get_data(out_dir):

    with tempfile.TemporaryDirectory() as tmpdir:
        print('[INFO] created temporary directory', tmpdir)
        print('[START] COMPOSING DATA')

        # Names and urls of the samples are so inconsistent that I have to list them manually
        samples = ['well7_5','well10','well09','well04','well06','well07','well03','well01OB',
                'well05','sagittal3','well01brain','well1_5','well2_5','sagittal1','spinalcord',
                'well11','sagittal2','well3_5','well08']
        
        n_cluster = []
        directories = []

        for sample in samples:
            print(f'[INFO] Get sample {sample}')

            if not os.path.exists(f'{out_dir}/{sample}'):
                os.makedirs(f'{out_dir}/{sample}')

            download_data(f'https://zenodo.org/records/8327576/files/{sample}_spatial.csv?download=1', 
                        f'{tmpdir}', f'{sample}_spatial.csv')
        
            sample_dir = f'{out_dir}/{sample}'
            directories.append(sample_dir)

            # Write out coordinates.tsv, labels.tsv and observations.tsv
            with open(f'{tmpdir}/{sample}_spatial.csv', 'r') as f_in, \
                open(f'{sample_dir}/labels.tsv', 'w') as f_out_labels,\
                open(f'{sample_dir}/coordinates.tsv', 'w') as f_out_coords, \
                open(f'{sample_dir}/observations.tsv', 'w') as f_out_obs :
                
                # skip first two lines
                headline_1 = f_in.readline()
                headline_2 = f_in.readline()
                
                f_out_obs.write(headline_1.replace(',','\t').replace('NAME', ''))

                clusters = []

                f_out_coords.write('\tx\ty\tz\n')
                f_out_labels.write('\tMain_molecular_cell_type\tSub_molecular_cell_type\tMain_molecular_tissue_region\
                                \tSub_molecular_tissue_region\tMolecular_spatial_cell_type\n')

                for l in f_in:
                    data = l.strip('\n').split(",")
                    f_out_labels.write(f'{data[0]}\t{data[4]}\t{data[5]}\t{data[6]}\t{data[7]}\t{data[8]}\n')
                    f_out_coords.write(f'{data[0]}\t{data[1]}\t{data[2]}\t{data[3]}\n')
                    f_out_obs.write(l.replace(',','\t'))
                    clusters.append(data[4])

            n_cluster.append(len(set(clusters)))

            download_data(f'https://zenodo.org/records/8327576/files/{sample}raw_expression_pd.csv?download=1', 
                        f'{tmpdir}', f'{sample}raw_expression_pd.csv')

            # Write counts.mtx
            df = pd.read_table(f'{tmpdir}/{sample}raw_expression_pd.csv', sep=',', index_col=0)
            features = df.index 
            df = df.transpose()
            mmwrite(f'{sample_dir}/counts.mtx', scipy.sparse.csr_matrix(df))

            # TODO no further information from the data for features.tsv
            # Write out features.tsv
            with open(f'{sample_dir}/features.tsv', 'w') as f_out_features :
                f_out_features.write('\tgene_version\n')
                for feature in features:
                    f_out_features.write(f'{feature}\tNA\n')

        ## Metadata files
        download_data(f'https://zenodo.org/records/8327576/files/metadata.csv?download=1', 
                    f'{tmpdir}', 'metadata.csv')  

        tmp_samples = []
        position = []
        patient = []

        # TODO data from mouse so no patient --> should patient be NA?
        # TODO position I assumed is the tissue type, which is here brain. Is that correct?
        with open(f'{tmpdir}/metadata.csv', 'r') as f_in:
            f_in.readline()
            f_in.readline()
            for l in f_in:
                data = l.strip('\n').split(",")
                if ( "_".join(data[0].split("_")[:-1]) not in tmp_samples):
                    tmp_samples.append("_".join(data[0].split("_")[:-1]))
                    position.append(data[7])
                    patient.append(data[2])

        sort_metadata = [tmp_samples.index(x) for x in samples]
        position = np.array(position)[sort_metadata]
        patient = np.array(patient)[sort_metadata]

        samples_df = pd.DataFrame({'patient': patient, 'sample': samples, 'position': position, 
                                'replicate': ['NA']*len(samples), 'directory': directories, 'n_clusters': n_cluster})
        samples_df.loc[
            :, ["patient", "sample", "position", "replicate", "directory", "n_clusters"]
        ].to_csv(f'{out_dir}/samples.tsv', sep="\t", index_label="")

        with open(f"{out_dir}/experiment.json", "w") as f:
            exp_info = {"technology": 'STARmap+'}
            json.dump(exp_info, f)

        print('[FINISH]')        

# In[]
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




