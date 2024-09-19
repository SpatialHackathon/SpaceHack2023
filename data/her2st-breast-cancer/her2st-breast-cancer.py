import argparse
import os
import shutil
from pathlib import Path
import json
import pandas as pd
import scipy
import tempfile
import requests

def download_large_file(url,target):
    # downloads large files from zenodo in chunks.
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(target, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192): 
                f.write(chunk)
            return target
            
def convert_data(out_dir):

    # This is a google drive link to an annotated data set, as per https://stagate.readthedocs.io/en/latest/T9_STARmap.html
    file_download_link = 'https://zenodo.org/api/records/3957257/files-archive'

    with tempfile.TemporaryDirectory() as temp_dir:

        temp_dir=Path(temp_dir)
        
        # download the folder 'temporary_file_host_spacehack'
        print(f'Downloading compressed data from zenodo into {temp_dir}...')
        download_large_file(file_download_link, temp_dir / 'archive.7z' )
        print('...done')

        print('Extracting files in command line:')
        os.system('7z e '+ str(temp_dir/'archive.7z')+' -o' + str(temp_dir/'extracted_root -aoa'))
    #     #extracts count-matrices.zip  images.zip  meta.zip  spot-selections.zip
        print(os.listdir(temp_dir))
        
        os.system('7z x '+ str(temp_dir/'extracted_root/meta.zip') +' -o' + str(temp_dir/'extracted_meta -pyUx44SzG6NdB32gY -aoa'))
        os.system('7z x '+ str(temp_dir/'extracted_root/spot-selections.zip') +' -o' + str(temp_dir/'extracted_spot-selections -pyUx44SzG6NdB32gY -aoa'))
        
        os.system('7z x '+ str(temp_dir/'extracted_root/images.zip') +' -o' + str(temp_dir/'extracted_images -pzNLXkYk3Q9znUseS -aoa'))
        os.system('7z x '+ str(temp_dir/'extracted_root/count-matrices.zip') +' -o' + str(temp_dir/'extracted_count-matrices -pzNLXkYk3Q9znUseS -aoa'))
        print('...done')

        print(os.listdir(temp_dir/'extracted_images'/'images'))

        # return

        samples = sorted([s[:2] for s in os.listdir(temp_dir/'extracted_meta')])

        print(f'Using annotated samples {samples}')
        
        out_dir = Path(out_dir)
        project_root = out_dir 

        os.makedirs(project_root, exist_ok=True)
        
        # create experiment.json
        experiment_json_dict = {"technology":"ST"}
        with open(project_root / 'experiment.json', 'w') as f:
            json.dump(experiment_json_dict, f)

        # create samples.tsv
        df_sample_tsv = pd.DataFrame(index = samples,)
        df_sample_tsv.to_csv(project_root/'samples.tsv',sep="\t",index_label="")
    
        
        for i,sample in enumerate(df_sample_tsv.index):
    
            print(f'Processing sample {i}: "{sample}"')

            os.system('gunzip '+ str(temp_dir/'extracted_spot-selections'/f'{sample}_selection.tsv.gz'))
    
            sample_output_folder = project_root / sample
    
            os.makedirs(sample_output_folder,exist_ok=True)

            df_coordinates = pd.read_csv(temp_dir /'extracted_meta'/f'{sample}_labeled_coordinates.tsv',sep='\t',index_col=0)

            df_selection = pd.read_csv(temp_dir /'extracted_spot-selections'/f'{sample}_selection.tsv',sep='\t',)
            df_selection.index = ((df_selection.x.astype(str)+'x'+df_selection.y.astype(str)))
            
            df_metadata = pd.merge(df_selection,df_coordinates,how='left',left_on=['new_x','new_y'],right_on=['x','y'],suffixes=['_selection','_coords'],)
            df_metadata.index=df_selection.index
            
            # Labels:
            df_labels_tsv = df_metadata[['label']].astype('category')
            print(f"Writing labels: {df_labels_tsv.head()}")
            df_labels_tsv.to_csv(sample_output_folder/'labels.tsv',sep="\t",index_label="")
    
            # Observations
            df_observations_tsv = df_metadata[['pixel_x_coords','pixel_y_coords','selected']].copy()
            df_observations_tsv.columns = ['x_pixel_on_H_E','y_pixel_on_H_E','selected']
            df_observations_tsv['selected'] = df_observations_tsv['selected'].astype('category')
            print(f"Writing observations: {df_observations_tsv.head()}")           
            df_observations_tsv.to_csv(sample_output_folder/'observations.tsv',sep="\t",index_label="")

    #         # Features

            print(os.listdir(temp_dir / 'extracted_count-matrices' / 'count-matrices'  ))
            os.system('gunzip '+ str(temp_dir/'extracted_count-matrices'/ 'count-matrices' / f'{sample}.tsv.gz'))
            df_counts = pd.read_csv(temp_dir / 'extracted_count-matrices' / 'count-matrices' / f'{sample}.tsv',sep='\t',index_col=0)
            df_counts = df_counts.reindex(df_metadata.index)
            df_features_tsv = df_counts.T[[]]
            print(f"Writing features: {df_features_tsv.head()}")
            df_features_tsv.to_csv(sample_output_folder/'features.tsv',sep="\t",index_label="")
        
            # Coordinates
            df_coordinates_tsv = df_coordinates[['x','y']]
            print(f"Coordinates: {df_coordinates_tsv.head()}")
            df_coordinates_tsv.to_csv(sample_output_folder/'coordinates.tsv',sep="\t",index_label="")
        
            # Matrix
            scipy.io.mmwrite(f"{sample_output_folder}/counts.mtx",df_counts.values)

            # Copy images:
            shutil.move(temp_dir/'extracted_images'/'images'/'HE'/f'{sample}.jpg',sample_output_folder/'H_E.jpg')
            
    
            print(f'{sample} ready.')
    
    print(f"{'her2st'} ready.")

def main():
    # Set up command-line argument parser
    parser = argparse.ArgumentParser(
        description="""An annotated STAGATE her2st breast cancer data set."""
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
