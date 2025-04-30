#!/usr/bin/env python

# Author_and_contribution: Jieran Sun; created script

import argparse

parser = argparse.ArgumentParser(description="Aggregate the results into a single tsv file. Notice that the folder structure should be input_folder/*/results_files")

parser.add_argument(
    "-i", "--input_folder", help="Input folder of method results.", required=True
)
parser.add_argument(
    "-o", "--out_file", help="Output file.", required=True
)

args = parser.parse_args()

output_file = args.out_file
input_folder = args.input_folder

from pathlib import Path
import pandas as pd

def get_combined_results(input_folder, target_file, output_file, prefix=None):

    """
    Aggregate the results from each folder into a single tsv file. The folder should have the structure of
    input_folder/*/results_files. The results_files should be tsv files with the same indexes. And for the 
    merged tsv, the column name will have its {folder name (i.e. *)} as the prefix. In the case where prefix is 
    defined, the prefix added to the column will be {folder name} without {prefix}.

    Args:
        input_folder: The folder containing the results of each method.
        target_file: The target file name inside each method folder.
        output_file: The output file name.
        prefix: The prefix of the folder name. If specified, the folder name will be used as the prefix
            for the columns of the dataframe.

    Returns:
        None, but write the output file into a csv in the respective path
    """
    input_folder = Path(input_folder)
    results_folders = [f for f in input_folder.iterdir() if f.is_dir() and not f.name.startswith((".", "_"))]

    if prefix is not None:
        results_folders = [f for f in results_folders if f.name.startswith(prefix)]

    combined_labels = []

    for result in results_folders:
        domain_file = result / target_file
        if domain_file.exists():  # Ensure the file exists
            folder_name = result.name
            if prefix is not None:
                folder_name = folder_name[len(prefix):]
            domain_df = pd.read_table(domain_file, sep="\t", index_col=0)
            domain_df = domain_df.add_prefix(f"{folder_name}_")
            combined_labels.append(domain_df)

    # Combine all the domain dataframes
    combined_df = pd.concat(combined_labels, axis=1, join="outer")

    # Write the combined dataframe to the output file
    combined_df.to_csv(output_file, sep="\t", index_label="")

# We assume the structure of {input_folder}/{methods}/{config}/{n_clusters}/results.tsv
# Generate all the potential folders
config_folders = [(method, config) 
                  for method in Path(input_folder).iterdir() if method.is_dir() and not method.name.startswith((".", "_"))
                  for config in method.iterdir() if config.is_dir() and config.name.startswith("config")
                ]

for _, config in config_folders:
    get_combined_results(input_folder=config, 
                         target_file="domains.tsv", 
                         output_file=config / "combined_nclusters.tsv",
                         prefix="cluster_")

for method, _ in config_folders:
    out_file = method / "combined_configs.tsv"
    if not out_file.exists():
        get_combined_results(input_folder=method, 
                             target_file="combined_nclusters.tsv", 
                             output_file=out_file,
                             prefix="config_")

# input folder is the final layer. It should be a sample folder
get_combined_results(input_folder=input_folder, 
                     target_file="combined_configs.tsv", 
                     output_file=output_file)