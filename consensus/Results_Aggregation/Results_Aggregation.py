#!/usr/bin/env python

# Author_and_contribution: Jieran Sun; created script

import argparse

# TODO adjust description
parser = argparse.ArgumentParser(description="Aggregate the results into a single tsv file. Notice that the folder structure should be input_folder/*/results_files")

parser.add_argument(
    "-i", "--input_folder", help="Input folder of method results.", required=True
)
parser.add_argument(
    "-f", "--file_name", help="the name of the result files to be merged", required=True
)
parser.add_argument(
    "-p", "--prefix", help="the prefix of the folders to be merged, the prefix will be used to be removed for column name prefix", required=False
)
parser.add_argument(
    "-o", "--out_file", help="Output file.", required=True
)

args = parser.parse_args()

target_file = args.file_name
output_file = args.out_file
input_folder = args.input_folder

if args.prefix is not None:
    prefix = args.prefix

from pathlib import Path
import pandas as pd

def get_combined_results(input_folder, target_file, output_file, prefix=None):

    input_folder = Path(input_folder)
    results_folders = [f for f in input_folder.iterdir() if f.is_dir() and not f.name.startswith([".", "_"])]

    if prefix is not None:
        results_folders = [f for f in results_folders if f.name.startswith(prefix)]

    combined_labels = []

    for result in results_folders:
        domain_file = result / target_file
        if domain_file.exists():  # Ensure the file exists
            folder_name = result.name
            if prefex is not None:
                folder_name = folder_name[len(prefix):]
            domain_df = pd.read_table(domain_file, sep="\t", index_col=0)
            domain_df = domain_df.add_prefix(f"{folder_name}_")
            combined_labels.append(domain_df)

    # Combine all the domain dataframes
    combined_df = pd.concat(combined_labels, axis=1, join="outer")

    # Write the combined dataframe to the output file
    combined_df.to_csv(output_file, sep="\t", index_label="")

# Generate all the potential folders
config_folders = [(method, config) 
                  for method in Path(input_folder).iterdir() if method.is_dir() and not method.name.startswith([".", "_"])
                  for config in method.iterdir() if config.is_dir() and config.name.startswith("config")
                ]



for _, config in config_folders:
    get_combined_results(input_folder=config, 
                         target_file="domains.tsv", 
                         output_file=config / "combined_nclusters.tsv",
                         prefix="cluster_")

for method, 