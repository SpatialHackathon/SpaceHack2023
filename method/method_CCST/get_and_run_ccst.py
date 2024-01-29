#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Liya Zaygerman, https://github.com/theinvisibleliya

import argparse

# global variables 
git_url = "https://github.com/xiaoyeye/CCST.git"  # The CCST repo
release_tag = "v1.0.1"  # Version used during SpaceHack2.0

parser = argparse.ArgumentParser(description="Method CCST, \
                                                https://github.com/xiaoyeye/CCST.git, \
                                                v1.0.1, \
                                                https://www.nature.com/articles/s43588-022-00266-5"
                                )

parser.add_argument(
    "-c", "--coordinates", help="Path to coordinates (as tsv).", required=True
)
parser.add_argument(
    "-m", "--matrix", help="Path to (transformed) counts (as mtx).", required=False
)
parser.add_argument(
    "-f", "--features", help="Path to features (as tsv).", required=True
)
parser.add_argument(
    "-o", "--observations", help="Path to observations (as tsv).", required=True
)
parser.add_argument(
    "-n",
    "--neighbors",
    help="Path to neighbor definitions. Square matrix (not necessarily symmetric) where each row contains the neighbors of this observation (as mtx).",
    required=False,
)
parser.add_argument("-d", "--out_dir", help="Output directory.", required=True)
parser.add_argument(
    "--dim_red",
    help="Reduced dimensionality representation (e.g. PCA).",
    required=False,
)
parser.add_argument("--image", help="Path to H&E staining.", required=False)
parser.add_argument(
    "--n_clusters", help="Number of clusters to return.", required=True, type=int
)
parser.add_argument(
    "--technology",
    help="The technology of the dataset (Visium, ST, imaging-based).",
    required=True,
)
parser.add_argument(
    "--seed", help="Seed to use for random operations.", required=True, type=int
)
parser.add_argument(
    "--config",
    help="Optional config file (json) used to pass additional parameters.",
    required=False,
)

args = parser.parse_args()

from pathlib import Path

out_dir = Path(args.out_dir)

# Output files
label_file = out_dir / "domains.tsv"
embedding_file = out_dir / "embedding.tsv"
# if additional output files are required write it also to out_dir

# Use these filepaths as input ...
coord_file = args.coordinates
feature_file = args.features
observation_file = args.observations

if args.neighbors is not None:
    neighbors_file = args.neighbors
if args.matrix is not None:
    matrix_file = args.matrix
if args.dim_red is not None:
    dimred_file = args.dim_red
if args.image is not None:
    image_file = args.image
if args.config is not None: # for CCST, we need the config file to add CCST parameters
    config_file = args.config

n_clusters = args.n_clusters
technology = args.technology
seed = args.seed

# TODO set the seed, if the method requires the seed elsewhere please pass it on
import random

random.seed(seed)
# np.random.seed(seed)
# torch.manual_seed(seed)

## Your code goes here
# we need the config to pass on the CCST-specific parameters
import json
try: 
    with open (args.config, "r") as c:
        config = json.load(c)
except:
    raise FileNotFoundError(f"config file not specified/not found")

# Since CCST is not really a package, but a command line tool, we generate a temporary directory and clone the CCST repo from github into it. Info: https://github.com/xiaoyeye/CCST
import tempfile
import subprocess
import os
import csv

def clone_and_process_repo(git_url, release_tag, output_dir):
    if config["data_type"] != "sc" and config["data_type"] != "nsc":
        raise ValueError(f"data_type must be specified in the .json config")
    with tempfile.TemporaryDirectory() as tmpdir:

        gitdir = Path(tmpdir) / "CCST"
        print(f"Created temporary directory at {tmpdir} to store the CCST repo")
        # Clone the repository
        subprocess.run(["git", "clone", git_url, gitdir], check=True)
        # Change to the repository directory
        subprocess.run(["git", "-C", gitdir, "checkout", release_tag], check=True)
        print(f"Got the CCST repo {release_tag}...")
        # get into the cloned CCST directory
        os.chdir(gitdir)

        # TODO: convert our files into what they want in the temp directory:
            # for sc
            # for non-sc

        """ CCST has 2 modes - "nsc" or "non-single-cell", for Visium
        and "sc" for Merfish. They write, that sc datasets other than
        Merfish can be used aswell. """

        if config["data_type"]=='nsc': # not single cell
            prepros_script_path = gitdir / "data_generation_ST.py"
            
            print("Running preprocessing for non single cell CCST...")
            command = ["python",prepros_script_path,
                       "--data_name",config["data_name"]]
            subprocess.run(command, check=True)
            # Idea 1: call their functions directly, omit run_CCST.py
            # from CCST_ST_utils import CCST_on_ST
            # args_ccst = # get from config
            # CCST_on_ST(args_ccst)
            
        elif config["data_type"]=='sc': # single-cell
            # TODO: convert the given tsv files into the right format for CCST
            prepros_script_path = gitdir / "data_generation_merfish.py"

            print("Running preprocessing for single cell CCST...")
            command = ["python", prepros_script_path,
                   "--data_name", config["data_name"]]
            # Idea 1: call their functions directly, omit run_CCST.py
            # from CCST_merfish_utils import CCST_on_MERFISH
            # args_ccst = # get from config
            # CCST_on_MERFISH(args_ccst)

        
        # # Idea 2: Run their offered run_CCST.py script, which separates
        # # "nsc" and "sc"
    
        # script_path = gitdir / "run_CCST.py"
        # # Run the Python script run_CCST.py from the CCST repo
        # print("Running the run_CCST script...")
        # command = ["python", script_path,
        #             "--data_path", tmpdir, # instead of data_path
        #            "--data_type", config["data_type"],
        #            "--data_name", config["data_name"],
        #            "--result_path", output_dir,
        #            "--embedding_data_path", output_dir]
    
        # subprocess.run(command, check=True)

    print('Temporary directory and file have been deleted.')

#######
    
# label_df = ...  # DataFrame with index (cell-id/barcode) and 1 column (label)
# embedding_df = None  # optional, DataFrame with index (cell-id/barcode) and n columns

######## Testing area, delete before the pull request #########
if __name__ == "__main__":
    out_dir.mkdir(parents=True, exist_ok=True)
    clone_and_process_repo(git_url, release_tag, out_dir)
    
##############################################################

## Write output
out_dir.mkdir(parents=True, exist_ok=True)

# TODO: save embeddings because CCST produces some - get them out of the temp directory
# label_df.columns = ["label"]
# label_df.to_csv(label_file, sep="\t", index_label="")
# if embedding_df is not None:
    # embedding_df.to_csv(embedding_file, sep="\t", index_label="")
