#!/usr/bin/env bash

# Create the stardust conda environment named stardust_env
# conda env create -f stardust.yml

# Activate the environment
# source activate stardust_env

# Install the required R packages

pip3 install torch==1.13.0 torchvision torchaudio --extra-index-url https://download.pytorch.org/whl/cpu --no-cache-dir #### CPU
pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv torch_geometric -f https://data.pyg.org/whl/torch-1.13.0+cpu.html --no-cache-dir ### CPU
