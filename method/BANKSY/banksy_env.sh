#!/bin/bash

# Create the BANKSY conda environment named banksy_env
conda env create -f banksy.yml

# Activate the environment
conda activate banksy_env

# Install the required R packages
Rscript -e "remotes::install_github('prabhakarlab/Banksy@v0.1.5', dependencies = TRUE)"



