#!/bin/bash

# Create the BANKSY conda environment named banksy_env
conda env create -f banksy.yml

# Activate the environment
conda activate banksy_env

# Install the required R packages
Rscript -e "remotes::install_github('prabhakarlab/Banksy', dependencies = TRUE, ref = 'b1a2c8bb2af06346f303637b9bba18faa1a1fe32')"



