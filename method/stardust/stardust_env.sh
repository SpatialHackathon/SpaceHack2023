#!/usr/bin/env bash

# Create the stardust conda environment named stardust_env
# conda env create -f stardust.yml

# Activate the environment
# source activate stardust_env

# Install the required R packages
Rscript -e "remotes::install_github('InfOmics/stardust', ref = 'f1b541704d4b4189b4daf4132289a084253349d9')"