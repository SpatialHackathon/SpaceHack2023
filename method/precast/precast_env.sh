#!/usr/bin/env bash

# Create the precast conda environment named drsc_env
# conda env create -f precast.yml

# Activate the environment
# conda activate precast_env

# Install the required R packages
Rscript -e "remotes::install_version(package = 'PRECAST', version = '1.6.3', repos = 'https://cran.uni-muenster.de/')"