#!/usr/bin/env bash

# Create the maple conda environment named maple_env
# conda env create -f maple.yml

# Activate the environment
# source activate maple_env

# Install the required R packages
Rscript -e "remotes::install_github('carter-allen/spruce', ref = '47b02300cc9a1d83213682bd78464115867d1763')"
Rscript -e "remotes::install_github('carter-allen/maple', ref = 'b173e89a7bc82c6ae09c7e0709d09ed22082172d')"