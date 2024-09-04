#!/usr/bin/env bash

# Create the spatialGE conda environment named spatialGE_env
# conda env create -f spatialGE.yml -n spatialGE_env

# Activate the environment
# conda activate spatialGE_env

# Install the required R packages
Rscript -e "remotes::install_github('FridleyLab/spatialGE@1.2.0.0000')"

