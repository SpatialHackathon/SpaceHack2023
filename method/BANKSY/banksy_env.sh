#!/usr/bin/env bash

# Create the BANKSY conda environment named banksy_env
# conda env create -f banksy.yml

# Activate the environment
# conda activate banksy_env

# Install the required R packages
Rscript -e "remotes::install_github('prabhakarlab/Banksy', dependencies = FALSE, ref = 'dbda6fde952e65f45409d9bca8e1f821746755cc')"

