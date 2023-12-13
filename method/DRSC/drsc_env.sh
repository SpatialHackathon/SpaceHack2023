#!/bin/bash

# Create the DR.SC conda environment named drsc_env
conda env create -f DRSC.yml

# Activate the environment
conda activate drsc_env

# Install the required R packages
Rscript -e "remotes::install_github('feiyoung/DR.SC', ref = 'faac8fc57c0c36f8bcb9f125cafa7886f2f05413')"
