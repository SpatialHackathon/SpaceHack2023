#!/bin/bash

# Create the DR.SC conda environment named drsc_env
conda env create -f DRSC.yml

# Activate the environment
conda activate drsc_env

# Install the required R packages
Rscript -e "remotes::install_version(package = 'DR.SC', version = '3.3', repos = 'https://cran.uni-muenster.de/')"