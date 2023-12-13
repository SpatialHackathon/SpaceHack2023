#!/bin/bash

# Create the maple conda environment named maple_env
conda env create -f maple.yml

# Activate the environment
source activate maple_env

# Install the required R packages
conda run -n maple_env R -e "install.packages('maple')"