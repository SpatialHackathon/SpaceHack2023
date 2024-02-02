#!/usr/bin/env bash

# Create the SC.MEB conda environment named scmeb_env
# conda env create -f SC.MEB.yml

# Activate the environment
# source activate scmeb_env

# Install the required R packages
Rscript -e "install.packages('SC.MEB')"

