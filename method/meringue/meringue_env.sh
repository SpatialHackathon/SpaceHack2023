#!/bin/bash

# Create the MERINGUE conda environment named scmeb_env
conda env create -f meringue.yml

# Activate the environment
source activate meringue_env

# Install the required R packages
Rscript -e "remotes::install_github('JEFworks-Lab/MERINGUE')"


