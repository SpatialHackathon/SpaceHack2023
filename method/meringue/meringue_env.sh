#!/usr/bin/env bash

# Create the MERINGUE conda environment named scmeb_env
# conda env create -f meringue.yml

# Activate the environment
# source activate meringue_env

# Install the required R packages
Rscript -e "remotes::install_github('JEFworks-Lab/MERINGUE', ref = 'ca9e2ccabd95680d9ca0b323a8a507c038f2ea13')"


