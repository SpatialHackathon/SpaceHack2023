#!/usr/bin/env bash

# Create the BASS conda environment named scmeb_env
# conda env create -f bass.yml -n bass_env

# Activate the environment
# source activate bass_env

# Install the required R packages
Rscript -e "remotes::install_github('xzhoulab/SPARK', ref = 'a8b4bf27b804604dfda53da42992f100b8e4e727', dependencies = TRUE)"
Rscript -e "remotes::install_github('zhengli09/BASS', ref = '37980c94a99f4b01ad5fa63555b3c5ab8af82b7e', dependencies = TRUE)"


