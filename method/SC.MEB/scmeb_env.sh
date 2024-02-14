#!/usr/bin/env bash

# Create the SC.MEB conda environment named scmeb_env
# conda env create -f SC.MEB.yml

# Activate the environment
# source activate scmeb_env

# Install the required R packages
# conda run -n scmeb_env R -e "install.packages('SC.MEB')"
Rscript -e "remotes::install_version(package = 'SC.MEB', version = '1.1', repos = 'https://cran.uni-muenster.de/')"

