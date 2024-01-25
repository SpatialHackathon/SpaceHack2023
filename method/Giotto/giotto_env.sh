#!/bin/bash

# Create the Giotto conda environment named giotto_env
conda env create -f giotto.yml

# Activate the environment
conda activate giotto_env

# Install the required R packages
Rscript -e "remotes::install_github('drieslab/GiottoUtils', dependencies = FALSE, ref = '7c8f0010de6c916228834823455f48ed5b3fa706')"
Rscript -e "remotes::install_github('drieslab/GiottoClass', dependencies = FALSE, ref = 'fca6eb3f5ee6e8e7e9cfe8a0bb82721107f4872d')"
#Rscript -e "remotes::install_github('drieslab/GiottoData', dependencies = FALSE, ref = '50606245a01f151c6c308f3282f7b3fd87c67027')"
Rscript -e "remotes::install_version('colorRamp2', version = '0.1.0', repos = 'https://cloud.r-project.org')"
Rscript -e "remotes::install_github('drieslab/GiottoVisuals', dependencies = FALSE, ref = '8a68d8840ba4724b9a6cbc223dc7d6ef6f88f050')"
Rscript -e "remotes::install_github('drieslab/Giotto', dependencies = FALSE, ref = 'bd438bad459ad0be42788ed749e95e009bc1e32b')"
Rscript -e "remotes::install_bitbucket(repo = 'qzhudfci/smfishhmrf-r', ref='2ab4825')"



