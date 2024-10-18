#!/usr/bin/env bash

# Create the SpatialARI conda environment named SpatialARI_env
# conda env create -f SpatialARI.yml -n SpatialARI_env

# Activate the environment
# conda activate SpatialARI_env

# Install the required R packages
Rscript -e "remotes::install_github('RoseYuan/ClusteringMetrics@5691a9e', dependencies=TRUE)"

