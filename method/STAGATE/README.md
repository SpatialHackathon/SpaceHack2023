Test script:

#!/bin/bash

# Please test env creation separately
conda_env=/home/jovyan/conda_envs/test_STAGATE # modify
script_path=/home/jovyan/scratch/SpaceHack2/userfolders/pcai/SpaceHack2023/method/STAGATE/method_STAGATE.py # modify
config_path=/home/jovyan/scratch/SpaceHack2/userfolders/pcai/SpaceHack2023/method/STAGATE/config/config_1.json # modify (optional)
test_dir=/home/jovyan/scratch/SpaceHack2/method_results/LIBD_DLPFC/Br5292_151507 # output_dir, can be modified

test_data=/home/jovyan/scratch/SpaceHack2/data/LIBD_DLPFC/Br5292_151507
test_data_processed=/home/jovyan/scratch/SpaceHack2/method_results/LIBD_DLPFC/Br5292_151507


# uncomment optional lines as needed else delete them
source /opt/conda/bin/activate $conda_env
$script_path \
-c $test_data/coordinates.tsv \
-m $test_data/counts.mtx \
-f $test_data_processed/log1p/hvg/features.tsv \
-o $test_data/observations.tsv \
-d $test_dir/STAGATE \
--n_clusters 7 \
--config $config_path \
--seed 42 
