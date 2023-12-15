#!/bin/bash

# Will only test 1 sample

# Please test env creation separately
technology=Visium # modify
n_clusters=5 # modify
data_dir=/home/jovyan/scratch/SpaceHack2/data/KPMP_KIDNEY # modify
out_dir=/home/jovyan/scratch/SpaceHack2/data/KPMP_KIDNEY/KPMP_test_out # output_dir, can be modified
git_repo=/home/jovyan/scratch/SpaceHack2/git_repo_clone
# git_repo=/home/jovyan/scratch/SpaceHack2/git_repo_clone
#git_repo=/home/jovyan/data/code/SpaceHack2023/

for sample_dir in $data_dir/*/ ; do
    
    sample=$(basename $sample_dir)
    sample_dir=${sample_dir%/}
    sample_out=$out_dir/$sample

    ## Methods
    source /opt/conda/bin/activate /home/jovyan/scratch/SpaceHack2/shared_envs/spagcn_test
    $git_repo/method/spaGCN/spaGCN.py \
    -c $sample_dir/coordinates.tsv \
    -m $sample_dir/counts.mtx  \
    -f $sample_dir/features.tsv \
    -o $sample_dir/observations.tsv \
    -d $out_dir/$sample/SpaGCN \
    --n_clusters $n_clusters \
    --technology $technology \
    --seed 42 \
    --config $git_repo/method/spaGCN/config/config_1.json

    #break
    
done
