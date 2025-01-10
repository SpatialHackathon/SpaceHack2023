import os
import json
from pathlib import Path
import pandas as pd

from shared.functions import get_git_directory, get_ncluster, get_sample_dirs, get_combined_domai

configfile: "path_config.yaml"
configfile: "excute_config.yaml"
# TODO should I set the base clusterings in the result folder? 
#configfile: "base_clusterings.yaml"

GIT_DIR = Path(get_git_directory(config))
DATASET_DIR = Path(config["dataset_dir"])
SEED = config["seed"]
datasets_selected = config["datasets_selected"]
consensus_algorithms = config["consensus_algorithms"]
# base_clusters = config["base_clusterings"]

def get_bc(wildcards):


def create_input_all(wildcards):
    files = []
    for dataset in datasets_selected:
        data_dir = DATASET_DIR / dataset
        if not data_dir.is_dir():
            continue

        for sample in get_sample_dirs(data_dir):
            files += [f"{sample}/consensus_{al}_{n}.tsv" 
                      for al in consensus_algorithms
                      for n  in config["n_clusters"][dataset]]

    return files

rule all:
    input:
        create_input_all,

rule consensus_calling:
    input:
        file=DATASET_DIR / "{dataset}/{sample}/combined_methods.tsv"
        # TODO How to structure this tsv file for proper instruction? Also allow no n_clust specify
        base_clusterings=DATASET_DIR / "{dataset}/{sample}/base_clustering_selected.tsv"
        script=GIT_DIR / "consensus/Consensus_{algorithm}.r"
    output:
        file=DATASET_DIR / "{dataset}/{sample}/consensus_{algorithm}_{nclust}.tsv",
    wildcard_constraints:
        algorithm="[a-zA-Z_-]+",
        nclust="[0-9_-]+",
    params:
        seed=SEED,
        columns=lambda wildcards: config[wildcards.dataset][wildcards.sample][wildcards.nclust]["columns"]
        n_clust=lambda wildcards: f"--n_clust {wildcards.nclust}" if wildcards.algorithm=="weighted" else ""
        lambda_var= lambda wildcards: f"--lambda {config['lambda']}" if wildcards.algorithm=="weighted" else ""
        jar_file=lambda wildcards: f"--jar_file {GIT_DIR}/consensus/Consensus_weighted/networkanalysis-1.3.0.jar" if wildcards.algorithm=="weighted" else ""
    conda:
        lambda wildcards: str(GIT_DIR / f"consensus/Consensus_{wildcards.algorithm}.yaml")
    shell:
        """
        ulimit -s unlimited
        {input.script} \
            -i {input.file} \
            -o {output.file} \
            -c {params.columns} \
            --seed {params.seed} \
            {params.n_clust} \
            {params.lambda_var} \
            {params.jar_file} \
        """
