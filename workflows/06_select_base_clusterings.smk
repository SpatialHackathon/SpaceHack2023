import os
import json
from pathlib import Path
import pandas as pd

from shared.functions import get_git_directory, get_ncluster, get_sample_dirs

configfile: "path_config.yaml"
configfile: "excute_config.yaml"

GIT_DIR = Path(get_git_directory(config))
DATASET_DIR = Path(config["dataset_dir"])
SEED = config["seed"]
datasets_selected = config["dataset_selected"]
selection_metrics = config["selection_metrics"]

def create_input_all(wildcards):
    files = []
    for dataset in datasets_selected:
        data_dir = DATASET_DIR / dataset
        if not data_dir.is_dir():
            continue

        files += [f"{sample}/consensus/base_clusterings/{s}/BC_rankings.tsv"
                  for sample in get_sample_dirs(data_dir)
                  for s in selection_metrics if s != "Manual_selection"]
        if "Manual_selection" in selection_metrics:
            files += [f"{sample}/consensus/base_clusterings/Manual_selection" 
                      for sample in get_sample_dirs(data_dir)]

    return files

rule all:
    input:
        create_input_all,

rule Cross_method_ARI:
    input:
        label_file=DATASET_DIR / "{dataset}/{sample}/combined_methods.tsv",
        script=GIT_DIR / "consensus/Cross_method_ARI/Cross_method_ARI.r"
    output:
        file=DATASET_DIR / "{dataset}/{sample}/consensus/Cross_method_ARI.tsv",
    wildcard_constraints:
        sample="[a-zA-Z0-9_-]+",
    conda:
        lambda wildcards: str(GIT_DIR / "Cross_method_ARI/Cross_method_ARI.yaml")
    shell:
        """
        ulimit -s unlimited
        {input.script} \
            -i {input.file} \
            -o {output.file} 
        """

rule Smoothness_entropy:
    input:
        label_file=DATASET_DIR / "{dataset}/{sample}/combined_methods.tsv",
        coordinate_file=DATASET_DIR / "{dataset}/{sample}/coordinates.tsv",
        script=GIT_DIR / "consensus/Smoothness_entropy/Smoothness_entropy.r"
    output:
        file=DATASET_DIR / "{dataset}/{sample}/consensus/Smoothness_entropy.tsv",
    wildcard_constraints:
        sample="[a-zA-Z0-9_-]+",
    params:
        n_neighbors=config["n_neighbors"]
        seed=SEED
    conda:
        lambda wildcards: str(GIT_DIR / "Smoothness_entropy/Smoothness_entropy.yaml")
    shell:
        """
        ulimit -s unlimited
        {input.script} \
            -i {input.file} \
            -o {output.file} \
            -c {input.coordinate_file} \
            -n {params.neighbors} \
            -s {params.SEED}
        """

rule auto_rank_BCs:
    input:
        label_file=DATASET_DIR / "{dataset}/{sample}/combined_methods.tsv",
        selection_metrics=DATASET_DIR / "{dataset}/{sample}/consensus/{s_metrics}.tsv",
        script=GIT_DIR / "consensus/BC_ranking/BC_ranking.r"
    output:
        file=DATASET_DIR / "{dataset}/{sample}/consensus/base_clusterings/{s_metrics}/BC_rankings.tsv",
    wildcard_constraints:
        s_metrics="[a-zA-Z0-9_-]+",
    params:
        max_percentage=0.9 if "max_percentage" not in config.keys() else config["max_percentage"]
    conda:
        lambda wildcards: str(GIT_DIR / f"consensus/BC_ranking/BC_ranking.yaml")
    shell:
        """
        ulimit -s unlimited
        {input.script} \
            -i {input.file} \
            -o {output.file} \
            --selection_metrics {input.selection_metrics} \
            --max_percentage {params.max_percentage}
        """

rule manual_select_BCs:
    input:
        directory=DATASET_DIR / "{dataset}/{sample}"
    output:
        directory==DATASET_DIR / "{dataset}/{sample}/consensus/base_clusterings/Manual_selection"
    shell:
        """
        if [ -d '{input.directory}' ]; then
            mkdir -p '{output.directory}';
            echo "Directory {output.directory} has been created.";
        else
            echo "Directory {input.directory} does not exist. Exiting.";
            exit 1;
        fi
        """