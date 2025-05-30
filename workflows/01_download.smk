import os
from shared.functions import get_git_directory

# workflow specific setting
configfile: "path_config.yaml"
configfile: "excute_config.yaml"

# Attach the specific github directory here
GIT_DIR = get_git_directory(config)

# Leave only datasets
DATASETS = config.pop("datasets")
datasets_selected = config["datasets_selected"]

# Get all the dataset folder
def get_all_input(wildcards):
    all_folder = []
    for dataset in datasets_selected:
        all_folder.append(config["DATASET_DIR"] + "/" + dataset)
    return all_folder


############## starting snakemake pipelines ##################

# Defining all output wanted from this snakemake
rule all:
    input:
        get_all_input,

rule download:
    output:
        dir=directory(config["DATASET_DIR"] + "/{dataset}"),
    conda:
        lambda wildcards: GIT_DIR + DATASETS[wildcards.dataset]["env"]
    params:
        script=lambda wildcards: GIT_DIR + DATASETS[wildcards.dataset]["script"],
    shell:
        "{params.script} -o {output.dir}"
