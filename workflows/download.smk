import os
from shared.functions import get_git_directory

# listed all the available datasets here
configfile: "example_configs/download_config.yaml"

print("Run Download Workflow")

# Attach the specific github directory here
GIT_DIR = get_git_directory(config)

# Leave only datasets
datasets = config.pop("datasets")

# Get all the dataset folder
def get_all_input(wildcards):
    all_folder = []
    for dataset in datasets:
        all_folder.append(config["results_dir"] + "/" + dataset)
    return all_folder

############## starting snakemake pipelines ##################

# Defining all output wanted from this snakemake
rule all:
    input:
        get_all_input,


rule download:
    output:
        dir=directory(config["results_dir"] + "/{dataset}"),
    conda:
        lambda wildcards: GIT_DIR + datasets[wildcards.dataset]["env"]
    params:
        script=lambda wildcards: GIT_DIR + datasets[wildcards.dataset]["script"],
    shell:
        "{params.script} -o {output.dir}"
