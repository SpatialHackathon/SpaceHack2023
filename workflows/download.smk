import os

from shared.functions import get_git_directory


configfile: "path_configs/datasets.yaml"


print("Run Download Workflow")

GIT_DIR = get_git_directory(config)

datasets = config.pop("datasets")


def get_all_input(wildcards):
    all_folder = []
    for dataset in datasets:
        all_folder.append(config["results_dir"] + "/" + dataset)
    return all_folder


rule all:
    input:
        get_all_input,


rule download:
    output:
        dir=directory(config["results_dir"] + "/{dataset}"),
    conda:
        lambda wildcards: GIT_DIR + "/" + datasets[wildcards.dataset]["env"]
    params:
        script=lambda wildcards: GIT_DIR + "/" + datasets[wildcards.dataset]["script"],
    shell:
        "{params.script} -o {output.dir}"
