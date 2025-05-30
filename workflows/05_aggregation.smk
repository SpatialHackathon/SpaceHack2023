import os
import json
from pathlib import Path
from shared.functions import check_files_in_folder, get_git_directory, get_sample_dirs

configfile: "path_config.yaml"
configfile: "excute_config.yaml"

GIT_DIR = Path(get_git_directory(config))

# Get all the methods and metrics that's being used
DATASET_DIR = Path(config["DATASET_DIR"])
datasets_selected= config["datasets_selected"]

def generate_input_files(data_dir):

    # Function to check for the existence of domains.tsv in a sample folder
    def has_domains(sample_path):
        return any(
            os.path.isfile(os.path.join(root, "domains.tsv"))
            for root, dirs, files in os.walk(sample_path)
        )

    result_files = []

    # sample directory
    aggregate_files = [
        f"{sample}/combined_methods.tsv"
        for sample in get_sample_dirs(data_dir)
        if any(
            os.path.isfile(os.path.join(root, "domains.tsv"))
            for root, dirs, files in os.walk(sample)
        )
    ]

    return aggregate_files


def generate_all_input(wildcards):
    all_input = []

    for dataset in datasets_selected:
        data_dir = DATASET_DIR / dataset
        if not data_dir.is_dir():
            continue

        all_input += generate_input_files(data_dir=data_dir)

    return all_input

rule all:
    input:
        generate_all_input,


rule aggregate_nclusters:
    input:
        results_folder=DATASET_DIR / "{dataset}" / "{sample}",
        script=GIT_DIR / "consensus" / "01_Results_Aggregation" / "Results_Aggregation.py",
    output:
        file=DATASET_DIR / "{dataset}" /  "{sample}" / "combined_methods.tsv",
    wildcard_constraints:
        dataset="[a-zA-Z0-9_-]+",
        sample="[a-zA-Z0-9_-]+",
    conda:
        GIT_DIR / "consensus" / "01_Results_Aggregation" / "Results_Aggregation.yaml"
    shell:
        """
        {input.script} \
            -i {input.results_folder} \
            -o {output.file} \
        """