import os
import json

from shared.functions import check_files_in_folder, get_git_directory, get_sample_dirs

configfile: "path_config.yaml"
configfile: "excute_config.yaml"

GIT_DIR = get_git_directory(config)

# Get all the methods and metrics that's being used
DATASET_DIR = config["dataset_dir"]
datasets_selected= config["datasets_selected"]
methods_selected = config["methods_selected"]


def generate_input_files(data_dir, methods):
    from pathlib import Path

    result_files = []
    # sample directory
    for sample_dir in get_sample_dirs(data_dir):
        # Check if ground truth is needed
        method_paths = [m for m in os.scandir(Path(sample_dir)) if m.name in methods]

        for method in method_list:
            if any(f.name.startswith("config") for f in method.iterdir()):
                res_fold = [
                    (config.name, ffolder) 
                    for config in os.scandir(method) if config.is_dir() and config.name.startswith("config")
                    for ffolder in os.scandir(config) if ffolder.is_dir() and ffolder.name.startswith("cluster")
                ]
            else:
                res_fold = [
                    (ffolder) for ffolder in os.scandir(method) if ffolder.is_dir() and ffolder.name.startswith("cluster")
                ]
            if len(res_fold) > 0:
                for folders in res_fold:
                    res_folder = folders[len(folders)-1] #Extract the folder information
                    if os.file.exist(res_folder / "domains.tsv"):
                        result_files.append(res_folder / "domains.tsv")

    return result_files


def generate_all_input(wildcards):
    all_input = []

        
    for metric in metrics_selected:
        for dataset in datasets_selected:
            data_dir = DATASET_DIR / dataset
            if not data_dir.is_dir():
                continue

            all_input += generate_metrics_results(
                data_dir=data_dir,
                metric_name=metric,
                methods=methods_selected,
                file_ext="txt",
            )

    return all_input

rule all:
    input:
        generate_all_input,


rule aggregate_nclusters:
    input:
        results_folder=DATASET_DIR + "/{dataset}/{sample}/{method_config}",
        script=GIT_DIR + "/consensus/Results_Aggregation.py",
    output:
        file=DATASET_DIR + "/{dataset}/{sample}/{method_config}/combined_nclusters.tsv",
    wildcard_constraints:
        dataset="[a-zA-Z0-9_-]+",
        sample="[a-zA-Z0-9_-]+",
        method_config="[a-zA-Z0-9_-]+(\/config_[a-zA-Z0-9_-]+)?",
    conda:
        GIT_DIR + "/consensus/Results_Aggregation.yaml"
    params:
        prefix="-p cluster_",
        file_name="domains.tsv"
    shell:
        """
        {input.script} \
            -i {input.results_folder} \
            -f {params.file_name}\
            {params.prefix} \
            -o {output.file} \
        """

rule aggregate_configs:
    input:
        results_folder=DATASET_DIR + "/{dataset}/{sample}/{method}",
        script=GIT_DIR + "/consensus/Results_Aggregation.py",
    output:
        file=DATASET_DIR + "/{dataset}/{sample}/{method}/combined_configs.tsv",
    wildcard_constraints:
        dataset="[a-zA-Z0-9_-]+",
        sample="[a-zA-Z0-9_-]+",
        method="[a-zA-Z0-9_-]+",
    conda:
        GIT_DIR + "/consensus/Results_Aggregation.yaml"
    params:
        prefix="-p config_",
        file_name="combined_nclusters.tsv"
    shell:
        """
        {input.script} \
            -i {input.results_folder} \
            -f {params.file_name}\
            {params.prefix} \
            -o {output.file} \
        """

rule aggregate_configs:
    input:
        results_folder=DATASET_DIR + "/{dataset}/{sample}/{method}",
        script=GIT_DIR + "/consensus/Results_Aggregation.py",
    output:
        file=DATASET_DIR + "/{dataset}/{sample}/{method}/combined_configs.tsv",
    wildcard_constraints:
        dataset="[a-zA-Z0-9_-]+",
        sample="[a-zA-Z0-9_-]+",
        method="[a-zA-Z0-9_-]+",
    conda:
        GIT_DIR + "/consensus/Results_Aggregation.yaml"
    params:
        prefix="-p config_",
        file_name="combined_nclusters.tsv"
    shell:
        """
        {input.script} \
            -i {input.results_folder} \
            -f {params.file_name}\
            {params.prefix} \
            -o {output.file} \
        """