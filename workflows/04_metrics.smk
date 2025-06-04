import os
import json
from pathlib import Path
from shared.functions import check_files_in_folder, get_git_directory, get_sample_dirs

configfile: "path_config.yaml"
configfile: "excute_config.yaml"

GIT_DIR = get_git_directory(config)

# Get all the methods and metrics that's being used
METRICS = config["metrics"]
DATASET_DIR = Path(config["DATASET_DIR"])
datasets_selected= config["datasets_selected"]
methods_selected = config["methods_selected"]
metrics_selected = config["metrics_selected"]


def generate_metrics_results(data_dir, metric_name, methods, file_ext):
    # getting metrics optargs.json file
    with open(GIT_DIR + METRICS[metric_name]["optargs"], "r") as file:
        opt = json.load(file)

    result_files = []
    # sample directory
    for sample_dir in get_sample_dirs(data_dir):
        sample_dir = Path(sample_dir)
        # Check if ground truth is needed
        if opt["groundtruth"] and not (sample_dir / "labels.tsv").exists():
            continue

        method_paths = [m for m in sample_dir.iterdir() if m.name in methods]

        for method in method_paths:
            method_path = method
            if any(f.name.startswith("config") for f in method_path.iterdir()):
                res_fold = [
                    (config.name, ffolder) 
                    for config in method_path.iterdir() if config.is_dir() and config.name.startswith("config")
                    for ffolder in config.iterdir() if ffolder.is_dir() and ffolder.name.startswith("cluster")
                ]
            else:
                res_fold = [
                    ffolder for ffolder in method.iterdir() if ffolder.is_dir() and ffolder.name.startswith("cluster")
                ]
            if len(res_fold) > 0:
                for folders in res_fold:
                    res_folder = Path(folders[-1]) if isinstance(folders, (tuple, list)) else folders
                    if (res_folder / "domains.tsv").exists():
                        if opt["embedding"] and not (res_folder / "embedding.tsv").exists():
                            continue
                        if opt["config_file"]:
                            for c in config["config_files"][metric_name].keys():
                                result_files.append(res_folder / metric_name / c / f"results.{file_ext}")
                        else:
                            result_files.append(res_folder / metric_name / f"results.{file_ext}")

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

def get_metric(wildcards):
    # Trim metric_config if it has config path to it
    metric = wildcards.metric_config
    if "config" in metric:
        metric = metric[: metric.find("/")]

    return metric


def get_sample_labels(wildcards):
    # getting metrics optargs.json file
    metric = get_metric(wildcards)
    with open(GIT_DIR + METRICS[metric]["optargs"], "r") as file:
        opt = json.load(file)

    if opt["groundtruth"]:
        samples_folder = DATASET_DIR / wildcards.dataset / wildcards.sample
        if not (samples_folder / "labels.tsv").exists():
            raise Exception("wrong optargs file (groundtruth)")

        return f"-g {samples_folder / 'labels.tsv'}"
    else:
        return ""


def get_method_embedding(wildcards):
    # getting metrics optargs.json file
    metric = get_metric(wildcards)
    with open(GIT_DIR + METRICS[metric]["optargs"], "r") as file:
        opt = json.load(file)

    if opt["embedding"]:
        method_config_folder = DATASET_DIR / wildcards.dataset / wildcards.sample / wildcards.method_config
        if not (method_config_folder / "embedding.tsv").exists():
            raise Exception("wrong optargs file (embedding)!")

        return f"-e {method_config_folder / 'embedding.tsv'}"
    else:
        return ""


def get_metric_config(wildcards):
    # getting metrics optargs.json file
    metric = get_metric(wildcards)
    with open(GIT_DIR + METRICS[metric]["optargs"], "r") as file:
        opt = json.load(file)

    if opt["config_file"]:
        config_key = wildcards.metric_config[wildcards.metric_config.find("/") + 1 :]
        if len(config) == 0:
            raise Exception("Wrong optargs or no config folder found")
        return f"-c {Path(GIT_DIR) / 'metric' / metric / config['config_files'][metric][config_key]}"
    else:
        return ""


def get_sample_coordinate(wildcards):
    # getting metrics optargs.json file
    metric = get_metric(wildcards)
    with open(GIT_DIR + METRICS[metric]["optargs"], "r") as file:
        opt = json.load(file)

    if "physical_coordinate" in opt.keys():
        if opt["physical_coordinate"]:
            return f"--coordinates {DATASET_DIR / wildcards.dataset / wildcards.sample / 'coordinates.tsv'}"
        else:
            return ""
    else:
        return ""


rule metric:
    input:
        domains=DATASET_DIR / "{dataset}/{sample}/{method_config}/{nclust}/domains.tsv",
        script=lambda wildcards: GIT_DIR + METRICS[get_metric(wildcards)]["script"],
    output:
        file=DATASET_DIR /
        "{dataset}/{sample}/{method_config}/{nclust}/{metric_config}/results.txt",
    wildcard_constraints:
        dataset="[a-zA-Z0-9_-]+",
        sample="[a-zA-Z0-9_-]+",
        method_config="[a-zA-Z0-9_-]+(/config_[a-zA-Z0-9_-]+)?",
        metric_config="[a-zA-Z0-9_-]+(/config_[a-zA-Z0-9_-]+)?",
    conda:
        lambda wildcards: GIT_DIR + METRICS[get_metric(wildcards)]["env"]
    params:
        sample_labels=get_sample_labels,
        embeddings=get_method_embedding,
        config=get_metric_config,
        physical_coordinate=get_sample_coordinate,
    shell:
        """
        {input.script} \
            -l {input.domains} \
            {params.sample_labels} \
            {params.embeddings} \
            {params.config} \
            {params.physical_coordinate} \
            -o {output.file}
        """
