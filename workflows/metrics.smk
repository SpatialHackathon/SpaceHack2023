import os
import json

from shared.functions import check_files_in_folder, get_git_directory, get_sample_dirs


# this specific pipeline setting
configfile: "example_configs/metrics_config.yaml"
# All methods and metrics available
configfile: "path_configs/metrics.yaml"
configfile: "path_configs/methods.yaml"


GIT_DIR = get_git_directory(config)

# Get all the methods and metrics that's being used
metrics = config["metrics"]
methods = list(config["methods"].keys())
DATASET_DIR = config["dataset_dir"]


def generate_metrics_results(data_dir, metrics_name, methods, file_ext):
    # getting metrics optargs.json file
    with open(GIT_DIR + metrics[metrics_name]["optargs"], "r") as file:
        opt = json.load(file)

    result_files = []
    # sample directory
    for sample_dir in get_sample_dirs(data_dir):
        # Check if ground truth is needed
        if opt["groundtruth"] and "labels.tsv" not in os.listdir(sample_dir):
            continue

        # Check all method results
        for method in methods:
            method_dir = os.path.join(sample_dir, method)
            if os.path.exists(method_dir):
                dirs_to_check = (
                    [method_dir]
                    if check_files_in_folder(method_dir, ["domains.tsv"])
                    else os.listdir(method_dir)
                )
                # method config directory
                for dir_to_check in dirs_to_check:
                    # Check if embedding is needed
                    if opt["embedding"] and "embedding.tsv" not in os.listdir(
                        os.path.join(method_dir, dir_to_check)
                    ):
                        continue

                    # Check if results exist
                    if check_files_in_folder(
                        os.path.join(method_dir, dir_to_check), ["domains.tsv"]
                    ):

                        # Metric config directory
                        config_files = (
                            config["config_files"][metrics_name].keys()
                            if opt["config_file"]
                            else [""]
                        )

                        # Generating final metric results path
                        for config_file_name in config_files:
                            result_files.append(
                                os.path.join(
                                    method_dir,
                                    dir_to_check,
                                    metrics_name,
                                    config_file_name,
                                    "results." + file_ext,
                                )
                            )
    return result_files


def generate_all_input(wildcards):
    all_input = []
    for metric in config["use_metrics"]:
        for dataset in config["datasets"]:
            data_dir = DATASET_DIR + "/" + dataset
            all_input += generate_metrics_results(
                data_dir=data_dir,
                metrics_name=metric,
                methods=methods,
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
    with open(GIT_DIR + metrics[metric]["optargs"], "r") as file:
        opt = json.load(file)

    if opt["groundtruth"]:
        samples_folder = os.path.join(DATASET_DIR, wildcards.dataset, wildcards.sample)
        if "labels.tsv" not in os.listdir(samples_folder):
            raise Exception("wrong optargs file (groundtruth)")

        return "-g " + os.path.join(samples_folder, "labels.tsv")
    else:
        return ""


def get_method_embedding(wildcards):
    # getting metrics optargs.json file
    metric = get_metric(wildcards)
    with open(GIT_DIR + metrics[metric]["optargs"], "r") as file:
        opt = json.load(file)

    if opt["embedding"]:
        method_config_folder = os.path.join(
            DATASET_DIR, wildcards.dataset, wildcards.sample, wildcards.method_config
        )
        if "embedding.tsv" not in os.listdir(method_config_folder):
            raise Exception("wrong optargs file (embedding)!")

        return "-e " + os.path.join(method_config_folder, "embedding.tsv")
    else:
        return ""


def get_metric_config(wildcards):
    # getting metrics optargs.json file
    metric = get_metric(wildcards)
    with open(GIT_DIR + metrics[metric]["optargs"], "r") as file:
        opt = json.load(file)

    if opt["config_file"]:
        config_key = wildcards.metric_config[wildcards.metric_config.find("/") + 1 :]
        if len(config) == 0:
            raise Exception("Wrong optargs or no config folder found")
        return (
            "-c "
            + GIT_DIR
            + "metric/"
            + metric
            + "/"
            + config["config_files"][metric][config_key]
        )
    else:
        return ""


def get_sample_coordinate(wildcards):
    # getting metrics optargs.json file
    metric = get_metric(wildcards)
    with open(GIT_DIR + metrics[metric]["optargs"], "r") as file:
        opt = json.load(file)

    if "physical_coordinate" in opt.keys():
        if opt["physical_coordinate"]:
            return (
                "--coordinates "
                + DATASET_DIR
                + f"/{wildcards.dataset}/{wildcards.sample}/coordinates.tsv"
            )
        else:
            return ""
    else:
        return ""


rule metric:
    input:
        domains=DATASET_DIR + "/{dataset}/{sample}/{method_config}/domains.tsv",
        script=lambda wildcards: GIT_DIR + metrics[get_metric(wildcards)]["script"],
    output:
        file=DATASET_DIR 
        + "/{dataset}/{sample}/{method_config}/{metric_config}/results.txt",
    wildcard_constraints:
        dataset="[a-zA-Z0-9_-]+",
        sample="[a-zA-Z0-9_-]+",
        method_config="[a-zA-Z0-9_-]+(\/config_[a-zA-Z0-9_-]+)?",
        metric_config="[a-zA-Z0-9_-]+(\/config_[a-zA-Z0-9_-]+)?",
    conda:
        lambda wildcards: GIT_DIR + metrics[get_metric(wildcards)]["env"]
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
