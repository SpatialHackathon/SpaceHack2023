import os

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


def generate_metrics_results(
    data_dir, metrics_name, methods, file_ext, configfiles=None
):
    result_files = []
    # sample directory
    for sample_dir in get_sample_dirs(data_dir):
        # method directory
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
                    if check_files_in_folder(
                        os.path.join(method_dir, dir_to_check), ["domains.tsv"]
                    ):
                        # Metric config directory
                        config_files = configfiles.keys() if configfiles else [""]
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
        if metric in config['configfiles'].keys():
            all_input += generate_metrics_results(
                config["data_dir"], metric, methods, file_ext="txt", configfiles=config['configfiles'][metric]
            )
        else:
            all_input += generate_metrics_results(
                config["data_dir"], metric, methods, file_ext="txt"
            )
    return all_input


rule all:
    input:
        generate_all_input,


def get_sample_labels(wildcards):
    samples_folder = os.path.join(config["data_dir"], wildcards.sample)
    if "labels.tsv" in os.listdir(samples_folder):
        return "-g " + os.path.join(samples_folder, "labels.tsv")
    else:
        return ""


def get_method_embedding(wildcards):
    method_config_folder = os.path.join(
        config["data_dir"], wildcards.sample, wildcards.method_config
    )
    if "embedding.tsv" in os.listdir(method_config_folder):
        return "-e " + os.path.join(method_config_folder, "embedding.tsv")
    else:
        return ""


rule metric:
    input:
        domains=config["data_dir"] + "/{sample}/{method_config}/domains.tsv",
    output:
        file=config["data_dir"] + "/{sample}/{method_config}/{metric}/results.txt",
    wildcard_constraints:
        sample="[a-zA-Z0-9_-]+",
        metric="[a-zA-Z0-9_-]+",
    conda:
        lambda wildcards: GIT_DIR + metrics[wildcards.metric]["env"]
    params:
        sample_labels=get_sample_labels,
        embeddings=get_method_embedding,
        script=lambda wildcards:GIT_DIR + metrics[wildcards.metric]["script"],
    shell:
        """
        {params.script} \
            -l {input.domains} \
            {params.sample_labels} \
            {params.embeddings} \
            -o {output.file}
        """
