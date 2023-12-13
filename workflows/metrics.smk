import os
from shared.functions import get_sample_dirs, check_files_in_folder

configfile: "config.yaml"

GIT_DIR = os.getenv("GIT_DIR", "/home/ubuntu/workspace/SpaceHack2023")

def generate_metrics_results(data_dir, metrics_name, methods, file_ext):
    result_files = []
    for sample_dir in get_sample_dirs(data_dir):
        for method in methods:
            method_dir = sample_dir + "/" + method
            if check_files_in_folder(method_dir, ["domains.tsv"]):
                result_files.append(method_dir + "/" + metrics_name + "/results." + file_ext)
            else:
                for config_dir in os.listdir(method_dir):
                    if check_files_in_folder(method_dir + "/" + config_dir, ["domains.tsv"]):
                        result_files.append(method_dir + "/" + config_dir + "/" + metrics_name + "/results." + file_ext)
    return(result_files)

def generate_all_ARI(wildcards):
    return generate_metrics_results(config["data_dir"], "ARI", config["ARI"]["methods"], file_ext = "txt")

rule all:
    input: generate_all_ARI

rule ARI:
    input:
        domains = config["data_dir"] + "/{sample}/{method_config}/domains.tsv",
        labels = config["data_dir"] + "/{sample}/labels.tsv"
    output:
        file = config["data_dir"] + "/{sample}/{method_config}/ARI/results.txt"
    conda:
        GIT_DIR + "/metric/ARI/ARI.yml"
    wildcard_constraints:
        sample="[a-zA-Z0-9_-]+"
    shell:
        """
        {GIT_DIR}/metric/ARI/ARI.py \
            -l {input.domains} \
            -g {input.labels} \
            -o {output.file}
        """
