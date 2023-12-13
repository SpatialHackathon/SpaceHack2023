import os
from shared.functions import get_sample_dirs

configfile: "config.yaml"

GIT_DIR = os.getenv("GIT_DIR", "/mnt/hack_data/code/SpaceHack2023")

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

def generate_all_XXX(wildcards):
    return generate_metrics_results(config["data_dir"], "XXX", config["XXX"]["methods"], file_ext = "json")

rule all:
    input: generate_all_XXX

rule xxx:
    input:
        file = config["data_dir"] + "/{sample_method_config}/domains.tsv"
    output:
        file = config["data_dir"] + "/{sample_method_config}/results.json"
    shell:
        """
        
        """