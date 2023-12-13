import os
from shared.functions import get_sample_dirs, check_files_in_folder, generate_metrics_results

configfile: "config.yaml"

GIT_DIR = os.getenv("GIT_DIR", "/home/ubuntu/workspace/SpaceHack2023")


def generate_all_ARI(wildcards):
    return generate_metrics_results(config["data_dir"], "ARI", config["ARI"]["methods"], file_ext = "txt")

def generate_all_Completeness(wildcards):
    return generate_metrics_results(config["data_dir"], "Completeness", config["Completeness"]["methods"], file_ext = "txt")

def generate_all_V_measure(wildcards):
    return generate_metrics_results(config["data_dir"], "V_measure", config["V_measure"]["methods"], file_ext = "txt", configfiles = config["V_measure"]["configfiles"])


rule all:
    input: generate_all_ARI, generate_all_Completeness, generate_all_V_measure

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
        /usr/bin/env python {GIT_DIR}/metric/ARI/ARI.py \
            -l {input.domains} \
            -g {input.labels} \
            -o {output.file}
        """

#rule Calinski_Harabasz: braucht embeddings

rule Completeness:
    input:
        domains = config["data_dir"] + "/{sample}/{method_config}/domains.tsv",
        labels = config["data_dir"] + "/{sample}/labels.tsv"
    output:
        file = config["data_dir"] + "/{sample}/{method_config}/Completeness/results.txt"
    conda:
        GIT_DIR + "/metric/Completeness/Completeness.yml"
    wildcard_constraints:
        sample="[a-zA-Z0-9_-]+"
    shell:
        """
        /usr/bin/env python {GIT_DIR}/metric/Completeness/Completeness.py \
            -l {input.domains} \
            -g {input.labels} \
            -o {output.file}
        """

rule V_measure:
    input:
        domains = config["data_dir"] + "/{sample}/{method_config}/domains.tsv",
        labels = config["data_dir"] + "/{sample}/labels.tsv"
    output:
        file = config["data_dir"] + "/{sample}/{method_config}/V_measure/{config_file_name}/results.txt"
    conda:
        GIT_DIR + "/metric/V_measure/V_measure.yml"
    wildcard_constraints:
        sample="[a-zA-Z0-9_-]+"
    params:
        configfile = lambda wildcards: config["V_measure"]["configfiles"][wildcards.config_file_name]
    shell:
        """
        /usr/bin/env python {GIT_DIR}/metric/V_measure/V_measure.py \
            -l {input.domains} \
            -g {input.labels} \
            -o {output.file} \
            --config {GIT_DIR}/metric/V_measure/config/{params.configfile}
        """
