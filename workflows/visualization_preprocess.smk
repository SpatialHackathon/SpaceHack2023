import os

from shared.functions import check_files_in_folder, get_git_directory, get_sample_dirs

configfile: "path_configs/datasets.yaml"

GIT_DIR = "/home/jovyan/scratch/SpaceHack2/userfolders/jsun/workflow/SpaceHack2023/"
DATASET_DIR = "/home/jovyan/scratch/SpaceHack2/userfolders/jsun/workflow/data"
DATASETS = config.pop("datasets")

def create_visualization(data_dir):
    input_files = []

    file_list = ["coordinates.tsv", "counts.mtx", "features.tsv", "observations.tsv", "qc"]
    if all([check_files_in_folder(sample_dir, file_list) for sample_dir in get_sample_dirs(data_dir)]):
        return [data_dir + "/.visualization/pp_report.pdf"]

    return input_files

def generate_all_input(wildcards):
    input_list = []
    for data_dir in get_sample_dirs(DATASET_DIR):
        input_list += create_visualization(data_dir)
    return input_list

# For QC setting
def get_opt(wildcards):
    import json
    from pathlib import Path

    dataset = Path(config["data_dir"]).name

    with open(GIT_DIR + DATASETS[dataset]["optargs"], "r") as file:
        opt = json.load(file)
    return opt

def get_technology(wildcards):
    import json

    with open(DATASET_DIR + f"/{wildcards.dataset}/experiment.json", "r") as file:
        info = json.load(file)
    return info["technology"]

##################### start visualization #####################

rule all:
    input:
        generate_all_input

rule qc_visualization:
    input:
        coordinates=DATASET_DIR  + "/{dataset}/{sample}/coordinates.tsv",
        matrix=DATASET_DIR  + "/{dataset}/{sample}/counts.mtx",
        features=DATASET_DIR  + "/{dataset}/{sample}/features.tsv",
        observations=DATASET_DIR  + "/{dataset}/{sample}/observations.tsv",
        obs_qc=DATASET_DIR + "/{dataset}/{sample}/qc/observations.tsv",
        fea_qc=DATASET_DIR + "/{dataset}/{sample}/qc/features.tsv",
        opt= lambda wildcards: GIT_DIR + DATASETS[wildcards.dataset]["optargs"],
    output:
        dir=directory(DATASET_DIR + "/{dataset}/{sample}/visualization/" ),
        file=DATASET_DIR + "/{dataset}/{sample}/visualization/pp_report_sample.pdf"
    params:
        technology=get_technology
    conda:
        GIT_DIR + "preprocessing/visualization/visualization.yml"
    wildcard_constraints:
        sample="[a-zA-Z0-9_-]+",
        dataset="[a-zA-Z0-9_-]+",
    shell:
        """
        python {GIT_DIR}preprocessing/visualization/qc_visualization.py \
          -c {input.coordinates} \
          -m {input.matrix} \
          -f {input.features} \
          -o {input.observations} \
          --observation_qc {input.obs_qc} \
          --feature_qc {input.fea_qc} \
          --opt {input.opt} \
          --technology {params.technology} \
          -d {output.dir}
        """

rule merge_pdf:
    input:
        lambda wildcards: [sample_dir + "/visualization/pp_report_sample.pdf" for sample_dir in get_sample_dirs(DATASET_DIR + "/" + wildcards.dataset)]
    output:
        dir=directory(DATASET_DIR + "/{dataset}/.visualization/"),
        file=DATASET_DIR + "/{dataset}/.visualization/pp_report.pdf"
    conda:
        GIT_DIR + "preprocessing/visualization/visualization.yml"
    wildcard_constraints:
        dataset="[a-zA-Z0-9_-]+"
    shell:
        """
        python {GIT_DIR}preprocessing/visualization/pdf_merge.py \
        -d {output.dir} \
        -p {input}
        """

# Visualzation file will be saved in a hidden folder `.visualization` to avoid 
# conflict with downstream analysis. To see or download those files, use:
# `find . -type f -path "*/.visualization/pp_report.pdf" -print0 | tar -czvf pp_reports.tar.gz --null -T -` to create a tar.gz file.
