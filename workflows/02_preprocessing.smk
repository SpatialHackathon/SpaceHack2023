import os
from shared.functions import check_files_in_folder, get_git_directory, get_sample_dirs

configfile: "path_config.yaml"
configfile: "excute_config.yaml"

GIT_DIR = get_git_directory(config)
DATASETS = config.pop("datasets")
datasets_selected = config.pop("datasets_selected")

# If all required input files are in the folder, generate the required output file for all sample folders
def create_input(file_list, input_file_name, data_dir):
    input_files = []
    for sample_dir in get_sample_dirs(data_dir):
        if check_files_in_folder(sample_dir, file_list):
            input_files.append(sample_dir + input_file_name)
    return input_files

def create_quality_control_input(wildcards):
    file_list = ["coordinates.tsv", "counts.mtx", "features.tsv", "observations.tsv"]
    all_qc_file = []
    for dataset in datasets_selected:
        data_dir = config["dataset_dir"] +  "/" + dataset
        if "experiment.json" in os.listdir(data_dir):
            all_qc_file += create_input(file_list, "/qc/counts.mtx", data_dir)
            all_qc_file += create_input(file_list, "/qc/features.tsv", data_dir)
            all_qc_file += create_input(file_list, "/qc/observations.tsv", data_dir)
            all_qc_file += create_input(file_list, "/qc/coordinates.tsv", data_dir)

    return all_qc_file

# Get the optargs.json file for QC input, if exists
def get_opt(wildcards):
    import json
    dataset = wildcards["dataset"]

    # default value
    opt = {"min_cells":1, "min_genes":1, "min_counts":1}

    # Check if customized value exist
    if "optargs" in DATASETS[dataset] and os.path.exists(GIT_DIR + DATASETS[dataset]["optargs"]):
        with open(GIT_DIR + DATASETS[dataset]["optargs"], "r") as file:
            opt_load = json.load(file)
            # Update opt values based on existing opt_load value
            opt.update({k: v for k, v in opt_load.items() if k in opt})

    return opt

####################### Preprocessing #######################
rule all:
    input:
        create_quality_control_input,

rule quality_control:
    input:
        coordinates=config["dataset_dir"] + "/{dataset}/{sample}/coordinates.tsv",
        matrix=config["dataset_dir"] + "/{dataset}/{sample}/counts.mtx",
        features=config["dataset_dir"] + "/{dataset}/{sample}/features.tsv",
        observations=config["dataset_dir"] + "/{dataset}/{sample}/observations.tsv",
    output:
        dir=directory(config["dataset_dir"] + "/{dataset}/{sample}/qc"),
        counts=config["dataset_dir"] + "/{dataset}/{sample}/qc/counts.mtx",
        features=config["dataset_dir"] + "/{dataset}/{sample}/qc/features.tsv",
        observations=config["dataset_dir"] + "/{dataset}/{sample}/qc/observations.tsv",
        coordinates=config["dataset_dir"] + "/{dataset}/{sample}/qc/coordinates.tsv",
    conda:
        GIT_DIR + "preprocessing/quality_control/qc_scanpy.yml"
    params:
        opt=lambda wildcards: get_opt(wildcards)
    shell:
        """
        python {GIT_DIR}preprocessing/quality_control/qc_scanpy.py \
          -c {input.coordinates} \
          -m {input.matrix} \
          -f {input.features} \
          -o {input.observations} \
          --min_genes {params.opt["min_genes"]}\
          --min_cells {params.opt["min_cells"]} \
          --min_counts {params.opt["min_counts"]} \
          -d {output.dir}
        """