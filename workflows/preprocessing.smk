import os
from shared.functions import check_files_in_folder, get_git_directory, get_sample_dirs

configfile: "example_configs/preprocessing_config.yaml"
configfile: "path_configs/datasets.yaml"
#config["data_dir"] = config["dataset_dir"] + config["use_datasets"][0]

GIT_DIR = get_git_directory(config)
DATASETS = config.pop("use_datasets")
ALL_DATASETS = config.pop("datasets")
NEIGHBORS_INFOS = config.pop("neighbors_infos")


# If all required input files are in the folder, generate the required output file for all sample folders
def create_input(file_list, input_file_name, data_dir):
    input_files = []
    for sample_dir in get_sample_dirs(data_dir):
        if check_files_in_folder(sample_dir, file_list):
            input_files.append(sample_dir + input_file_name)
    return input_files


def create_neighbors_input(wildcards):
    file_list = ["coordinates.tsv", "counts.mtx", "features.tsv", "observations.tsv"]
    all_neighbors_files = []
    for neighbors_method in NEIGHBORS_INFOS.keys():
        all_neighbors_files += create_input(
            file_list, "/" + neighbors_method + "/spatial_connectivities.mtx"
        )
    # print(all_neighbors_files)
    return all_neighbors_files


def create_quality_control_input(wildcards):
    file_list = ["coordinates.tsv", "counts.mtx", "features.tsv", "observations.tsv"]
    all_qc_file = []
    for dataset in DATASETS:
        data_dir = config["dataset_dir"] +  "/" + dataset
        if "experiment.json" in os.listdir(data_dir):
            all_qc_file += create_input(file_list, "/qc/counts.mtx", data_dir)
            all_qc_file += create_input(file_list, "/qc/features.tsv", data_dir)
            all_qc_file += create_input(file_list, "/qc/observations.tsv", data_dir)
            all_qc_file += create_input(file_list, "/qc/coordinates.tsv", data_dir)

    return all_qc_file


def create_transformation_log1p_input(wildcards):
    file_list = ["coordinates.tsv", "counts.mtx", "features.tsv", "observations.tsv"]
    return create_input(file_list, "/log1p/counts.mtx")


def create_selection_hvg_input(wildcards):
    file_list = ["coordinates.tsv", "features.tsv", "observations.tsv"]
    return create_input(file_list, "/log1p/hvg/features.tsv")


def create_dimensionality_reduction_pca_input(wildcards):
    file_list = ["coordinates.tsv", "features.tsv", "observations.tsv"]
    return create_input(
        file_list, "/log1p/hvg/pca_" + config["n_pcs"] + "/dimensionality_reduction.tsv"
    )


# Get the optargs.json file for QC input
def get_opt(wildcards):
    import json
    from pathlib import Path
    dataset = wildcards["dataset"]
    if "optargs" in ALL_DATASETS[dataset] and os.path.exists(GIT_DIR + ALL_DATASETS[dataset]["optargs"]):
        with open(GIT_DIR + ALL_DATASETS[dataset]["optargs"], "r") as file:
            opt = json.load(file)
    else:
        opt = {"min_cells":1, "min_genes":1, "min_counts":1}
    return opt


####################### Preprocessing #######################
rule all:
    input:
#        create_neighbors_input,
#        create_transformation_log1p_input,
#        create_selection_hvg_input,
#        create_dimensionality_reduction_pca_input,
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
        min_genes=lambda wildcards: (
            f"--min_genes {get_opt(wildcards)['min_genes']}"
            if "min_genes" in get_opt(wildcards).keys()
            else ""
        ),
        min_cells=lambda wildcards: (
            f"--min_cells {get_opt(wildcards)['min_cells']}"
            if "min_cells" in get_opt(wildcards).keys()
            else ""
        ),
        min_counts=lambda wildcards: (
            f"--min_counts {get_opt(wildcards)['min_counts']}"
            if "min_counts" in get_opt(wildcards).keys()
            else ""
        ),
    shell:
        """
        python {GIT_DIR}preprocessing/quality_control/qc_scanpy.py \
          -c {input.coordinates} \
          -m {input.matrix} \
          -f {input.features} \
          -o {input.observations} \
          {params.min_genes}\
          {params.min_cells} \
          {params.min_counts} \
          -d {output.dir}
        """


rule neighbors:
    input:
        coordinates=config["dataset_dir"] + "/{dataset}/{sample}/qc/coordinates.tsv",
        matrix=config["dataset_dir"] + "/{dataset}/{sample}/qc/counts.mtx",
        features=config["dataset_dir"] + "/{dataset}/{sample}/qc/features.tsv",
        observations=config["dataset_dir"] + "/{dataset}/{sample}/qc/observations.tsv",
    output:
        file=config["dataset_dir"]
        + "/{dataset}/{sample}/{neighbors_method}/spatial_connectivities.mtx",
        outdir=directory(config["dataset_dir"] + "/{dataset}/{sample}/{neighbors_method}"),
    conda:
        lambda wildcards: GIT_DIR + NEIGHBORS_INFOS[wildcards.neighbors_method]["env"]
    wildcard_constraints:
        neighbors_method="(?!qc$)[a-zA-Z0-9_-]+",
    params:
        script=lambda wildcards: GIT_DIR
        + NEIGHBORS_INFOS[wildcards.neighbors_method]["script"],
    shell:
        """
        python {params.script} \
          -c {input.coordinates} \
          -m {input.matrix} \
          -f {input.features} \
          -o {input.observations} \
          -d {output.outdir}
        """


rule transformation_log1p:
    input:
        coordinates=config["dataset_dir"] + "/{dataset}/{sample}/qc/coordinates.tsv",
        matrix=config["dataset_dir"] + "/{dataset}/{sample}/qc/counts.mtx",
        features=config["dataset_dir"] + "/{dataset}/{sample}/qc/features.tsv",
        observations=config["dataset_dir"] + "/{dataset}/{sample}/qc/observations.tsv",
    output:
        dir=directory(config["dataset_dir"] + "/{dataset}/{sample}/log1p"),
        file=config["dataset_dir"] + "/{dataset}/{sample}/log1p/counts.mtx",
    conda:
        GIT_DIR + "preprocessing/transformation/log1p.yml"
    shell:
        """
        python {GIT_DIR}preprocessing/transformation/log1p.py \
          -c {input.coordinates} \
          -m {input.matrix} \
          -f {input.features} \
          -o {input.observations} \
          -d {output.dir}
        """


rule selection_hvg:
    input:
        coordinates=config["dataset_dir"] + "/{dataset}/{sample}/qc/coordinates.tsv",
        matrix=config["dataset_dir"] + "/{dataset}/{sample}/log1p/counts.mtx",
        features=config["dataset_dir"] + "/{dataset}/{sample}/qc/features.tsv",
        observations=config["dataset_dir"] + "/{dataset}/{sample}/qc/observations.tsv",
    output:
        dir=directory(config["dataset_dir"] + "/{dataset}/{sample}/log1p/hvg"),
        file=config["dataset_dir"] + "/{dataset}/{sample}/log1p/hvg/features.tsv",
    conda:
        GIT_DIR + "preprocessing/feature_selection/highly_variable_genes_scanpy.yml"
    shell:
        """
        python {GIT_DIR}preprocessing/feature_selection/highly_variable_genes_scanpy.py \
          -c {input.coordinates} \
          -m {input.matrix} \
          -f {input.features} \
          -o {input.observations} \
          -d {output.dir}
        """


rule dimensionality_reduction_pca:
    input:
        coordinates=config["dataset_dir"] + "/{dataset}/{sample}/qc/coordinates.tsv",
        matrix=config["dataset_dir"] + "/{dataset}/{sample}/log1p/counts.mtx",
        features=config["dataset_dir"] + "/{dataset}/{sample}/log1p/hvg/features.tsv",
        observations=config["dataset_dir"] + "/{dataset}/{sample}/qc/observations.tsv",
    params:
        n_components=config["n_pcs"],
        seed=config["seed"]
    output:
        folder=directory(
            config["dataset_dir"] + "/{dataset}/{sample}/log1p/hvg/pca_" + config["n_pcs"]
        ),
        file=config["dataset_dir"]
        + "/{dataset}/{sample}/log1p/hvg/pca_"
        + config["n_pcs"]
        + "/dimensionality_reduction.tsv",
    conda:
        GIT_DIR + "preprocessing/dimensionality_reduction/PCA.yml"
    shell:
        """
        python {GIT_DIR}preprocessing/dimensionality_reduction/PCA.py \
          -c {input.coordinates} \
          -m {input.matrix} \
          -f {input.features} \
          -o {input.observations} \
          -n {params.n_components} \
          -d {output.folder} \
          --seed {params.seed}
        """
