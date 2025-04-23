import os
import json
from pathlib import Path
import pandas as pd

from shared.functions import get_git_directory, get_ncluster, get_sample_dirs, get_combined_domain


# script specific setting
configfile: "example_configs/methods_config.yaml"
# All methods available
configfile: "path_configs/methods.yaml"

GIT_DIR = Path(get_git_directory(config))
DATASET_DIR = Path(config["dataset_dir"])
SEED = config["seed"]

methods = config.pop("methods")

# Find the technology of the datasets from their experiments.json
def get_technology(path):
    import json

    with open(path / "experiment.json", "r") as file:
        info = json.load(file)
    return info["technology"]

# Generates desired output based on no. of sample and config (output:domains.tsv)
def create_input(method, data_dir):
    from pathlib import Path
    input_files = []
    sample_dirs = get_sample_dirs(data_dir)
    dataset = data_dir.name

    for sample_dir in sample_dirs:

        sample_dir = Path(sample_dir)
        # Get cluster numbers for sweeping
        if dataset in config["n_clusters"].keys():
            n_clusters = config["n_clusters"][dataset]
        else:
            n_clusters = [get_ncluster(data_dir / "samples.tsv", sample_dir.name)]

        if method in config["config_files"].keys():
            for config_file_name in config["config_files"][method].keys():
                for n in n_clusters:
                    input_files.append(
                        sample_dir / method / config_file_name / f"cluster_{n}" / "domains.tsv"
                        )
                if len(n_clusters) > 1:
                    input_files.append(sample_dir / method / config_file_name / "combined_domains.tsv")
        else:
            for n in n_clusters:
                input_files.append(sample_dir / method / f"cluster_{n}" / "domains.tsv")
            if len(n_clusters) > 1:
                input_files.append(sample_dir / method / "combined_domains.tsv")
    return input_files


# For each method included, create all desirable outcome locations, because this function is
# defined on "use_methods" only, the script will only run the methods in that session in config file
def create_input_all(wildcards):
    files = []
    for dataset in config["datasets"]:
        data_dir = DATASET_DIR / dataset

        # Check if the dataset has been downloaded
        if not data_dir.is_dir():
            continue

        tech = get_technology(data_dir)
        for method in config["use_methods"]:
            if method in ["GraphST", "BayesSpace"] and tech != "Visium":
                continue
            files += create_input(method, data_dir)
    return files


rule all:
    input:
        create_input_all,


def get_sample_image(wildcards):
    # Using schema options:
    with open(GIT_DIR / methods[wildcards.method]["optargs"], "r") as file:
        opt = json.load(file)

    if opt["image"]:
        files = ["H_E.tiff", "H_E.png"]
        for file in files:
            image = DATASET_DIR / wildcards.dataset / wildcards.sample / file
            if os.path.isfile(image):
                return "--image " + image
            elif file == "H_E.png":
                return ""
    else:
        return ""


def get_config_file(wildcards):
    current_config = (
        GIT_DIR
        / "method"
        / wildcards.method
        / config["config_files"][wildcards.method][wildcards.config_file_name]
    )
    return current_config


##########################################################
# requirements


# Find if the method has an additional shell scripts for installation
def get_requirements(wildcards):
    if methods[wildcards.method].get("env_additional") is not None:
        return f"{wildcards.method}_requirements.info"
    else:
        return []


# if additional scripts are found, go through this process before generating the results
rule installation_requirements:
    params:
        install_script=lambda wildcards: GIT_DIR / methods[wildcards.method]["env_additional"],
    output:
        "{method}_requirements.info",
    conda:
        lambda wildcards: str(GIT_DIR / methods[wildcards.method]["env"])
    shell:
        """
        {params.install_script} && touch {output}
        """


##########################################################
# methods


# Get optargs options based on optargs files
def get_optargs(wildcards):
    with open(GIT_DIR / methods[wildcards.method]["optargs"], "r") as file:
        opt = json.load(file)
    return opt


# Get matrix in the input session
def get_matrix_input(wildcards):
    opt = get_optargs(wildcards)

    matrix_input = []
    # Find preprocessing steps
    match opt["matrix"]:
        case "counts":
            matrix_input = DATASET_DIR / wildcards.dataset / wildcards.sample / "qc/counts.mtx"
        case "transform":
            matrix_input = DATASET_DIR / wildcards.dataset / wildcards.sample / "log1p/counts.mtx"
        case "dimensionality_reduction":
            matrix_input = (
                DATASET_DIR / wildcards.dataset / wildcards.sample / "log1p/hvg/pca_35/dimensionality_reduction.tsv"
            )

    if matrix_input == []:
        raise (ValueError("no valid matrix option! Check your optargs.json file!"))

    return matrix_input


# Get features
def get_feature_input(wildcards):
    opt = get_optargs(wildcards)

    # feature input option
    if opt["integrated_feature_selection"]:
        feature_input = (
            DATASET_DIR /wildcards.dataset/wildcards.sample/"log1p/hvg/features.tsv"
        )
    else:
        feature_input = DATASET_DIR /wildcards.dataset/wildcards.sample/"qc/features.tsv"

    return feature_input


# Get neighbors
def get_neighbor_input(wildcards):
    opt = get_optargs(wildcards)

    neighbor_input = []
    # feature input option
    if opt["neighbors"]:
        neighbor_input = (
            DATASET_DIR /wildcards.dataset/wildcards.sample/"delaunay_triangulation/spatial_connectivities.mtx"
        )

    return neighbor_input


rule method_with_config:
    input:
        coordinates=DATASET_DIR / "{dataset}/{sample}/qc/coordinates.tsv",
        observations=DATASET_DIR / "{dataset}/{sample}/qc/observations.tsv",
        requirements=get_requirements,
        matrix=get_matrix_input,
        features=get_feature_input,
        neighbors=get_neighbor_input,
        script=lambda wildcards: GIT_DIR / methods[wildcards.method]["script"],
    output:
        dir=directory(DATASET_DIR / "{dataset}/{sample}/{method}/{config_file_name}/cluster_{n_clusters}"),
        file=DATASET_DIR / "{dataset}/{sample}/{method}/{config_file_name}/cluster_{n_clusters}/domains.tsv",
    params:
        matrix=lambda wildcards: (
            "-m "
            if get_optargs(wildcards)["matrix"] != "dimensionality_reduction"
            else "--dim_red "
        ),
        neighbors=lambda wildcards: "-n " if get_optargs(wildcards)["neighbors"] else "",
        technology=lambda wildcards: get_technology(DATASET_DIR / wildcards.dataset),
        seed=SEED,
        config_file=get_config_file,
        image=get_sample_image,
    conda:
        lambda wildcards: str(GIT_DIR / methods[wildcards.method]["env"])
    benchmark:
        DATASET_DIR / "{dataset}/{sample}/{method}/{config_file_name}/cluster_{n_clusters}/benchmark_method.txt"
    wildcard_constraints:
        config_file_name="config_[a-zA-Z0-9_-]+",
        n_clusters="\\d+"
    shell:
        """
        ulimit -s unlimited
        {input.script} \
            -c {input.coordinates} \
            {params.matrix}{input.matrix} \
            -f {input.features} \
            -o {input.observations} \
            -d {output.dir} \
            {params.image} \
            {params.neighbors}{input.neighbors} \
            --n_clusters {wildcards.n_clusters} \
            --technology "{params.technology}" \
            --seed {params.seed} \
            --config {params.config_file}
        """


rule method_without_config:
    input:
        coordinates=DATASET_DIR / "{dataset}/{sample}/qc/coordinates.tsv",
        observations=DATASET_DIR / "{dataset}/{sample}/qc/observations.tsv",
        requirements=get_requirements,
        matrix=get_matrix_input,
        features=get_feature_input,
        neighbors=get_neighbor_input,
    output:
        dir=directory(DATASET_DIR / "{dataset}/{sample}/{method}/cluster_{n_clusters}"),
        file=DATASET_DIR / "{dataset}/{sample}/{method}/cluster_{n_clusters}/domains.tsv",
    params:
        matrix=lambda wildcards: (
            "-m "
            if get_optargs(wildcards)["matrix"] != "dimensionality_reduction"
            else "--dim_red "
        ),
        neighbors=lambda wildcards: "-n " if get_optargs(wildcards)["neighbors"] else "",
        #n_clusters=lambda wildcards: get_ncluster(
        #    DATASET_DIR / wildcards.dataset / "samples.tsv", wildcards.sample
        #),
        technology=lambda wildcards: get_technology(DATASET_DIR / wildcards.dataset),
        seed=SEED,
        image=get_sample_image,
        script=lambda wildcards: GIT_DIR / methods[wildcards.method]["script"],
    conda:
        lambda wildcards: str(GIT_DIR / methods[wildcards.method]["env"])
    benchmark:
        DATASET_DIR / "{dataset}/{sample}/{method}/cluster_{n_clusters}/benchmark_method.txt"
    wildcard_constraints:
        method="[a-zA-Z0-9_-]+",
        n_clusters="\\d+"
    shell:
        """
        ulimit -s unlimited
        {params.script} \
            -c {input.coordinates} \
            {params.matrix}{input.matrix} \
            -f {input.features} \
            -o {input.observations} \
            -d {output.dir} \
            {params.image} \
            {params.neighbors}{input.neighbors} \
            --n_clusters {wildcards.n_clusters} \
            --technology {params.technology} \
            --seed {params.seed}
        """


rule get_combined_clusters:
    input:
        expand(DATASET_DIR / "{{dataset}}/{{sample}}/{{method}}/{{config_file_name}}/cluster_{n_clusters}/domains.tsv",
               n_clusters=lambda wildcards: config["n_clusters"][wildcards.dataset])
    output:
        combined_file=DATASET_DIR / "{dataset}/{sample}/{method}/{config_file_name}/combined_domains.tsv"
    run:
        # Get the directory path where cluster files are located
        results_path = Path(output.combined_file).parent
        
        # Find all cluster_* directories containing domains.tsv
        results = [f for f in results_path.iterdir() if f.name.startswith("cluster_")]
        
        combined_labels = []
        
        for result in results:
            domain_file = result / "domains.tsv"
            if domain_file.exists():  # Ensure the file exists
                folder_name = result.name
                domain_df = pd.read_table(domain_file, sep="\t", index_col=0)
                domain_df.columns = [folder_name]
                combined_labels.append(domain_df)
        
        # Combine all the domain dataframes
        combined_df = pd.concat(combined_labels, axis=1)
        
        # Write the combined dataframe to the output file
        combined_df.to_csv(output.combined_file, sep="\t", index_label="")