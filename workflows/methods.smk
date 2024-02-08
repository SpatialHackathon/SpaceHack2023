import os

from shared.functions import get_git_directory, get_ncluster, get_sample_dirs

configfile: "example_configs/methods_config.yaml"


GIT_DIR = get_git_directory(config)
SEED = config["seed"]

methods = config.pop("methods")

# Find the technology of the datasets from their experiments.json
def get_technology(path):
    import json
    from pathlib import Path

    with open(Path(path) / "experiment.json", "r") as file:
        info = json.load(file)
    return info["technology"]


TECHNOLOGY = get_technology(config["data_dir"])


# Generates desired output based on no. of sample and config (output:domains.tsv)
def create_input(method):
    input_files = []
    sample_dirs = get_sample_dirs(config["data_dir"])
    for sample_dir in sample_dirs:
        if method in config["config_files"].keys():
            for config_file_name in config["config_files"][method].keys():
                input_files.append(
                    sample_dir + "/" + method + "/" + config_file_name + "/domains.tsv"
                )
        else:
            input_files.append(sample_dir + "/" + method + "/domains.tsv")
    return input_files


# For each method included, create all desirable outcome locations, because this function is
# defined on "use_methods" only, the script will only run the methods in that session in config file
def create_input_all(wildcards):
    files = []
    for method in config["use_methods"]:
        files += create_input(method)
    return files

rule all:
    input:
        create_input_all,


def get_sample_image(wildcards):
    files = ["H_E.tiff", "H_E.png"]
    for file in files:
        image = config["data_dir"] + "/" + wildcards.sample + "/" + file
        if os.path.isfile(image):
            return "--image " + image
    return ""


def get_config_file(wildcards):
    current_config = (
        GIT_DIR
        + "method/"
        + wildcards.method
        + "/"
        + config["config_files"][wildcards.method][wildcards.config_file_name]
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
        install_script=lambda wildcards: methods[wildcards.method]["env_additional"],
    output:
        temp("{method}_requirements.info"),
    conda:
        lambda wildcards: GIT_DIR + methods[wildcards.method]["env"]
    shell:
        """
        {params.install_script} && touch {output}
        """


##########################################################
# methods


rule method_with_config:
    input:
        coordinates=config["data_dir"] + "/{sample}/coordinates.tsv",
        matrix=config["data_dir"] + "/{sample}/log1p/counts.mtx",
        features=config["data_dir"] + "/{sample}/log1p/hvg/features.tsv",
        observations=config["data_dir"] + "/{sample}/observations.tsv",
        neighbors=config["data_dir"]
        + "/{sample}/delaunay_triangulation/spatial_connectivities.mtx",
        dim_red=config["data_dir"]
        + "/{sample}/log1p/hvg/pca_20/dimensionality_reduction.tsv",
        requirements=get_requirements,
    output:
        dir=directory(config["data_dir"] + "/{sample}/{method}/{config_file_name}"),
        file=config["data_dir"] + "/{sample}/{method}/{config_file_name}/domains.tsv",
    params:
        n_clusters=lambda wildcards: get_ncluster(
            config["data_dir"] + "/samples.tsv", wildcards.sample
        ),
        technology=TECHNOLOGY,
        seed=SEED,
        configfile=get_config_file,
        image=get_sample_image,
        script=lambda wildcards: GIT_DIR + methods[wildcards.method]["script"],
    conda:
        lambda wildcards: GIT_DIR + methods[wildcards.method]["env"]
    wildcard_constraints:
        config_file_name="config_[a-zA-Z0-9_-]+",
    shell:
        """
        {params.script} \
            -c {input.coordinates} \
            -m {input.matrix} \
            -f {input.features} \
            -o {input.observations} \
            -n {input.neighbors} \
            -d {output.dir} \
            {params.image} \
            --dim_red {input.dim_red} \
            --n_clusters {params.n_clusters} \
            --technology {params.technology} \
            --seed {params.seed} \
            --config {params.configfile}
        """


rule method_without_config:
    input:
        coordinates=config["data_dir"] + "/{sample}/coordinates.tsv",
        matrix=config["data_dir"] + "/{sample}/log1p/counts.mtx",
        features=config["data_dir"] + "/{sample}/log1p/hvg/features.tsv",
        observations=config["data_dir"] + "/{sample}/observations.tsv",
        neighbors=config["data_dir"]
        + "/{sample}/delaunay_triangulation/spatial_connectivities.mtx",
        dim_red=config["data_dir"]
        + "/{sample}/log1p/hvg/pca_20/dimensionality_reduction.tsv",
        requirements=get_requirements,
    output:
        dir=directory(config["data_dir"] + "/{sample}/{method}"),
        file=config["data_dir"] + "/{sample}/{method}/domains.tsv",
    params:
        n_clusters=lambda wildcards: get_ncluster(
            config["data_dir"] + "/samples.tsv", wildcards.sample
        ),
        technology=TECHNOLOGY,
        seed=SEED,
        image=get_sample_image,
        script=lambda wildcards: GIT_DIR + methods[wildcards.method]["script"],
    conda:
        lambda wildcards: GIT_DIR + methods[wildcards.method]["env"]
    wildcard_constraints:
        method="[a-zA-Z0-9_-]+",
    shell:
        """
        {params.script} \
            -c {input.coordinates} \
            -m {input.matrix} \
            -f {input.features} \
            -o {input.observations} \
            -n {input.neighbors} \
            -d {output.dir} \
            {params.image} \
            --dim_red {input.dim_red} \
            --n_clusters {params.n_clusters} \
            --technology {params.technology} \
            --seed {params.seed}
        """
