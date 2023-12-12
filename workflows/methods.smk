import os
from shared.functions import get_sample_dirs

configfile: "config.yaml"

GIT_DIR = os.getenv("GIT_DIR", "/mnt/hack_data/code/SpaceHack2023")

def create_spaGCN_input(wildcards):
    input_files = []
    sample_dirs = get_sample_dirs(config["data_dir"])
    for sample_dir in sample_dirs:
        for config_file_name in config["config_files"]["spaGCN"].keys():
            input_files.append(sample_dir + "/spaGCN/" + config_file_name + "/domains.tsv")
    return input_files

def create_BayesSpace_input(wildcards):
    input_files = []
    sample_dirs = get_sample_dirs(config["data_dir"])
    for sample_dir in sample_dirs:
        input_files.append(sample_dir + "/BayesSpace/domains.tsv")
    return input_files

rule all:
    input: create_spaGCN_input, create_BayesSpace_input

rule method_spaGCN:
    input:
        coordinates = config["data_dir"] + "/{sample}/coordinates.tsv",
        matrix = config["data_dir"] + "/{sample}/log1p/counts.mtx",
        features = config["data_dir"] + "/{sample}/log1p/hvg/features.tsv",
        observations = config["data_dir"] + "/{sample}/observations.tsv",
        image = config["data_dir"] + "/{sample}/H_E.tiff",
    output:
        dir = directory(config["data_dir"] + "/{sample}/spaGCN/{config_file_name}"),
        file = config["data_dir"] + "/{sample}/spaGCN/{config_file_name}/domains.tsv",
    params:
        n_clusters = "7",
        technology = "Visium",
        seed = "42",
        configfile = lambda wildcards: config["config_files"]["spaGCN"][wildcards.config_file_name]
    conda:
        GIT_DIR + "/method/spaGCN/spaGCN.yml"
    shell:
        """
        {GIT_DIR}/method/spaGCN/spaGCN.py \
            -c {input.coordinates} \
            -m {input.matrix} \
            -f {input.features} \
            -o {input.observations} \
            -d {output.dir} \
            --image {input.image} \
            --n_clusters {params.n_clusters} \
            --technology {params.technology} \
            --seed {params.seed} \
            --config {GIT_DIR}/method/spaGCN/{params.configfile}
        """


rule method_BayesSpace:
    input:
        coordinates = config["data_dir"] + "/{sample}/coordinates.tsv",
        features = config["data_dir"] + "/{sample}/log1p/hvg/features.tsv",
        observations = config["data_dir"] + "/{sample}/observations.tsv",
        dim_red = config["data_dir"] + "/{sample}/log1p/hvg/pca_20/dimensionality_reduction.tsv",
    output:
        dir = directory(config["data_dir"] + "/{sample}/BayesSpace"),
        file = config["data_dir"] + "/{sample}/BayesSpace/domains.tsv",
    params:
        n_clusters = "7",
        technology = "Visium",
        seed = "42"
    conda:
        GIT_DIR + "/method/BayesSpace/BayesSpace.yml"
    shell:
        """
        {GIT_DIR}/method/BayesSpace/BayesSpace.r \
            -c {input.coordinates} \
            -f {input.features} \
            -o {input.observations} \
            -d {output.dir} \
            --dim_red {input.dim_red} \
            --n_clusters {params.n_clusters} \
            --technology {params.technology} \
            --seed {params.seed}
        """

