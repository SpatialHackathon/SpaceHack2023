import os
from shared.functions import get_sample_dirs, check_files_in_folder

configfile: "config.yaml"

GIT_DIR = os.getenv("GIT_DIR", "/mnt/hack_data/code/SpaceHack2023")

def create_input(file_list, input_file_name):
    input_files = []
    for sample_dir in get_sample_dirs(config["data_dir"]):
        if check_files_in_folder(sample_dir, file_list):
            input_files.append(sample_dir + input_file_name)
    return input_files

def create_neighbors_delaunay_triangulation_input(wildcards):
    file_list = ["coordinates.tsv", "counts.mtx", "features.tsv", "observations.tsv"]
    return create_input(file_list, "/delaunay_triangulation.mtx")

def create_transformation_log1p_input(wildcards):
    file_list = ["coordinates.tsv", "counts.mtx", "features.tsv", "observations.tsv"]
    return create_input(file_list, "/log1p/counts.mtx")

def create_selection_hvg_input(wildcards):
    file_list = ["coordinates.tsv", "features.tsv", "observations.tsv"]
    return create_input(file_list, "/log1p/hvg/features.tsv")

def create_dimensionality_reduction_pca_input(wildcards):
    file_list = ["coordinates.tsv", "features.tsv", "observations.tsv"]
    return create_input(file_list, "/log1p/hvg/pca_" + config["n_pcs"] + "/dimensionality_reduction.tsv")

rule all:
    input:
        create_neighbors_delaunay_triangulation_input,
        create_transformation_log1p_input,
        create_selection_hvg_input,
        create_dimensionality_reduction_pca_input,

rule neighbors_delaunay_triangulation:
    input:
        coordinates = config["data_dir"] + "/{sample}/coordinates.tsv",
        matrix = config["data_dir"] + "/{sample}/counts.mtx",
        features = config["data_dir"] + "/{sample}/features.tsv",
        observations = config["data_dir"] + "/{sample}/observations.tsv",
    output:
        file = config["data_dir"] + "/{sample}/delaunay_triangulation.mtx",
    conda:
        GIT_DIR + "/preprocessing/neighbors/delaunay_triangulation.yml"
    shell:
        """
        {GIT_DIR}/preprocessing/neighbors/delaunay_triangulation.py \
          -c {input.coordinates} \
          -m {input.matrix} \
          -f {input.features} \
          -o {input.observations} \
          -d {output.file}
        """

rule transformation_log1p:
    input:
        coordinates = config["data_dir"] + "/{sample}/coordinates.tsv",
        matrix = config["data_dir"] + "/{sample}/counts.mtx",
        features = config["data_dir"] + "/{sample}/features.tsv",
        observations = config["data_dir"] + "/{sample}/observations.tsv",
    output:
        dir=directory(config["data_dir"] + "/{sample}/log1p"),
        file = config["data_dir"] + "/{sample}/log1p/counts.mtx",
    conda:
        GIT_DIR + "/preprocessing/transformation/log1p.yml"
    shell:
        """
        {GIT_DIR}/preprocessing/transformation/log1p.py \
          -c {input.coordinates} \
          -m {input.matrix} \
          -f {input.features} \
          -o {input.observations} \
          -d {output.dir}
        """

rule selection_hvg:
    input:
        coordinates = config["data_dir"] + "/{sample}/coordinates.tsv",
        matrix = config["data_dir"] + "/{sample}/log1p/counts.mtx",
        features = config["data_dir"] + "/{sample}/features.tsv",
        observations = config["data_dir"] + "/{sample}/observations.tsv",
    output:
        dir = directory(config["data_dir"] + "/{sample}/log1p/hvg"),
        file = config["data_dir"] + "/{sample}/log1p/hvg/features.tsv",
    conda:
        GIT_DIR + "/preprocessing/feature_selection/highly_variable_genes_scanpy.yml"
    shell:
        """
        {GIT_DIR}/preprocessing/feature_selection/highly_variable_genes_scanpy.py \
          -c {input.coordinates} \
          -m {input.matrix} \
          -f {input.features} \
          -o {input.observations} \
          -d {output.dir}
        """

rule dimensionality_reduction_pca:
    input:
        coordinates = config["data_dir"] + "/{sample}/coordinates.tsv",
        matrix = config["data_dir"] + "/{sample}/log1p/counts.mtx",
        features = config["data_dir"] + "/{sample}/log1p/hvg/features.tsv",
        observations = config["data_dir"] + "/{sample}/observations.tsv",
    params:
        n_components=config["n_pcs"],
    output:
        folder = directory(config["data_dir"] + "/{sample}/log1p/hvg/pca_" + config["n_pcs"]),
        file = config["data_dir"] + "/{sample}/log1p/hvg/pca_" + config["n_pcs"] + "/dimensionality_reduction.tsv",
    conda:
        GIT_DIR + "/preprocessing/dimensionality_reduction/PCA.yml"
    shell:
        """
        {GIT_DIR}/preprocessing/dimensionality_reduction/PCA.py \
          -c {input.coordinates} \
          -m {input.matrix} \
          -f {input.features} \
          -o {input.observations} \
          -n {params.n_components} \
          -d {output.folder}
        """

