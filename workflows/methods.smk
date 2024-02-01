import os

from shared.functions import get_git_directory, get_ncluster, get_sample_dirs


configfile: "path_configs/methods.yaml"


GIT_DIR = get_git_directory(config)
SEED = config["seed"]

TECHNOLOGY = config["technology"]

methods = config["methods"]


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


def get_requirements(wildcards):
    if wildcards.method == "BANKSY":
        return "BANKSY_requirements.info"
    if wildcards.method == "DRSC":
        return "DRSC_requirements.info"
    if wildcards.method == "maple":
        return "maple_requirements.info"
    if wildcards.method == "meringue":
        return "meringue_requirements.info"
    if wildcards.method == "precast":
        return "precast_requirements.info"
    if wildcards.method == "scMEB":
        return "scMEB_requirements.info"
    return []


rule BANKSY_requirements:
    output:
        temp("BANKSY_requirements.info"),
    conda:
        GIT_DIR + methods["BANKSY"]["env"]
    shell:
        """
        conda run Rscript -e \"remotes::install_github('prabhakarlab/Banksy', dependencies = TRUE, ref = 'b1a2c8bb2af06346f303637b9bba18faa1a1fe32')\"
        touch BANKSY_requirements.info
        """


rule DRSC_requirements:
    output:
        temp("DRSC_requirements.info"),
    conda:
        GIT_DIR + methods["DRSC"]["env"]
    shell:
        """
        conda run Rscript -e \"remotes::install_version(package = 'DR.SC', version = '3.3', repos = 'https://cran.uni-muenster.de/')\"
        touch DRSC_requirements.info
        """


rule maple_requirements:
    output:
        temp("maple_requirements.info"),
    conda:
        GIT_DIR + methods["maple"]["env"]
    shell:
        """
        conda run Rscript -e \"remotes::install_github('carter-allen/maple', ref = 'b173e89a7bc82c6ae09c7e0709d09ed22082172d')\"
        touch maple_requirements.info
        """


rule meringue_requirements:
    output:
        temp("meringue_requirements.info"),
    conda:
        GIT_DIR + methods["meringue"]["env"]
    shell:
        """
        conda run Rscript -e \"remotes::install_github('JEFworks-Lab/MERINGUE', ref = 'ca9e2ccabd95680d9ca0b323a8a507c038f2ea13')\"
        touch meringue_requirements.info
        """


rule precast_requirements:
    output:
        temp("precast_requirements.info"),
    conda:
        GIT_DIR + methods["precast"]["env"]
    shell:
        """
        conda run Rscript -e \"remotes::install_version(package = 'PRECAST', version = '1.6.3', repos = 'https://cran.uni-muenster.de/')\"
        touch precast_requirements.info
        """


rule scMEB_requirements:
    output:
        temp("scMEB_requirements.info"),
    conda:
        GIT_DIR + methods["scMEB"]["env"]
    shell:
        """
        conda run Rscript -e \"if(!require(SC.MEB)) install.packages('SC.MEB',repos = 'https://cran.r-project.org/')\"
        touch scMEB_requirements.info
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
        + "/{sample}/delaunay_traingulation/spatial_connectivities.mtx",
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
        + "/{sample}/delaunay_traingulation/spatial_connectivities.mtx",
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
