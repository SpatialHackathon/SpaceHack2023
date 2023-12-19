import os
from shared.functions import get_sample_dirs, get_ncluster

configfile: "config.yaml"

GIT_DIR = os.getenv("GIT_DIR", "/home/ubuntu/workspace/SpaceHack2023")
TECHNOLOGY = config["technology"]
SEED = "42"

methods_info = dict()
methods_info["BANKSY"] = dict()
methods_info["BANKSY"]["script"] = "Rscript " + GIT_DIR + "/method/BANKSY/banksy.r"
methods_info["BANKSY"]["env"] = GIT_DIR + "/method/BANKSY/banksy.yml"
methods_info["BayesSpace"] = dict()
methods_info["BayesSpace"]["script"] = "Rscript " + GIT_DIR + "/method/BayesSpace/BayesSpace.r"
methods_info["BayesSpace"]["env"] = GIT_DIR + "/method/BayesSpace/BayesSpace.yml"
methods_info["DRSC"] = dict()
methods_info["DRSC"]["script"] = "Rscript " + GIT_DIR + "/method/DRSC/DRSC.r"
methods_info["DRSC"]["env"] = GIT_DIR + "/method/DRSC/DRSC.yml"
methods_info["GraphST"] = dict()
methods_info["GraphST"]["script"] = "python " + GIT_DIR + "/method/GraphST/method_GraphST.py"
methods_info["GraphST"]["env"] = GIT_DIR + "/method/GraphST/GraphST.yml"
methods_info["maple"] = dict()
methods_info["maple"]["script"] = "Rscript " + GIT_DIR + "/method/maple/maple.r"
methods_info["maple"]["env"] = GIT_DIR + "/method/maple/maple.yml"
methods_info["meringue"] = dict()
methods_info["meringue"]["script"] = "Rscript " + GIT_DIR + "/method/meringue/meringue.r"
methods_info["meringue"]["env"] = GIT_DIR + "/method/meringue/meringue.yml"
methods_info["precast"] = dict()
methods_info["precast"]["script"] = "Rscript " + GIT_DIR + "/method/precast/precast.r"
methods_info["precast"]["env"] = GIT_DIR + "/method/precast/precast.yml"
methods_info["scMEB"] = dict()
methods_info["scMEB"]["script"] = "Rscript " + GIT_DIR + "/method/SC.MEB/SC.MEB.r"
methods_info["scMEB"]["env"] = GIT_DIR + "/method/SC.MEB/SC.MEB.yml"
methods_info["SCAN_IT"] = dict()
methods_info["SCAN_IT"]["script"] = "python " + GIT_DIR + "/method/SCAN-IT/method_scanit.py"
methods_info["SCAN_IT"]["env"] = GIT_DIR + "/method/SCAN-IT/scanit.yml"
methods_info["scanpy"] = dict()
methods_info["scanpy"]["script"] = "python " + GIT_DIR + "/method/scanpy/method_scanpy.py"
methods_info["scanpy"]["env"] = GIT_DIR + "/method/scanpy/scanpy_env.yaml"
methods_info["SOTIP"] = dict()
methods_info["SOTIP"]["script"] = "python " + GIT_DIR + "/method/SOTIP/method_sotip.py"
methods_info["SOTIP"]["env"] = GIT_DIR + "/method/SOTIP/sotip.yml"
methods_info["SpaceFlow"] = dict()
methods_info["SpaceFlow"]["script"] = "python " + GIT_DIR + "/method/SpaceFlow/method_spaceflow.py"
methods_info["SpaceFlow"]["env"] = GIT_DIR + "/method/SpaceFlow/spaceflow.yml"
methods_info["spaGCN"] = dict()
methods_info["spaGCN"]["script"] = "python " + GIT_DIR + "/method/spaGCN/spaGCN.py"
methods_info["spaGCN"]["env"] = GIT_DIR + "/method/spaGCN/spaGCN.yml"
methods_info["STAGATE"] = dict()
methods_info["STAGATE"]["script"] = "python " + GIT_DIR + "/method/STAGATE/method_STAGATE.py"
methods_info["STAGATE"]["env"] = GIT_DIR + "/method/STAGATE/STAGATE.yml"


def create_input(method):
    input_files = []
    sample_dirs = get_sample_dirs(config["data_dir"])
    for sample_dir in sample_dirs:
        if method in config["config_files"].keys():
            for config_file_name in config["config_files"][method].keys():
                input_files.append(sample_dir + "/" + method + "/" + config_file_name + "/domains.tsv")
        else:
            input_files.append(sample_dir + "/" + method + "/domains.tsv")
    return input_files


def create_input_all(wildcards):
    files = []
    for method in config["methods"]:
        files += create_input(method)
    return files

rule all:
    input:
        create_input_all


def get_sample_image(wildcards):
    files = ["H_E.tiff", "H_E.png"]
    for file in files:
        image = config["data_dir"] + "/" + wildcards.sample + "/" + file
        if os.path.isfile(image):
            return "--image " + image
    return ""

def get_config_file(wildcards):
    current_config = GIT_DIR + "/method/" + wildcards.method + "/" + config["config_files"][wildcards.method][wildcards.config_file_name]
    return current_config

##########################################################
# requirements


def get_requirements(wildcards):
    if wildcards.method == "BANKSY": return "BANKSY_requirements.info"
    if wildcards.method == "DRSC": return "DRSC_requirements.info"
    if wildcards.method == "maple": return "maple_requirements.info"
    if wildcards.method == "meringue": return "meringue_requirements.info"
    if wildcards.method == "precast": return "precast_requirements.info"
    if wildcards.method == "scMEB": return "scMEB_requirements.info"
    return []

rule BANKSY_requirements:
    output:
        temp("BANKSY_requirements.info")
    conda:
        GIT_DIR + "/method/BANKSY/banksy.yml"
    shell:
        """
        conda run Rscript -e \"remotes::install_github('prabhakarlab/Banksy', dependencies = TRUE, ref = 'b1a2c8bb2af06346f303637b9bba18faa1a1fe32')\"
        touch BANKSY_requirements.info
        """

rule DRSC_requirements:
    output:
        temp("DRSC_requirements.info")
    conda:
        GIT_DIR + "/method/DRSC/DRSC.yml"
    shell:
        """
        conda run Rscript -e \"remotes::install_version(package = 'DR.SC', version = '3.3', repos = 'https://cran.uni-muenster.de/')\"
        touch DRSC_requirements.info
        """

rule maple_requirements:
    output:
        temp("maple_requirements.info")
    conda:
        GIT_DIR + "/method/maple/maple.yml"
    shell:
        """
        conda run Rscript -e \"remotes::install_github('carter-allen/maple', ref = 'b173e89a7bc82c6ae09c7e0709d09ed22082172d')\"
        touch maple_requirements.info
        """

rule meringue_requirements:
    output:
        temp("meringue_requirements.info")
    conda:
        GIT_DIR + "/method/meringue/meringue.yml"
    shell:
        """
        conda run Rscript -e \"remotes::install_github('JEFworks-Lab/MERINGUE', ref = 'ca9e2ccabd95680d9ca0b323a8a507c038f2ea13')\"
        touch meringue_requirements.info
        """

rule precast_requirements:
    output:
        temp("precast_requirements.info")
    conda:
        GIT_DIR + "/method/precast/precast.yml"
    shell:
        """
        conda run Rscript -e \"remotes::install_version(package = 'PRECAST', version = '1.6.3', repos = 'https://cran.uni-muenster.de/')\"
        touch precast_requirements.info
        """

rule scMEB_requirements:
    output:
        temp("scMEB_requirements.info")
    conda:
        GIT_DIR + "/method/SC.MEB/SC.MEB.yml"
    shell:
        """
        conda run Rscript -e \"if(!require(SC.MEB)) install.packages('SC.MEB',repos = 'https://cran.r-project.org/')\"
        touch scMEB_requirements.info
        """



##########################################################
# methods

rule method_with_config:
    input:
        coordinates = config["data_dir"] + "/{sample}/coordinates.tsv",
        matrix = config["data_dir"] + "/{sample}/log1p/counts.mtx",
        features = config["data_dir"] + "/{sample}/log1p/hvg/features.tsv",
        observations = config["data_dir"] + "/{sample}/observations.tsv",
        neighbors = config["data_dir"] + "/{sample}/delaunay_triangulation.mtx",
        dim_red = config["data_dir"] + "/{sample}/log1p/hvg/pca_20/dimensionality_reduction.tsv",
        requirements = get_requirements
    output:
        dir = directory(config["data_dir"] + "/{sample}/{method}/{config_file_name}"),
        file = config["data_dir"] + "/{sample}/{method}/{config_file_name}/domains.tsv",
    params:
        n_clusters = lambda wildcards: get_ncluster(config["data_dir"] + "/samples.tsv", wildcards.sample),
        technology = TECHNOLOGY,
        seed = SEED,
        configfile = get_config_file,
        image = get_sample_image,
        script = lambda wildcards: methods_info[wildcards.method]["script"]
    conda:
        lambda wildcards: methods_info[wildcards.method]["env"]
    wildcard_constraints:
        config_file_name="config_[a-zA-Z0-9_-]+"
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
        coordinates = config["data_dir"] + "/{sample}/coordinates.tsv",
        matrix = config["data_dir"] + "/{sample}/log1p/counts.mtx",
        features = config["data_dir"] + "/{sample}/log1p/hvg/features.tsv",
        observations = config["data_dir"] + "/{sample}/observations.tsv",
        neighbors = config["data_dir"] + "/{sample}/delaunay_triangulation.mtx",
        dim_red = config["data_dir"] + "/{sample}/log1p/hvg/pca_20/dimensionality_reduction.tsv",
        requirements = get_requirements
    output:
        dir = directory(config["data_dir"] + "/{sample}/{method}"),
        file = config["data_dir"] + "/{sample}/{method}/domains.tsv",
    params:
        n_clusters = lambda wildcards: get_ncluster(config["data_dir"] + "/samples.tsv", wildcards.sample),
        technology = TECHNOLOGY,
        seed = SEED,
        image = get_sample_image,
        script = lambda wildcards: methods_info[wildcards.method]["script"]
    conda:
        lambda wildcards: methods_info[wildcards.method]["env"]
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
