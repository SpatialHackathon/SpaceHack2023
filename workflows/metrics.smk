import os
from shared.functions import get_sample_dirs, check_files_in_folder

#configfile: "config.yaml"

GIT_DIR = os.getenv("GIT_DIR", "/home/ubuntu/workspace/SpaceHack2023")
methods = ["BANKSY", "BayesSpace", "DRSC", "meringue", "scMEB", "scanpy", "spaGCN", "STAGATE"]


def generate_metrics_results(data_dir, metrics_name, methods, file_ext, configfiles=None):
    result_files = []
    for sample_dir in get_sample_dirs(data_dir):
        for method in methods:
            method_dir = os.path.join(sample_dir, method)
            if os.path.exists(method_dir):
                dirs_to_check = [method_dir] if check_files_in_folder(method_dir, ["domains.tsv"]) else os.listdir(method_dir)
                for dir_to_check in dirs_to_check:
                    if check_files_in_folder(os.path.join(method_dir, dir_to_check), ["domains.tsv"]):
                        config_files = configfiles.keys() if configfiles else [""]
                        for config_file_name in config_files:
                            result_files.append(os.path.join(method_dir, dir_to_check, metrics_name, config_file_name, "results." + file_ext))
    return result_files


metrics_info = {}
metrics_info["ARI"] = {}
metrics_info["ARI"]["script"] = "python " + GIT_DIR + "/metric/ARI/ARI.py"
metrics_info["ARI"]["env"] = GIT_DIR + "/metric/ARI/ARI.yml"
metrics_info["Calinski-Harabasz"] = {}
metrics_info["Calinski-Harabasz"]["script"] = "python " + GIT_DIR + "/metric/Calinski-Harabasz/Calinski-Harabasz.py"
metrics_info["Calinski-Harabasz"]["env"] = GIT_DIR + "/metric/Calinski-Harabasz/Calinski-Harabasz.yml"
metrics_info["CHAOS"] = {}
metrics_info["CHAOS"]["script"] = "Rscript " + GIT_DIR + "/metric/CHAOS/CHAOS.r"
metrics_info["CHAOS"]["env"] = GIT_DIR + "/metric/CHAOS/CHAOS.yml"
metrics_info["cluster-specific-silhouette"] = {}
metrics_info["cluster-specific-silhouette"]["script"] = "Rscript " + GIT_DIR + "/metric/cluster-specific-silhouette/cluster-specific-silhouette.r"
metrics_info["cluster-specific-silhouette"]["env"] = GIT_DIR + "/metric/cluster-specific-silhouette/cluster-specific-silhouette.yml"
metrics_info["Completeness"] = {}
metrics_info["Completeness"]["script"] = "python " + GIT_DIR + "/metric/Completeness/Completeness.py"
metrics_info["Completeness"]["env"] = GIT_DIR + "/metric/Completeness/Completeness.yml"
metrics_info["Davies-Bouldin"] = {}
metrics_info["Davies-Bouldin"]["script"] = "python " + GIT_DIR + "/metric/Davies-Bouldin/Davies-Bouldin.py"
metrics_info["Davies-Bouldin"]["env"] = GIT_DIR + "/metric/Davies-Bouldin/Davies-Bouldin.yml"
metrics_info["domain-specific-f1"] = {}
metrics_info["domain-specific-f1"]["script"] = "Rscript " + GIT_DIR + "/metric/domain-specific-f1/domain-specific-f1.r"
metrics_info["domain-specific-f1"]["env"] = GIT_DIR + "/metric/domain-specific-f1/domain-specific-f1.yml"
metrics_info["Entropy"] = {}
metrics_info["Entropy"]["script"] = "python " + GIT_DIR + "/metric/Entropy/Entropy.py"
metrics_info["Entropy"]["env"] = GIT_DIR + "/metric/Entropy/Entropy.yml"
metrics_info["FMI"] = {}
metrics_info["FMI"]["script"] = "python " + GIT_DIR + "/metric/FMI/FMI.py"
metrics_info["FMI"]["env"] = GIT_DIR + "/metric/FMI/FMI.yml"
metrics_info["Homogeneity"] = {}
metrics_info["Homogeneity"]["script"] = "python " + GIT_DIR + "/metric/Homogeneity/Homogeneity.py"
metrics_info["Homogeneity"]["env"] = GIT_DIR + "/metric/Homogeneity/Homogeneity.yml"
metrics_info["jaccard"] = {}
metrics_info["jaccard"]["script"] = "python " + GIT_DIR + "/metric/jaccard/jaccard.py"
metrics_info["jaccard"]["env"] = GIT_DIR + "/metric/jaccard/jaccard.yaml"
#metrics_info["LISI"] = {}
#metrics_info["LISI"]["script"] = "python " + GIT_DIR + "/metric/LISI/LISI.py"
#metrics_info["LISI"]["env"] = GIT_DIR + "/metric/LISI/LISI.yml"
metrics_info["MCC"] = {}
metrics_info["MCC"]["script"] = "python " + GIT_DIR + "/metric/MCC/MCC.py"
metrics_info["MCC"]["env"] = GIT_DIR + "/metric/MCC/MCC.yaml"
metrics_info["NMI"] = {}
metrics_info["NMI"]["script"] = "Rscript " + GIT_DIR + "/metric/NMI/NMI.r"
metrics_info["NMI"]["env"] = GIT_DIR + "/metric/NMI/NMI.yml"
metrics_info["PAS"] = {}
metrics_info["PAS"]["script"] = "Rscript " + GIT_DIR + "/metric/PAS/PAS.r"
metrics_info["PAS"]["env"] = GIT_DIR + "/metric/PAS/PAS.yml"
#metrics_info["V_measure"] = {}
#metrics_info["V_measure"]["script"] = "" + GIT_DIR + "/metric/"
#metrics_info["V_measure"]["env"] = GIT_DIR + "/metric/"


def generate_all_input(wildcards):
    all_input = []
    for metrics in metrics_info.keys():
        all_input += generate_metrics_results(config["data_dir"], metrics, methods, file_ext = "txt")
    return all_input


rule all:
    input: generate_all_input

def get_sample_labels(wildcards):
    samples_folder = os.path.join(config["data_dir"], wildcards.sample)
    if ("labels.tsv" in os.listdir(samples_folder)): return "-g " + os.path.join(samples_folder, "labels.tsv")
    else: return ""

def get_method_embedding(wildcards):
    method_config_folder = os.path.join(config["data_dir"], wildcards.sample, wildcards.method_config)
    if ("embedding.tsv" in os.listdir(method_config_folder)): return "-e " + os.path.join(method_config_folder, "embedding.tsv")
    else: return ""


rule metrics:
    input:
        domains = config["data_dir"] + "/{sample}/{method_config}/domains.tsv"
    output:
        file = config["data_dir"] + "/{sample}/{method_config}/{metric}/results.txt"
    wildcard_constraints:
        sample="[a-zA-Z0-9_-]+",
        metric="[a-zA-Z0-9_-]+"
    conda:
        lambda wildcards: metrics_info[wildcards.metric]["env"]
    params:
        sample_labels = get_sample_labels,
        embeddings = get_method_embedding,
        script = lambda wildcards: metrics_info[wildcards.metric]["script"]
    shell:
        """
        {params.script} \
            -l {input.domains} \
            {params.sample_labels} \
            {params.embeddings} \
            -o {output.file}
        """
