import os
import json
from pathlib import Path
import pandas as pd

from shared.functions import get_git_directory, get_ncluster, get_sample_dirs

configfile: "path_config.yaml"
configfile: "excute_config.yaml"

GIT_DIR = Path(get_git_directory(config))
DATASET_DIR = Path(config["DATASET_DIR"])
SEED = config["SEED"]
datasets_selected = config["datasets_selected"]
consensus_algorithms = config["consensus_algorithms"]
n_clust_con = config.get("n_clust_consensus", {})
n_bcs = config["bc_numbers"]
bc_selection_metrics = config["selection_criteria"]
cme = config["cross_method_entropy"]

def create_input_all(wildcards):
    files = []
    for dataset in datasets_selected:
        data_dir = DATASET_DIR / dataset
        if not data_dir.is_dir():
            continue
        if dataset not in n_clust_con:
            if dataset in config.get("n_clusters", {}):
                n_clust_con[dataset] = config["n_clusters"][dataset]
            else:
                n_clust_con[dataset] = [get_ncluster(data_dir / "samples.tsv", sample_dir.name)]
        files += [f"{sample}/consensus/BC_{bc_num}/{b}/cluster_{n_clu}/consensus_{a}.tsv" 
                    for sample in get_sample_dirs(data_dir)
                    for a in consensus_algorithms
                    for b in bc_selection_metrics
                    for bc_num in n_bcs
                    for n_clu in n_clust_con[dataset]]
        if cme:
            ent_files = ["cross_meth_ent", "selected_bcs"]
            files += [f"{sample}/consensus/BC_{bc_num}/{b}/cluster_{n_clu}/{r}.tsv" 
                        for sample in get_sample_dirs(data_dir)
                        for b in bc_selection_metrics
                        for bc_num in n_bcs
                        for n_clu in n_clust_con[dataset]
                        for r in ent_files]
    return files

rule all:
    input:
        create_input_all,

rule consensus_calling:
    input:
        file=DATASET_DIR / "{dataset}/{sample}/combined_methods.tsv",
        base_clusterings=DATASET_DIR / "{dataset}/{sample}/consensus/base_clusterings/{bc_metrics}/BC_rankings.tsv",
        script=GIT_DIR / "consensus/03_Consensus_{algorithm}/Consensus_{algorithm}.r",
    output:
        file=DATASET_DIR / "{dataset}/{sample}/consensus/BC_{n_bcs}/{bc_metrics}/cluster_{n_clu}/consensus_{algorithm}.tsv",
    wildcard_constraints:
        algorithm="[a-zA-Z0-9_-]+",
        bc_metrics="[a-zA-Z0-9_-]+",
        n_clu="[0-9_-]+",
        n_bcs="[0-9_-]+",
    params:
        seed=SEED,
        lambda_var= lambda wildcards: f"--lambda {config['lambda']}" if wildcards.algorithm=="weighted" and config.get('lambda') is not None else "",
        jar_file=lambda wildcards: f"--jar_file {GIT_DIR}/consensus/03_Consensus_weighted/networkanalysis-1.3.0.jar" if wildcards.algorithm=="weighted" else "",
    conda:
        lambda wildcards: str(GIT_DIR / f"consensus/03_Consensus_{wildcards.algorithm}/Consensus_{wildcards.algorithm}.yaml")
    shell:
        """
        ulimit -s unlimited
        {input.script} \
            -i {input.file} \
            -o {output.file} \
            --seed {params.seed} \
            --base_clusterings {input.base_clusterings} \
            --n_clusters {wildcards.n_clu} \
            --n_bcs {wildcards.n_bcs} \
            {params.lambda_var} \
            {params.jar_file} \
        """

rule cross_method_entropy:
    input:
        file=DATASET_DIR / "{dataset}/{sample}/combined_methods.tsv",
        base_clusterings=DATASET_DIR / "{dataset}/{sample}/consensus/base_clusterings/{bc_metrics}/BC_rankings.tsv",
        script=GIT_DIR / "consensus/03_Cross_method_entropy/Cross_method_entropy.r",
    output:
        file=DATASET_DIR / "{dataset}/{sample}/consensus/BC_{n_bcs}/{bc_metrics}/cluster_{n_clu}/cross_meth_ent.tsv",
        selected_bc=DATASET_DIR / "{dataset}/{sample}/consensus/BC_{n_bcs}/{bc_metrics}/cluster_{n_clu}/selected_bcs.tsv",
    wildcard_constraints:
        bc_metrics="[a-zA-Z0-9_-]+",
        n_bcs="[0-9_-]+",
    conda:
        lambda wildcards: str(GIT_DIR / f"consensus/03_Cross_method_entropy/Cross_method_entropy.yaml")
    shell:
        """
        ulimit -s unlimited
        {input.script} \
            -i {input.file} \
            -o {output.file} \
            --BC_output {output.selected_bc} \
            --BC_ranking {input.base_clusterings} \
            --n_clusters {wildcards.n_clu} \
            --n_bcs {wildcards.n_bcs}
        """
