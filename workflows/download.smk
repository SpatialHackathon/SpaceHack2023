import os

configfile: "config.yaml"

GIT_DIR = os.getenv("GIT_DIR", "/mnt/hack_data/code/SpaceHack2023")

rule all:
    output:
        dir=directory(config['data_dir'])
    conda:
        GIT_DIR + "/data/" + config["env"]
    shell:
        "{GIT_DIR}/data/{config[script]} -o {output.dir}"
