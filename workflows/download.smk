import os

configfile: "config.yaml"

rule all:
    output:
        dir=directory(config['data_dir'])
    conda:
        config["env"]
    shell:
        "{config[script]} -o {output.dir}"
