import os

configfile: "config.yaml"

GIT_DIR = os.getenv("GIT_DIR", "/home/ubuntu/workspace/SpaceHack2023")

rule all:
    input:
        dir=directory(config['results_dir'] + "/" + config['dataset'])

rule libd_dlpfc:
    output:
        dir=directory(config['results_dir'] + "/libd_dlpfc")
    conda:
        GIT_DIR + "/data/libd_dlpfc/libd_dlpfc.yml"
    shell:
        "Rscript {GIT_DIR}/data/libd_dlpfc/libd_dlpfc.r -o {output.dir}"

rule xenium_mouse_brain_SergioSalas:
    output:
        dir=directory(config['results_dir'] + "/xenium-mouse-brain-SergioSalas")
    conda:
        GIT_DIR + "/data/xenium-mouse-brain-SergioSalas/environment.yml"
    shell:
        "python {GIT_DIR}/data/xenium-mouse-brain-SergioSalas/prepare_data.py -o {output.dir}"
