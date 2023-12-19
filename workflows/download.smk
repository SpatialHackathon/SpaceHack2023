import os

configfile: "config.yaml"

GIT_DIR = os.getenv("GIT_DIR", "/home/ubuntu/workspace/SpaceHack2023")

def get_all_input(wildcards):
    all_folder = []
    for dataset in config["datasets"]:
        all_folder.append(config['results_dir'] + "/" + dataset)
    return all_folder

rule all:
    input:
        get_all_input


datasets = dict()
datasets["abc_atlas_wmb_thalamus"] = dict()
datasets["abc_atlas_wmb_thalamus"]["script"] = "python " + GIT_DIR + "/data/abc_atlas_wmb_thalamus/abc_atlas_wmb_thalamus.py"
datasets["abc_atlas_wmb_thalamus"]["env"] = GIT_DIR + "/data/abc_atlas_wmb_thalamus/abc_atlas_wmb_thalamus.yml"
datasets["libd_dlpfc"] = dict()
datasets["libd_dlpfc"]["script"] = "Rscript " + GIT_DIR + "/data/libd_dlpfc/libd_dlpfc.r"
datasets["libd_dlpfc"]["env"] = GIT_DIR + "/data/libd_dlpfc/libd_dlpfc.yml"
datasets["mouse_brain_sagittal_anterior"] = dict()
datasets["mouse_brain_sagittal_anterior"]["script"] = "python " + GIT_DIR + "/data/mouse_brain_sagittal_anterior/mouse_brain_sagittal_anterior.py"
datasets["mouse_brain_sagittal_anterior"]["env"] = GIT_DIR + "/data/mouse_brain_sagittal_anterior/mouse_brain_sagittal_anterior.yml"
datasets["mouse_brain_sagittal_posterior"] = dict()
datasets["mouse_brain_sagittal_posterior"]["script"] = "python " + GIT_DIR + "/data/mouse_brain_sagittal_posterior/mouse_brain_sagittal_posterior.py"
datasets["mouse_brain_sagittal_posterior"]["env"] = GIT_DIR + "/data/mouse_brain_sagittal_posterior/mouse_brain_sagittal_posterior.yml"
datasets["mouse_kidney_coronal"] = dict()
datasets["mouse_kidney_coronal"]["script"] = "python " + GIT_DIR + "/data/mouse_kidney_coronal/mouse_kidney_coronal.py"
datasets["mouse_kidney_coronal"]["env"] = GIT_DIR + "/data/mouse_kidney_coronal/mouse_kidney_coronal.yml"
datasets["osmfish_Ssp"] = dict()
datasets["osmfish_Ssp"]["script"] = "python " + GIT_DIR + "/data/osmfish_Ssp/osmfish_Ssp.py"
datasets["osmfish_Ssp"]["env"] = GIT_DIR + "/data/osmfish_Ssp/osmfish_Ssp.yml"
datasets["osmfish_Ssp"] = dict()
datasets["osmfish_Ssp"]["script"] = "python " + GIT_DIR + "/data/osmfish_Ssp/osmfish_Ssp.py"
datasets["osmfish_Ssp"]["env"] = GIT_DIR + "/data/osmfish_Ssp/osmfish_Ssp.yml"
datasets["SEA_AD_data"] = dict()
datasets["SEA_AD_data"]["script"] = "python " + GIT_DIR + "/data/SEA_AD_data/SEA_AD_data.py"
datasets["SEA_AD_data"]["env"] = GIT_DIR + "/data/SEA_AD_data/SEA_AD_data.yml"
datasets["sotip_simulation"] = dict()
datasets["sotip_simulation"]["script"] = "python " + GIT_DIR + "/data/sotip_simulation/sotip_simulation.py"
datasets["sotip_simulation"]["env"] = GIT_DIR + "/data/sotip_simulation/sotip.yml"
datasets["spatialDLPFC"] = dict()
datasets["spatialDLPFC"]["script"] = "Rscript " + GIT_DIR + "/data/spatialDLPFC/spatialDLPFC.r"
datasets["spatialDLPFC"]["env"] = GIT_DIR + "/data/spatialDLPFC/spatialDLPFC.yml"
datasets["STARmap_2018_mouse_cortex"] = dict()
datasets["STARmap_2018_mouse_cortex"]["script"] = "python " + GIT_DIR + "/data/STARmap-2018-mouse-cortex/STARmap-2018-mouse-cortex.py"
datasets["STARmap_2018_mouse_cortex"]["env"] = GIT_DIR + "/data/STARmap-2018-mouse-cortex/environment.yml"
datasets["visium_breast_cancer_SEDR"] = dict()
datasets["visium_breast_cancer_SEDR"]["script"] = "python " + GIT_DIR + "/data/visium_breast_cancer_SEDR/visium_breast_cancer_SEDR.py"
datasets["visium_breast_cancer_SEDR"]["env"] = GIT_DIR + "/data/visium_breast_cancer_SEDR/visium_breast_cancer_SEDR.yml"
datasets["visium_chicken_heart"] = dict()
datasets["visium_chicken_heart"]["script"] = "python " + GIT_DIR + "/data/visium_chicken_heart/chicken_heart.py"
datasets["visium_chicken_heart"]["env"] = GIT_DIR + "/data/visium_chicken_heart/chicken_heart.yml"
datasets["xenium_breast_cancer"] = dict()
datasets["xenium_breast_cancer"]["script"] = "python " + GIT_DIR + "/data/xenium-breast-cancer/xenium-breast-cancer.py"
datasets["xenium_breast_cancer"]["env"] = GIT_DIR + "/data/xenium-breast-cancer/xenium-breast-cancer.yml"
datasets["xenium_mouse_brain_SergioSalas"] = dict()
datasets["xenium_mouse_brain_SergioSalas"]["script"] = "python " + GIT_DIR + "/data/xenium-mouse-brain-SergioSalas/xenium-mouse-brain-SergioSalas.py"
datasets["xenium_mouse_brain_SergioSalas"]["env"] = GIT_DIR + "/data/xenium-mouse-brain-SergioSalas/environment.yml"


rule download:
    output:
        dir=directory(config['results_dir'] + "/{dataset}")
    conda:
        lambda wildcards: datasets[wildcards.dataset]["env"]
    params:
        script = lambda wildcards: datasets[wildcards.dataset]["script"]
    shell:
        "{params.script} -o {output.dir}"
