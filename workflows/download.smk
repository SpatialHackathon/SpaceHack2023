import os

print("Run Download Workflow")

configfile: "config.yaml"

print(config)

GIT_DIR = os.getenv("GIT_DIR", "/home/ubuntu/workspace/SpaceHack2023")

print(GIT_DIR)

def get_all_input(wildcards):
    all_folder = []
    for dataset in config["datasets"]:
        all_folder.append(config['results_dir'] + "/" + dataset)
    return all_folder

rule all:
    input:
        get_all_input


datasets_info = dict()
datasets_info["abc_atlas_wmb_thalamus"] = dict()
datasets_info["abc_atlas_wmb_thalamus"]["script"] = "python " + GIT_DIR + "/data/abc_atlas_wmb_thalamus/abc_atlas_wmb_thalamus.py"
datasets_info["abc_atlas_wmb_thalamus"]["env"] = GIT_DIR + "/data/abc_atlas_wmb_thalamus/abc_atlas_wmb_thalamus.yml"
datasets_info["libd_dlpfc"] = dict()
datasets_info["libd_dlpfc"]["script"] = "Rscript " + GIT_DIR + "/data/libd_dlpfc/libd_dlpfc.r"
datasets_info["libd_dlpfc"]["env"] = GIT_DIR + "/data/libd_dlpfc/libd_dlpfc.yml"
datasets_info["mouse_brain_sagittal_anterior"] = dict()
datasets_info["mouse_brain_sagittal_anterior"]["script"] = "python " + GIT_DIR + "/data/mouse_brain_sagittal_anterior/mouse_brain_sagittal_anterior.py"
datasets_info["mouse_brain_sagittal_anterior"]["env"] = GIT_DIR + "/data/mouse_brain_sagittal_anterior/mouse_brain_sagittal_anterior.yml"
datasets_info["mouse_brain_sagittal_posterior"] = dict()
datasets_info["mouse_brain_sagittal_posterior"]["script"] = "python " + GIT_DIR + "/data/mouse_brain_sagittal_posterior/mouse_brain_sagittal_posterior.py"
datasets_info["mouse_brain_sagittal_posterior"]["env"] = GIT_DIR + "/data/mouse_brain_sagittal_posterior/mouse_brain_sagittal_posterior.yml"
datasets_info["mouse_kidney_coronal"] = dict()
datasets_info["mouse_kidney_coronal"]["script"] = "python " + GIT_DIR + "/data/mouse_kidney_coronal/mouse_kidney_coronal.py"
datasets_info["mouse_kidney_coronal"]["env"] = GIT_DIR + "/data/mouse_kidney_coronal/mouse_kidney_coronal.yml"
datasets_info["osmfish_Ssp"] = dict()
datasets_info["osmfish_Ssp"]["script"] = "python " + GIT_DIR + "/data/osmfish_Ssp/osmfish_Ssp.py"
datasets_info["osmfish_Ssp"]["env"] = GIT_DIR + "/data/osmfish_Ssp/osmfish_Ssp.yml"
datasets_info["osmfish_Ssp"] = dict()
datasets_info["osmfish_Ssp"]["script"] = "python " + GIT_DIR + "/data/osmfish_Ssp/osmfish_Ssp.py"
datasets_info["osmfish_Ssp"]["env"] = GIT_DIR + "/data/osmfish_Ssp/osmfish_Ssp.yml"
datasets_info["SEA_AD_data"] = dict()
datasets_info["SEA_AD_data"]["script"] = "python " + GIT_DIR + "/data/SEA_AD_data/SEA_AD_data.py"
datasets_info["SEA_AD_data"]["env"] = GIT_DIR + "/data/SEA_AD_data/SEA_AD_data.yml"
datasets_info["sotip_simulation"] = dict()
datasets_info["sotip_simulation"]["script"] = "python " + GIT_DIR + "/data/sotip_simulation/sotip_simulation.py"
datasets_info["sotip_simulation"]["env"] = GIT_DIR + "/data/sotip_simulation/sotip.yml"
datasets_info["spatialDLPFC"] = dict()
datasets_info["spatialDLPFC"]["script"] = "Rscript " + GIT_DIR + "/data/spatialDLPFC/spatialDLPFC.r"
datasets_info["spatialDLPFC"]["env"] = GIT_DIR + "/data/spatialDLPFC/spatialDLPFC.yml"
datasets_info["STARmap_2018_mouse_cortex"] = dict()
datasets_info["STARmap_2018_mouse_cortex"]["script"] = "python " + GIT_DIR + "/data/STARmap-2018-mouse-cortex/STARmap-2018-mouse-cortex.py"
datasets_info["STARmap_2018_mouse_cortex"]["env"] = GIT_DIR + "/data/STARmap-2018-mouse-cortex/environment.yml"
datasets_info["visium_breast_cancer_SEDR"] = dict()
datasets_info["visium_breast_cancer_SEDR"]["script"] = "python " + GIT_DIR + "/data/visium_breast_cancer_SEDR/visium_breast_cancer_SEDR.py"
datasets_info["visium_breast_cancer_SEDR"]["env"] = GIT_DIR + "/data/visium_breast_cancer_SEDR/visium_breast_cancer_SEDR.yml"
datasets_info["visium_chicken_heart"] = dict()
datasets_info["visium_chicken_heart"]["script"] = "python " + GIT_DIR + "/data/visium_chicken_heart/chicken_heart.py"
datasets_info["visium_chicken_heart"]["env"] = GIT_DIR + "/data/visium_chicken_heart/chicken_heart.yml"
datasets_info["xenium_breast_cancer"] = dict() # not working atm
datasets_info["xenium_breast_cancer"]["script"] = "python " + GIT_DIR + "/data/xenium-breast-cancer/xenium-breast-cancer.py"
datasets_info["xenium_breast_cancer"]["env"] = GIT_DIR + "/data/xenium-breast-cancer/xenium-breast-cancer.yml"
datasets_info["xenium_mouse_brain_SergioSalas"] = dict()
datasets_info["xenium_mouse_brain_SergioSalas"]["script"] = "python " + GIT_DIR + "/data/xenium-mouse-brain-SergioSalas/xenium-mouse-brain-SergioSalas.py"
datasets_info["xenium_mouse_brain_SergioSalas"]["env"] = GIT_DIR + "/data/xenium-mouse-brain-SergioSalas/environment.yml"


rule download:
    output:
        dir=directory(config['results_dir'] + "/{dataset}")
    conda:
        lambda wildcards: datasets_info[wildcards.dataset]["env"]
    params:
        script = lambda wildcards: datasets_info[wildcards.dataset]["script"]
    shell:
        "{params.script} -o {output.dir}"
