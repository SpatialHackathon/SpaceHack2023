## Workflow modification

* git_dir and data_dir/result_dir in every config file

## How to run snakemake

Running snakemake: download -> preprocessing -> methods -> metrics

* dry run: `snakemake -s <process>.smk -nf`

* actual run `snakemake -s <process>.smk --cores <n_of_cores> --use-conda --ri`
     * `ri`: in case you use keyboard interruption to quit the previous job. This will make sure snakemake rereun those incomplete job.

* Try not to kill snakemake when it's installing conda packages.

* If you're using a server or any HPC environment to run the workflow, it's recommended to use customized Snakemake profiles for job scheduling. You can find [HPC-specific Snakemake profile here](https://github.com/Snakemake-Profiles).

## Example usage

0. Use `excute_config_test.yaml` and `path_config_test.yaml` as your `excute_config.yaml` and `path_config_test.yaml` (just rename those files).
1. Download data

```
snakemake -s 01_download.smk --cores <n_of_cores> --use-conda --ri
```

2. Preprocess the data

```
snakemake -s 02_preprocessing.smk --cores <n_of_cores> --use-conda --ri
```

3. Execute method

```
snakemake -s 03_methods.smk --cores <n_of_cores> --use-conda --ri
```

4. Calculate metric

```
snakemake -s 04_metrics.smk --cores <n_of_cores> --use-conda --ri
```

5. Aggregate all the results

```
snakemake -s 05_aggregation.smk --cores <n_of_cores> --use-conda --ri
```

6. Create consensus

```
snakemake -s 06_consensus.smk --cores <n_of_cores> --use-conda --ri
```

## execute_config.yaml

You can use the file excute_config_test.yaml as a template for the execution of the workflow. The config follows the following structure: 

```
###### Universal parameters #######
# Directories, modify based on your own
GIT_DIR: path/for/github/repo
DATASET_DIR: path/to/datasets
SEED: 2023     # for the individual methods

###### Dataset selected for excutation #######
datasets_selected:
  - "list_datasets_you_want_to_analyse_and_are_placed_in_DATASET_DIR"

### Not used in this project
  - "list_datasets_you_do_not_want_to_use"

###### Methods selected for excutation #######
methods_selected:
## Native Implementation Done 
  - "list_methods_you_want_to_consider_for_the_consensus"

# If some datasets specify number of clusters. Add it here
n_clusters:
  visium_hd_cancer_colon: [5, 6, 7, 9, 11, 14]

###### Metrics selected for excutation #######
metrics_selected:
  - "list_implemented_metrices_that_you_want_to_consider_to_analyse_methods"

###### Base clustering selection parameters #######
# As used by scanpy (sc.pp.neighbors()).
selection_criteria:
  - "Cross_method_ARI"
  - "Smoothness_entropy"
  - "Manual_selection"
n_neighbors: 6

###### Consensus Clustering parameters #######
bc_numbers: [8]     # number of base clustering results
consensus_algorithms:
  - "lca"
  - "kmode"
  - "weighted"
# In case you need to re-define desired cluster number. Do it here. Otherwise n_clust value would be used
n_clust_consensus:
  abc_atlas_wmb_thalamus: [16, 19, 20, 21, 24, 28, 32]

# For weighted clustering
lambda: null

# For cross-method entropy
cross_method_entropy: true
```


## path_config.yaml

You can use the file path_config_test.yaml as a template for the execution of the workflow. The config follows the following structure: 

```
# The yaml file follows the following structure

# * categories (datasets/methods/metrics)
#   - {name}
#     - env: path/to/conda/env/.yaml
#     - script: path/to/script/.{py|r}
#     - env_additional: (optional)path/to/installation/script/.sh
#     - optargs: path/to/input/parameters/.json

# * config_files (for methods/metrics)
#   - {name} # MUST BE THE SAME AS THE METHOD/METRIC NAME
#     - {config_name}: path/to/config
#     - script: path to the excutation script
#     - env_additional: Only for certain methods, need installation shell script (.sh)
#     - optargs: optional arguments file (for input control/quality control)

# Notice for new addition:
# - name must be the same as the folder name!
# - All identation is 2 spaces!
# - When adding methods/metrics, remember to also add config_files if avaliable! 
# - Comment out configs that you don't want to run.
```

### You can generate the path_config also with generate_path_config.sh







