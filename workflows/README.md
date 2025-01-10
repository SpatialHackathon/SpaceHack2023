## Workflow modification

* git_dir and data_dir/result_dir in every config file

## How to run snakemake

Running snakemake: download -> preprocessing -> methods -> metrics

* dry run: `snakemake -s <process>.smk -nf`

* actual run `snakemake -s <process>.smk --cores <n_of_cores> --use-conda --ri`
     * `ri`: in case you use keyboard interruption to quit the previous job. This will make sure snakemake rereun those incomplete job.

* Try not to kill snakemake when it's installing conda packages.

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





