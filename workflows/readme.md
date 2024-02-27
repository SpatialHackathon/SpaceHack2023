## Workflow modification

* git_dir and data_dir/result_dir in every config file

## How to run snakemake

Running snakemake: download -> preprocessing -> methods -> metrics

* dry run: `snakemake -s <process>.smk -nf`

* actual run `snakemake -s <process>.smk --cores <n_of_cores> --use-conda --ri`
     * `ri`: in case you use keyboard interruption to quit the previous job. This will make sure snakemake rereun those incomplete job.

* Try not to kill snakemake when it's installing conda packages.