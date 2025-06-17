# Consensus modules

## Implementing a new consensus method

To implement a new method follow the [Contribution guide](../CONTRIBUTING.md) and make sure you adopt all the necessary conventions specified in this document.

For examples, have a look at [select base clusterings](./02_BC_ranking/) and
[calculate consensus](./03_Consensus_kmode/).

## Layout and interface

Consensus module contains 3 steps:

1. Aggregate all label results into a single tsv file; see [here](./01_Results_Aggregation/).
2. Select the base-clusterings for consensus, either **automatically** or **manually**.
3. Run the consensus algorithm to obtain the final consensus labels.

### Base-clusterings selection

#### Manual Selection

For manual selection, create a TSV file specifying which base clusterings to use. Name the file **`BC_rankings.tsv`** and place it in results of the dataset: `/dataset/consensus/base_clusterings/Manual_selection/`.


The file should have number of clusters as column headers, and clustering label names as values. Each row corresponds to a method result.

```
7                          8
method1_default_7_label    method1_default_8_label  
method2_default_7_label    method2_default_8_label
method3_default_7_label    method3_default_8_label
```

#### Automatic Selection

The base-clusterings step requires 2 files (see templates). Replace `{consensus_BC}` in the file names with your method name, and place the files in the `consensus` folder.
* `{consensus_BC}.yaml`: a conda recipe defining the dependencies of the method module script following the format:

```yaml
channels:
  - r
  - conda-forge
dependencies:
  - r-base=4.4.2
  - r-optparse=1.7.5
```
    
* `{consensus_BC}.py/.r`: method module script. 
   * Check the TODOs in the `consensus_BC.py` or `consensus_BC.r` [template](../templates/).
   * The command line arguments can be modified. Further arguments can be passed using the `../workflows/excute_config.yaml` files.
   * see further instruction below.
   
##### Input Format

* `Aggregated Labels File (-i, --input_file)`: Path to a TSV file containing the aggregated labels for observations. Index: Observation ID or barcode. Columns: Clustering results named using the pattern `{method}_{config}_{n_clusters}_label`.

* Include any additional files required for selecting base clusterings.

##### Output Format

The script generates the following output file in the specified output directory and file name (`-o, --output_file`):

   * Contains selected clustering label names for the specified number of clusters.
   * Format: TSV with numbers of clusters as column headers and method configurations as rows (same format as manual annotation).
   
   
### Consensus calculation

Consensus calculation requires 2 files (see templates). Replace `{consensus}` in the file names with your method name, and place the files in the `consensus` folder.

* `{consensus}.yaml`: a conda recipe defining the dependencies of the method module script.

* `{consensus}.py/.r`: method module script. 
   * Check the TODOs in the `consensus.py` or `consensus.r` [template](../templates/).
   * The command line arguments are fixed and should not be modified. Further arguments can be passed using the `../workflows/excute_config.yaml` files.
   * see further instruction below.
   
#### Input Format

* `Aggregated Labels File (-i, --input_file)`: Path to a TSV file containing the aggregated labels for observations. Index: Observation ID or barcode. Columns: Clustering results named using the pattern `{method}_{config}_{n_clusters}_label`. Output from the first `Aggregation` step.
* `Base Clusterings File (-b, --base_clusterings)`: Path to a TSV file containing the chosen base clusterings for consensus calculation. Index: Method and config (e.g.,`scanpy_default_10_label`). Columns: Number of clusters. Output from the `Base-clusterings selection` step.

Parameters:

* `--n_clusters`: Number of clusters to return.
* `--n_bcs`: Number of base clustering results feed into the algorithm.
* `--seed`: Seed for random operations.

#### Output Format

The script generates the following output file in the specified output directory and file name (`-o, --output_file`):

   * Contains labels for observations.
   * Format: TSV with observation IDs as the index and a single label column.

## Example usage of module scripts (Testing)

```sh
Rscript consensus_BC.r -i combined_methods.tsv -o BC_rankings.tsv 
# TODO: add any parameters required
    
Rscript consensus.r -i combined_methods.tsv -b BC_rankings.tsv -o consensus.tsv \
    --n_clusters 8 --n_bcs 5 --seed 42
```

## Add to workflow

* Please request one of the organisers to add your algorithm scripts to the `06_select_base_clusterings.smk` and/or `07_consensus.smk` file.
* Add your consensus algorithm to the excute_config.yaml under `Consensus Clustering parameters`.
