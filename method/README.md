# Method modules

## Implementing a new method

To implement a new method follow the [Contribution guide](../contributing.md) and make sure you adopt all the necessary conventions specified in this document.

For examples have a look 
[here for a method in Python](./spaGCN/) or 
[here for a method in R](./BayesSpace/).

## Layout and interface

Method modules require 4 files (see templates). '{method}' in the file names should be
replaced by the name of your module and all files placed in a subfolder of the same name.

* `{method}.yml`: a conda recipe defining the dependencies of the method module script following the format:
    
```yaml
channels:
- conda-forge
dependencies:
- anndata=0.10.3
- gitpython=3.1.40
```
* `{method}_optargs.json`: defining optional arguments for the workflow following the format:
    
```json
{
    "matrix": "counts",
    "integrated_feature_selection": false,
    "image": true,
    "neighbors": false,
    "config_file": true
}
```

Entries are:

```
matrix: 
description: What input does the method take
type: string
enum:
    - counts
    - transform
    - dimensionality_reduction
    # - counts_or_transform

integrated_feature_selection:
description: Can the method use existing feature selections?
type: boolean

image:
description: Can the method use H&E images?
type: boolean

neighbors:
description: Can the method use existing neighbor definitions?
type: boolean

config_file:
description: Does the method take an additional config file?
type: boolean
```

* `config/config_default.json`: defining crucial parameters for the method that might need to be varied. 
The user can define multiple configs that can be tested with the workflow. 
The parameters in this file can be adjusted depending on the requirements of the methods. 
For example for conST:

```json
{
    "k": 10,
    "min_cells": 3,
    "use_img": false,
    "using_mask": false,
    "refinement": false,
    "source1": "https://github.com/ys-zong/conST/blob/main/conST_cluster.ipynb",
    "source2": "https://github.com/ys-zong/conST/blob/main/src/utils_func.py#L51"
}
```
    


* `{method}.py/.r`: method module script. 
   * Check the TODOs in the `method.py` or `method.r` [template]({{ repo_branch_url }}/templates/).
   * The command line arguments are fixed and should not be modified. Further arguments can be passed using the `config/config_{name}.json` files.
   * see further instruction below.

### Input Format

* `Coordinates File (-c, --coordinates)`: Path to a TSV file containing spatial coordinates. Index: Observation ID or barcode. Columns: x, y (and optionally z).
* `Features File (-f, --features)`: Path to a TSV file with rows representing features (e.g., genes). Index: Feature ID or name.
* `Observations File (-o, --observations)`: Path to a TSV file with rows representing observations (e.g., cells). Index: Observation ID or barcode.

Optional Files:

* `Matrix File (-m, --matrix)`: Path to a counts matrix in Matrix Market (MTX) format.
* `Neighbors File (-n, --neighbors)`: Path to a square matrix defining neighbors for each observation.
* `Dimensionality Reduction File (--dim_red)`: Path to reduced-dimensionality data (e.g., PCA) in TSV format.
* `Image File (--image)`: Path to an H&E stained image.
* `Config File (--config)`: Path to an optional JSON configuration file.

Parameters:

* `--n_clusters`: Number of clusters to return.
* `--technology`: Technology of the dataset (e.g., Visium, ST).
* `--seed`: Seed for random operations.

### Output Format

The script generates the following output files in the specified output directory (`-d, --out_dir`):

1. Domains File (`domains.tsv`):
   * Contains labels for observations.
   * Format: TSV with observation IDs as the index and a single label column.
2. Embedding File (`embedding.tsv`):
   * Optional output containing reduced-dimensionality representations.
   * Format: TSV with observation IDs as the index and n columns for embedding dimensions.

## Example usage of module scripts (Testing)

```sh
python method.py -c coordinates.tsv -f features.tsv -o observations.tsv \
    -m counts.mtx -d output_dir --n_clusters 5 --technology Visium --seed 42
```

## Add to workflow

* Add your method to the excute_config.yaml under `Methods selected for execution`.
* Add your method scripts to the path_config.yaml under `methods`.
