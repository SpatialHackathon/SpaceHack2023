# Metric modules

## Implementing a new dataset module

To implement a new metric follow the [Contribution guide](../contributing.md) and make sure you adopt all the necessary conventions specified in this document.

For examples have a look 
[here for a method in Python]({{ repo_branch_url }}/metric/ARI/) or 
[here for a method in R]({{ repo_branch_url }}/metric/LISI/).

## Metric module layout and interface

Metric modules require 3 files (see templates). '{metric}' in the file names should be
replaced by the name of your module and all files placed in a subfolder of the same name.

* `{metric}.yml`: dependencies of the metric module script following the format:
```yaml
channels:
- conda-forge
dependencies:
- anndata=0.10.3
- gitpython=3.1.40
```

* `{metric}_optargs.json`: defining optional arguments for the workflow following the format:
```json
{
    "groundtruth": true,   # Does the metric need groundtruth labels? (boolean)
    "embedding": false,    # Does the metric need embeddings? (boolean)
    "config_file": true    # Does the metric take an additional config file? (boolean)
}
```

It is optional to add when your metric requires this:

```
physical_coordinate:
description: Does the metric take physcial coordination of the sample?
type: boolean
```

* `{metric}.py/.r`: metric module script. 
   * Check the TODOs in the `metric.py` or `metric.r` [template]({{ repo_branch_url }}/templates/).
   * The command line arguments are fixed and should not be modified.
   * see further instruction below.


### Input Format

* `Labels File (-l, --labels)`: Path to a file containing cluster labels. Format: Text file where each row corresponds to a label for a specific observation.

Optional Files:

* `Ground Truth File (-g, --ground_truth)`: Path to a file containing ground truth labels. Use this for metrics requiring true labels for comparison.
* `Embedding File (-e, --embedding)`: Path to a file containing latent space embeddings. Useful for metrics that do not rely on ground truth labels.
* `Config File (-c, --config)`: Path to an optional JSON file with additional parameters for metric calculation.

### Output Format

The script writes the calculated metric to the specified output file (`-o, --out_file`) in scientific notation with five decimal places.

### Example usage of module scripts (Testing)

```sh
python metric.py -l labels.txt -g ground_truth.txt -o result.txt
```

### Add to workflow

* Add your metric to the excute_config.yaml under `Metrics selected for execution`.
* Add your metric scripts to the path_config.yaml under `metrics`.
