# SpaceHack - contributing modules (data, methods, metrics)

Our workflow is set up to allow everyone to contribute "modules" in their preferred programming language (.. as long as that is either R or Python). A module can either be a dataset, a computational method, or an evaluation metric.
![image](https://github.com/SpatialHackathon/SpaceHack2023/assets/114547/7c002916-0a90-4fe7-8745-489313bc0192)

This repository contains some templates and examples of how to implement your module so that it interfaces seamlessly with other modules in the workflow. For example, if you want to implement a new method, you do not need to worry about input data or evaluation metrics as long as you follow the template for reading input and writing output - if you correctly adhere to the input and output guidelines, you should be able to interface with our default data modules and default evaluation metrics modules. The default modules are:
 - data: LIBD Visium DLPFC dataset (4 samples, each with 3 replicates)
 - methods: BayesSpace and SpaGCN
 - evaluation metrics: ARI and V

# How to contribute a module (method/dataset/metric)

Module contribution will be managed via GitHub. The steps to contribute a module are:
 1. Create or claim a **GitHub issue** from the [SpaceHack issue board.](https://github.com/SpatialHackathon/SpaceHack2023/issues) that describes the module you want to implement. There are currently 90 issues to claim, but if you come up with a new idea, please **create** a new issue, add the appropriate **tags**, and **assign** the task to yourself.
 2. Add **metadata** to our metadata [spreadsheet](https://docs.google.com/spreadsheets/d/1QCeAF4yQG4bhZSGPQwwVBj_XF7ADY_2mK5xivAIfHsc/edit). Please fill in as much as you can as metadata is helpful! If you feel the need, please also add new columns or add additional notes. The metadata should be added to the appropriate tabs:
    - [datasets](https://docs.google.com/spreadsheets/d/1QCeAF4yQG4bhZSGPQwwVBj_XF7ADY_2mK5xivAIfHsc/edit#gid=1453488771)
    - [computational methods](https://docs.google.com/spreadsheets/d/1QCeAF4yQG4bhZSGPQwwVBj_XF7ADY_2mK5xivAIfHsc/edit#gid=0)
    - [evaluation metrics](https://docs.google.com/spreadsheets/d/1QCeAF4yQG4bhZSGPQwwVBj_XF7ADY_2mK5xivAIfHsc/edit#gid=4776337)
    - [simulations and technical evaluation](https://docs.google.com/spreadsheets/d/1QCeAF4yQG4bhZSGPQwwVBj_XF7ADY_2mK5xivAIfHsc/edit#gid=640974611)
 3. Now you are ready to create a new git **[branch](https://learngitbranching.js.org/)**. Try to give your new branch an intuitive prefix such as `data_...`, `method_...` or `metric_..`. You can create a new branch in several ways: (i) [create a branch directly from the issue board](https://docs.github.com/en/issues/tracking-your-work-with-issues/creating-a-branch-for-an-issue) and then `git checkout` that branch, or (ii) via the command line:
```
# clone the template repository
git clone https://github.com/SpatialHackathon/SpaceHack2023.git
# create and switch to a new branch for your e.g. method "X"
git branch method_x_naveedishaque # try to make the branch name unique!
git checkout method_x_naveedishaque
# link the branch to the issue via the issue board: https://docs.github.com/en/issues/tracking-your-work-with-issues/linking-a-pull-request-to-an-issue
```

Modify the files, filenames, and code in `template/`, referring to the examples in the `data`, `method`, or `metric` subfolder. If your method requires a specific type or preprocessing, please reach out to the main developers!

**We support python and R, but we will explain everything in python.**

## Data

Data modules requires 3 files (see templates).

* `data.yml`: dependencies of the data module script following the format:
```
channels:
  - conda-forge
dependencies:
  - anndata=0.10.3
  - gitpython=3.1.40
```
* `data_optargs.json`: defining optional arguments for the workflow following the format:
```
{
    "min_cells" : 10,   # Minimum number of cell expressed required for a gene to pass filtering (int)
    "min_genes" : 20,   # Minimum number of genes expressed required for a cell to pass filtering (int)
    "min_counts": 30    # Minimum number of counts required for a cell to pass filtering (int)
}
```
* `data.py/.r`: data module script. 
   * Check the TODOs in the data.py or data.r template.
   * see further instruction below.

### Input Format

* `features_df`: DataFrame with rows representing features (e.g., genes) and columns representing additional metadata. Index: Feature ID or name.
* `observations_df`: DataFrame with rows representing observations (e.g., cells) and columns representing additional metadata. Index: Observation ID or barcode.
* `coordinates_df`: DataFrame with rows representing observations and columns (x, y, optionally z) for spatial coordinates. Index: Observation ID or barcode.
* `counts`: Matrix (2D array, e.g. .mtx file) with dimensions (#observations x #features). Matches the order of features_df and observations_df.

Optional Input Data

* `labels_df`: DataFrame with observation IDs as the index and a single column (label).
* `img`: Path to an optional image file (e.g., H&E stained image).

### Output Format

The output directory structure is organized as follows:

```
<out_dir>/
|___ sample_1/  (Sample name is user-defined)
|     |___ coordinates.tsv
|     |___ features.tsv
|     |___ observations.tsv
|     |___ counts.mtx  (Matrix Market format using `scipy.io.mmwrite`)
|     |___ labels.tsv  (Optional)
|     |___ H_E.(tiff/png/...)  (Optional)
|     |___ H_E.json  (Optional, required if H&E image is provided)
|
|___ sample_2/
|     |___ ...
|___ samples.tsv  (Metadata for all samples)
|___ experiment.json  (Contains metadata such as technology)
```

### Example usage of data scripts (Testing)

```
python data.py -o /path/to/output
```

### Add to workflow

* Add your data to the excute_config.yaml under `Dataset selected for excutation`.
* Add your data scripts to the path_config.yaml under `datasets`.


## Method

Method modules requires 3 files (see templates).

* `method.yml`: dependencies of the data module script following the format:
```
channels:
  - conda-forge
dependencies:
  - anndata=0.10.3
  - gitpython=3.1.40
```
* `method_optargs.json`: defining optional arguments for the workflow following the format:
```
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


* `method.py/.r`: method module script. 
   * Check the TODOs in the method.py or method.r template.
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


### Example usage of module scripts (Testing)

```
python method.py -c coordinates.tsv -f features.tsv -o observations.tsv \
    -m counts.mtx -d output_dir --n_clusters 5 --technology Visium --seed 42
```

### Add to workflow

* Add your method to the excute_config.yaml under `Methods selected for excutation`.
* Add your data scripts to the path_config.yaml under `methods`.


## Metric

Metric modules requires 3 files (see templates).

* `metric.yml`: dependencies of the data module script following the format:
```
channels:
  - conda-forge
dependencies:
  - anndata=0.10.3
  - gitpython=3.1.40
```
* `metric_optargs.json`: defining optional arguments for the workflow following the format:
```
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

* `metric.py/.r`: metric module script. 
   * Check the TODOs in the metric.py or metric.r template.
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

```
python metric.py -l labels.txt -g ground_truth.txt -o result.txt
```

### Add to workflow

* Add your method to the excute_config.yaml under `Metrics selected for excutation`.
* Add your data scripts to the path_config.yaml under `metrics`.

## Final steps

* Create a [pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request?tool=cli)
* Mark the code as `reads for review`, one of our developers will check it and merge your contributed module into the GitHub main branch!


### License

We have adopted the "MIT No Attribution" (MIT-0) License. It is currently attributed to the "SpaceHack organizers", but please also make sure to add your name to your contributions. More on MIT-0 [here](https://github.com/aws/mit-0)
