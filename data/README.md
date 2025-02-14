# Dataset modules

## Implementing a new dataset

To implement a new dataset follow the [Contribution guide](../contributing.md) and make sure you adopt all the necessary conventions specified in this document.

For an example have a look [here]({{ repo_branch_url }}/data/libd_dlpfc/).


## Layout and interface

Data modules require 3 files (see templates). '{data}' in the file names should be
replaced by the name of your module and all files placed in a subfolder of the same name.

* `{data}.yml`: a conda recipe defining the dependencies of the data module script following the format:
     
```yaml
channels:
- conda-forge
dependencies:
- anndata=0.10.3
- gitpython=3.1.40
```

* `{data}_optargs.json`: defining optional arguments for the workflow following the format:
     
```json
{
     "min_cells" : 10,   # Minimum number of cell expressed required for a gene to pass filtering (int)
     "min_genes" : 20,   # Minimum number of genes expressed required for a cell to pass filtering (int)
     "min_counts": 30    # Minimum number of counts required for a cell to pass filtering (int)
}
```

* `{data}.py/.r`: data module script. 
   * Check the TODOs in the `data.py` or `data.r` [template]({{ repo_branch_url }}/templates/).
   * The command line arguments are fixed and should not be modified.
   * see further instruction below.

### Input Format

* `features_df`: DataFrame with rows representing features (e.g., genes) and columns representing additional metadata. Index: Feature ID or name.
* `observations_df`: DataFrame with rows representing observations (e.g., cells) and columns representing additional metadata. Index: Observation ID or barcode.
* `coordinates_df`: DataFrame with rows representing observations and columns (x, y, optionally z) for spatial coordinates. Index: Observation ID or barcode.
* `counts`: Matrix (2D array, e.g. .mtx file) with dimensions (#observations x #features). Matches the order of features_df and observations_df.

Optional Input Data

* `labels_df`: DataFrame with observation IDs as the index and a single column (label).
* `img`: Path to an optional image file (e.g., H&E stained image).


### File structure 

```
dataset
├── sample_1
│   ├── coordinates.tsv
│   ├── counts.mtx
│   ├── features.tsv
│   ├── observations.tsv
│   ├── labels.tsv (optional)
│   └── H_E.{tiff,json} (optional)
├── sample_2
│   ├── …
│   └── …
├── experiment.json
└── samples.tsv
```

Keep the headers of all files the same - these are important for interfacing!

The index of the observations in the tsv-files, depending on the technology, could be a barcode, cell-ID, or similar.

#### `coordinates.tsv`
```
     x    y
AAAA 1234 9876
ATAC 1357 9753
CAAG 3579 7531
```

#### `counts.mtx`

This should be in MatrixMarket format.

#### `features.tsv`/`observations.tsv`

```
     row  col  selected
AAAA 1    1    true
ATAC 2    3    true
CAAG 5    2    false
```

The column `selected` is used for subsetting but is optional. `row` and `col` is needed in `observations.tsv` for bead-array based methods (Visium/ST).


#### `labels.tsv` 
Annotations of the ground truth domains

```
     label     label_confidence
AAAA Domain1   True
ATAC Domain1   True
CAAG Domain2   False
```
The column `label_confidence` is optional and used to indicate those cells 
and/or labels that are ground truth, if not all labels are high enough
confidence to be considered ground truth.


#### `image.tiff`
Images can be added in any format as appropriate (does not have to be tiff). If an image is available, please also add a json with relevant metadata (e.g. scale, but this might evolve during the hackathon)

#### `experiment.json`
Currently only technology (e.g. Visium, ST, MERSCOPE, MERFISH, Stereo-seq, Slide-seq, Xenium, STARmap, STARmap+, osmFISH, seqFISH) but more fields might be added.

#### `samples.tsv`
Sample directory and all relevant metadata, e.g. patient, replicate, slice, … and if applicable #clusters


## Example usage of data scripts (Testing)

```sh
python data.py -o /path/to/output
```

## Add to workflow 
<!-- TODO: update? -->

* Add your data to the excute_config.yaml under `Dataset selected for execution`.
* Add your data scripts to the path_config.yaml under `datasets`.
