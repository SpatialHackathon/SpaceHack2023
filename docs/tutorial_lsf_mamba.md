# SACCELERATOR Tutorial: Running the Pipeline with LSF (bsub) and Mamba

This tutorial covers:
1. Running the pipeline on two example datasets (`libd_dlpfc` and `xenium-breast-cancer`)
2. Managing Snakemake conda environments with proper names (avoiding hashed names)
3. Running the pipeline on new datasets (Xenium Prime FFPE Human Ovarian Cancer & Xenium Fresh Frozen Human Ovarian Cancer)

---

## Part 1: Running the Pipeline on Example Datasets

### 1.1 Prerequisites

- Access to an HPC cluster with LSF job scheduler (`bsub`)
- [Mamba](https://mamba.readthedocs.io/) installed (or Micromamba)
- Snakemake ≥ 7.x installed in a base environment
- Git

### 1.2 Clone the Repository

```bash
git clone https://github.com/SpatialHackathon/SACCELERATOR.git
cd SACCELERATOR
```

### 1.3 Set Up a Snakemake Environment

```bash
mamba create -n snakemake -c conda-forge -c bioconda snakemake mamba
mamba activate snakemake
```

### 1.4 Configure the Pipeline

You need to edit two config files in the `workflows/` directory.

#### `workflows/excute_config.yaml`

```yaml
###### Universal parameters #######
GIT_DIR: /path/to/your/SACCELERATOR    # Absolute path to your cloned repo
DATASET_DIR: /path/to/your/data         # Where downloaded data will be stored
SEED: 2023

###### Dataset selected for execution #######
datasets_selected:
  - "libd_dlpfc"
  - "xenium-breast-cancer"

###### Methods selected for execution #######
methods_selected:
  - "spaGCN"
  - "scanpy"
  - "seurat"

# Cluster numbers for datasets
n_clusters:
  libd_dlpfc: [7, 9, 11]
  xenium-breast-cancer: [5, 7, 9]

###### Metrics selected for execution #######
metrics_selected:
  - "ARI"
  - "V_measure"

###### Base clustering selection parameters #######
selection_criteria:
  - "Cross_method_ARI"
n_neighbors: 6

###### Consensus Clustering parameters #######
bc_numbers: [8]
consensus_algorithms:
  - "lca"
n_clust_consensus: {}

# For weighted clustering
lambda: null

# For cross-method entropy
cross_method_entropy: true
```

#### `workflows/path_config.yaml`

Use the provided `path_config.yaml` as-is, or regenerate it:

```bash
cd workflows
bash generate_path_config.sh
# Then rename the output:
mv generated_path_config.yml path_config.yaml
```

Make sure paths in `path_config.yaml` reference the correct relative locations for your selected datasets/methods/metrics.

### 1.5 Create a Snakemake Profile for LSF

Create a directory for your LSF profile:

```bash
mkdir -p ~/.config/snakemake/lsf
```

Create `~/.config/snakemake/lsf/config.yaml`:

```yaml
executor: cluster-generic
cluster-generic-submit-cmd: >
  bsub
  -J {rule}_{wildcards}
  -q normal
  -n {threads}
  -R "rusage[mem={resources.mem_mb}]"
  -M {resources.mem_mb}
  -W {resources.time_min}
  -o logs/{rule}_{wildcards}.out
  -e logs/{rule}_{wildcards}.err
cluster-generic-status-cmd: "bjobs -noheader -o 'stat' {}"
cluster-generic-cancel-cmd: "bkill {}"
jobs: 50
latency-wait: 60
use-conda: true
conda-frontend: mamba
rerun-incomplete: true
default-resources:
  mem_mb: 8000
  time_min: 120
  threads: 1
```

> **Note**: For Snakemake <8, use the older profile format. For Snakemake 8+, see the [snakemake-executor-plugin-lsf](https://github.com/snakemake/snakemake-executor-plugin-lsf) plugin.

Alternatively, for **Snakemake 8+** with the LSF executor plugin:

```bash
pip install snakemake-executor-plugin-lsf
```

Then your profile (`~/.config/snakemake/lsf/config.yaml`):

```yaml
executor: lsf
jobs: 50
latency-wait: 60
use-conda: true
conda-frontend: mamba
rerun-incomplete: true
default-resources:
  mem_mb: 8000
  runtime: 120
  lsf_queue: normal
```

Create a logs directory:

```bash
mkdir -p logs
```

### 1.6 Run the Pipeline Steps

Navigate to the workflows directory:

```bash
cd /path/to/your/SACCELERATOR/workflows
```

Run each step sequentially. Use the `--profile lsf` flag to submit jobs via bsub:

```bash
# Step 1: Download data
snakemake -s 01_download.smk --profile lsf

# Step 2: Preprocessing (QC, normalization, etc.)
snakemake -s 02_preprocessing.smk --profile lsf

# Step 3: Run clustering methods
snakemake -s 03_methods.smk --profile lsf

# Step 4: Calculate metrics
snakemake -s 04_metrics.smk --profile lsf

# Step 5: Aggregate results
snakemake -s 05_aggregation.smk --profile lsf

# Step 6: Select base clusterings
snakemake -s 06_select_base_clusterings.smk --profile lsf

# Step 7: Consensus clustering
snakemake -s 07_consensus.smk --profile lsf
```

#### Dry Run First

Always do a dry run before submitting:

```bash
snakemake -s 01_download.smk --profile lsf -n
```

---

## Part 2: Conda Environment Naming — Avoiding Hashed Environment Names

### The Problem

By default, Snakemake creates conda environments using a hash of the YAML file content. This results in environment directories like:

```
.snakemake/conda/a1b2c3d4e5f6/
```

These are hard to inspect, debug, or reuse across projects.

### Solution: Use Named Environments with `--conda-prefix` and Pre-created Envs

There are several approaches to get properly named, reusable environments:

#### Approach A: Shared `--conda-prefix` Directory (Recommended for Reuse)

Set a shared prefix so all Snakemake runs (across projects) share the same conda env cache:

```bash
# Set a shared location for conda envs
export SNAKEMAKE_CONDA_PREFIX=/path/to/shared/conda_envs

snakemake -s 03_methods.smk --profile lsf --conda-prefix $SNAKEMAKE_CONDA_PREFIX
```

This ensures environments are reused across runs. The hashed names persist, but they are cached and shared.

You can also add this to your profile config:

```yaml
# In ~/.config/snakemake/lsf/config.yaml
conda-prefix: /path/to/shared/conda_envs
```

#### Approach B: Pre-create Named Environments and Use `--conda-prefix` with Symlinks

Create a wrapper script that pre-creates environments with meaningful names and symlinks them:

```bash
#!/bin/bash
# File: setup_named_envs.sh
# Pre-create named conda environments from SACCELERATOR YAML files

GIT_DIR="/path/to/your/SACCELERATOR"
ENV_DIR="/path/to/shared/conda_envs"

mkdir -p "$ENV_DIR"

# Dataset environments
for yml in "$GIT_DIR"/data/*//*.yml "$GIT_DIR"/data/*/*.yaml; do
    [ -f "$yml" ] || continue
    name=$(basename "$(dirname "$yml")")
    echo "Creating env: data_${name}"
    mamba env create -f "$yml" -n "saccelerator_data_${name}" --yes 2>/dev/null || \
    mamba env update -f "$yml" -n "saccelerator_data_${name}" --prune
done

# Method environments
for yml in "$GIT_DIR"/method/*/*.yml "$GIT_DIR"/method/*/*.yaml; do
    [ -f "$yml" ] || continue
    name=$(basename "$(dirname "$yml")")
    echo "Creating env: method_${name}"
    mamba env create -f "$yml" -n "saccelerator_method_${name}" --yes 2>/dev/null || \
    mamba env update -f "$yml" -n "saccelerator_method_${name}" --prune
done

# Metric environments
for yml in "$GIT_DIR"/metric/*/*.yml "$GIT_DIR"/metric/*/*.yaml; do
    [ -f "$yml" ] || continue
    name=$(basename "$(dirname "$yml")")
    echo "Creating env: metric_${name}"
    mamba env create -f "$yml" -n "saccelerator_metric_${name}" --yes 2>/dev/null || \
    mamba env update -f "$yml" -n "saccelerator_metric_${name}" --prune
done

# Preprocessing environments
for yml in "$GIT_DIR"/preprocessing/*/*.yml "$GIT_DIR"/preprocessing/*/*/*.yml; do
    [ -f "$yml" ] || continue
    name=$(basename "$yml" .yml)
    echo "Creating env: preproc_${name}"
    mamba env create -f "$yml" -n "saccelerator_preproc_${name}" --yes 2>/dev/null || \
    mamba env update -f "$yml" -n "saccelerator_preproc_${name}" --prune
done

echo "All environments created. List them with: mamba env list"
```

#### Approach C: Use `envmodules` Directive (Alternative for HPC)

If your cluster has environment modules, you can modify the Snakemake rules to use `envmodules:` instead of `conda:`. However, this requires modifying the `.smk` files, which is not recommended for this workflow.

#### Approach D: Create Environments Before Running, Then Use `--conda-create-envs-only`

First, let Snakemake create the environments without running the pipeline:

```bash
snakemake -s 03_methods.smk --profile lsf --conda-create-envs-only
```

This pre-creates all required environments. Once done, subsequent runs will reuse them automatically. Combined with `--conda-prefix`, this gives you persistent, reusable environments.

#### How to Map Hashed Names to Tools

To find out which hash corresponds to which tool:

```bash
# List all conda envs created by Snakemake
ls .snakemake/conda/

# Check what's in each env
for dir in .snakemake/conda/*/; do
    if [ -f "${dir}environment.yaml" ]; then
        echo "=== $dir ==="
        cat "${dir}environment.yaml"
        echo ""
    fi
done
```

Or use Snakemake's built-in listing:

```bash
snakemake -s 03_methods.smk --list-conda-envs
```

---

## Part 3: Running on New Datasets — Xenium Ovarian Cancer

This section shows how to add and run the pipeline on two new 10x Genomics datasets:
- [Xenium Prime FFPE Human Ovarian Cancer](https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-ovarian-cancer)
- [Xenium Comparison Fresh Frozen Human Ovarian Cancer](https://www.10xgenomics.com/datasets/xenium-comparison-fresh-frozen-human-ovarian-cancer)

### 3.1 Create a Data Module for the New Datasets

Each dataset in SACCELERATOR needs:
1. A download/processing script (Python or R)
2. A conda environment YAML
3. An optargs JSON file
4. An entry in `path_config.yaml`
5. An entry in `excute_config.yaml`

#### Directory Structure

```
data/
  xenium-ovarian-cancer-ffpe/
    xenium-ovarian-cancer-ffpe.py
    xenium-ovarian-cancer-ffpe.yml
    xenium-ovarian-cancer-ffpe_optargs.json
  xenium-ovarian-cancer-ff/
    xenium-ovarian-cancer-ff.py
    xenium-ovarian-cancer-ff.yml
    xenium-ovarian-cancer-ff_optargs.json
```

#### 3.2 Create the Conda Environment YAML

`data/xenium-ovarian-cancer-ffpe/xenium-ovarian-cancer-ffpe.yml`:

```yaml
channels:
  - conda-forge
dependencies:
  - python=3.10
  - anndata>=0.10
  - pandas>=2.0
  - scipy>=1.10
  - requests>=2.28
  - pip
  - pip:
    - spatialdata-io>=0.1
    - pypdl>=0.5
```

Use the same for the FF dataset (or adjust as needed):

`data/xenium-ovarian-cancer-ff/xenium-ovarian-cancer-ff.yml`:

```yaml
channels:
  - conda-forge
dependencies:
  - python=3.10
  - anndata>=0.10
  - pandas>=2.0
  - scipy>=1.10
  - requests>=2.28
  - pip
  - pip:
    - spatialdata-io>=0.1
    - pypdl>=0.5
```

#### 3.3 Create the Optargs JSON

`data/xenium-ovarian-cancer-ffpe/xenium-ovarian-cancer-ffpe_optargs.json`:

```json
{
    "min_cells": 1,
    "min_genes": 1,
    "min_counts": 1
}
```

Same for the FF dataset.

#### 3.4 Create the Download/Processing Script

The script must:
- Accept `-o/--out_dir` as argument
- Download raw data from 10x Genomics
- Output the standardized format:
  - `{out_dir}/experiment.json` — technology metadata
  - `{out_dir}/samples.tsv` — sample table
  - `{out_dir}/{sample_name}/counts.mtx` — cell×gene count matrix (Market Matrix format)
  - `{out_dir}/{sample_name}/features.tsv` — gene metadata
  - `{out_dir}/{sample_name}/observations.tsv` — cell metadata
  - `{out_dir}/{sample_name}/coordinates.tsv` — spatial coordinates (columns: x, y)
  - `{out_dir}/{sample_name}/labels.tsv` — ground truth labels (if available)

**`data/xenium-ovarian-cancer-ffpe/xenium-ovarian-cancer-ffpe.py`**:

```python
#!/usr/bin/env python

"""
Download and process Xenium Prime FFPE Human Ovarian Cancer dataset.
Source: https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-ovarian-cancer
"""

import argparse
import json
import os
import shutil
import tempfile

import pandas as pd
import scipy.io
from pypdl import Downloader
from spatialdata_io import xenium


# Update these URLs from the 10x Genomics dataset page.
# Go to the dataset page and find the direct download links for the output bundle.
LINKS = {
    # "https://cf.10xgenomics.com/samples/xenium/.../_outs.zip": "md5_checksum",
    # Add the actual download URL once you retrieve it from the 10x dataset page
}


def download_links(links, temp_dir):
    """Download all files from the links dict."""
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 6.1; Win64; x64; rv:47.0) Gecko/20100101 Firefox/47.0"
    }
    dl = Downloader(headers=headers)
    for link, checksum in links.items():
        print(f"Downloading {link}")
        file = dl.start(
            url=link,
            file_path=temp_dir,
            segments=10,
            display=True,
            multithread=True,
            block=True,
            retries=3,
        )
        if checksum and not file.validate_hash(checksum, "md5"):
            raise ValueError(f"File {file} is corrupted")


def process_xenium_output(xenium_path, out_path, sample_name):
    """Convert Xenium output to SACCELERATOR format."""
    print(f"Processing {xenium_path} -> {out_path}/{sample_name}")

    # Read using spatialdata-io
    sdata = xenium(
        xenium_path,
        cells_boundaries=False,
        nucleus_boundaries=False,
        cells_as_circles=False,
        cells_labels=False,
        nucleus_labels=False,
        transcripts=False,
        morphology_mip=False,
        morphology_focus=False,
    )
    adata = sdata["table"]

    complete_path = os.path.join(out_path, sample_name)
    os.makedirs(complete_path, exist_ok=True)

    # Observations
    obs = adata.obs.copy()
    obs["selected"] = "true"
    obs.to_csv(f"{complete_path}/observations.tsv", sep="\t", index_label="")

    # Features
    features = adata.var.copy()
    features["selected"] = "true"
    features.to_csv(f"{complete_path}/features.tsv", sep="\t", index_label="")

    # Coordinates
    coords = pd.DataFrame(adata.obsm["spatial"], columns=["x", "y"])
    coords.index = adata.obs.index
    coords.to_csv(f"{complete_path}/coordinates.tsv", sep="\t", index_label="")

    # Count matrix (cells x genes in MMwrite)
    scipy.io.mmwrite(f"{complete_path}/counts.mtx", adata.X)

    return adata.shape[0]  # return number of cells


def write_experiment_json(out_path):
    experiment = {
        "technology": "Xenium",
        "species": "human",
        "tissue": "ovary",
        "is_3D": False,
    }
    with open(os.path.join(out_path, "experiment.json"), "w") as f:
        json.dump(experiment, f, indent=2)


def write_samples_tsv(out_path, samples_info):
    df = pd.DataFrame(samples_info)
    df.to_csv(f"{out_path}/samples.tsv", sep="\t", index_label=False)


def main():
    parser = argparse.ArgumentParser(
        description="Download Xenium Prime FFPE Human Ovarian Cancer dataset."
    )
    parser.add_argument(
        "-o", "--out_dir", help="Output directory to write files to.", required=True
    )
    args = parser.parse_args()

    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    with tempfile.TemporaryDirectory() as temp_dir:
        # Download
        download_links(LINKS, temp_dir)

        # Unzip and process
        for file in os.listdir(temp_dir):
            if file.endswith(".zip"):
                sample_name = "FFPE_ovarian"
                sample_path = os.path.join(temp_dir, file.replace(".zip", ""))
                shutil.unpack_archive(os.path.join(temp_dir, file), sample_path)
                process_xenium_output(sample_path, out_dir, sample_name)

    # Write metadata
    write_experiment_json(out_dir)
    write_samples_tsv(out_dir, {
        "patient": ["ovarian_ffpe"],
        "sample": ["FFPE_ovarian"],
        "position": [0],
        "replicate": [0],
        "directory": ["FFPE_ovarian"],
        "n_clusters": [7],  # Adjust based on expected domains
    })


if __name__ == "__main__":
    main()
```

> **Important**: You need to get the actual download URLs from the 10x Genomics dataset pages. Visit:
> - https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-ovarian-cancer
> - https://www.10xgenomics.com/datasets/xenium-comparison-fresh-frozen-human-ovarian-cancer
>
> Look for the "Output Files" section and copy the direct download links for the `*_outs.zip` bundle. Update the `LINKS` dictionary accordingly.

Create a similar script for `xenium-ovarian-cancer-ff` adjusting the URLs and sample names.

#### 3.5 Make Scripts Executable

```bash
chmod +x data/xenium-ovarian-cancer-ffpe/xenium-ovarian-cancer-ffpe.py
chmod +x data/xenium-ovarian-cancer-ff/xenium-ovarian-cancer-ff.py
```

#### 3.6 Register the Datasets in Config Files

**Add to `workflows/path_config.yaml`** under the `datasets:` section:

```yaml
  xenium-ovarian-cancer-ffpe:
    env: data/xenium-ovarian-cancer-ffpe/xenium-ovarian-cancer-ffpe.yml
    script: data/xenium-ovarian-cancer-ffpe/xenium-ovarian-cancer-ffpe.py
    optargs: data/xenium-ovarian-cancer-ffpe/xenium-ovarian-cancer-ffpe_optargs.json
  xenium-ovarian-cancer-ff:
    env: data/xenium-ovarian-cancer-ff/xenium-ovarian-cancer-ff.yml
    script: data/xenium-ovarian-cancer-ff/xenium-ovarian-cancer-ff.py
    optargs: data/xenium-ovarian-cancer-ff/xenium-ovarian-cancer-ff_optargs.json
```

**Update `workflows/excute_config.yaml`**:

```yaml
datasets_selected:
  - "xenium-ovarian-cancer-ffpe"
  - "xenium-ovarian-cancer-ff"

methods_selected:
  - "spaGCN"
  - "scanpy"
  - "CellCharter"
  - "BANKSY"

n_clusters:
  xenium-ovarian-cancer-ffpe: [5, 7, 9, 11]
  xenium-ovarian-cancer-ff: [5, 7, 9, 11]
```

### 3.7 Run the Pipeline

```bash
cd /path/to/your/SACCELERATOR/workflows

# Step 1: Download and format data
snakemake -s 01_download.smk --profile lsf

# Step 2: Preprocessing
snakemake -s 02_preprocessing.smk --profile lsf

# Step 3: Run methods
snakemake -s 03_methods.smk --profile lsf

# Step 4: Metrics
snakemake -s 04_metrics.smk --profile lsf

# Step 5: Aggregate
snakemake -s 05_aggregation.smk --profile lsf

# Step 6: Select base clusterings
snakemake -s 06_select_base_clusterings.smk --profile lsf

# Step 7: Consensus
snakemake -s 07_consensus.smk --profile lsf
```

### 3.8 Getting 10x Genomics Download URLs

The 10x Genomics dataset pages require agreeing to their terms. To get the direct URLs:

1. Go to https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-ovarian-cancer
2. Click "Download" and accept terms
3. Copy the link for "Feature / cell matrix (HDF5)" or the full "Output file" zip
4. For Xenium data, the key files you need are in the `*_outs.zip`:
   - `cells.csv.gz` or `cells.parquet` (cell positions and metadata)
   - `cell_feature_matrix/` (count matrix in MEX format)
   - `gene_panel.json` (panel info)

Alternatively, use `curl` with the direct S3/CloudFront URLs:

```bash
# Example pattern (URLs change per dataset):
# https://cf.10xgenomics.com/samples/xenium/3.0.0/Xenium_Prime_FFPE_Human_Ovarian_Cancer/Xenium_Prime_FFPE_Human_Ovarian_Cancer_outs.zip
```

---

## Troubleshooting

### Common Issues

1. **"Directory does not exist" errors**: Make sure `DATASET_DIR` exists and the download step completed successfully.

2. **Conda environment creation fails**: Use `--conda-create-envs-only` to create envs first without running rules.

3. **LSF job memory issues**: Increase `mem_mb` in your profile or add per-rule resources in the Snakefiles.

4. **Method incompatible with technology**: Some methods (e.g., BayesSpace, GraphST) only support certain technologies. Check the `optargs.json` → `"technology"` field for each method.

5. **Snakemake locked**: If a previous run failed mid-execution:
   ```bash
   snakemake -s <step>.smk --unlock
   ```

### Useful Commands

```bash
# Check what would run (dry run)
snakemake -s 03_methods.smk --profile lsf -n

# Force re-run of incomplete jobs
snakemake -s 03_methods.smk --profile lsf --rerun-incomplete

# Run on a single specific rule/target
snakemake -s 03_methods.smk --profile lsf /path/to/data/xenium-ovarian-cancer-ffpe/FFPE_ovarian/spaGCN/config_default/cluster_7/domains.tsv

# Show the DAG (dependency graph)
snakemake -s 03_methods.smk --dag | dot -Tpng > dag.png
```

---

## Summary

| Step | Command |
|------|---------|
| Download data | `snakemake -s 01_download.smk --profile lsf` |
| Preprocess | `snakemake -s 02_preprocessing.smk --profile lsf` |
| Run methods | `snakemake -s 03_methods.smk --profile lsf` |
| Compute metrics | `snakemake -s 04_metrics.smk --profile lsf` |
| Aggregate | `snakemake -s 05_aggregation.smk --profile lsf` |
| Select base clusterings | `snakemake -s 06_select_base_clusterings.smk --profile lsf` |
| Consensus | `snakemake -s 07_consensus.smk --profile lsf` |

Key flags:
- `--profile lsf` — use your LSF profile for job submission
- `--conda-frontend mamba` — use Mamba for faster env creation
- `--conda-prefix /shared/path` — share envs across projects
- `--conda-create-envs-only` — pre-create all environments without running
- `-n` — dry run
- `--rerun-incomplete` or `--ri` — re-run interrupted jobs
