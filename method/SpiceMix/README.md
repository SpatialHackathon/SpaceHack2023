Here's an example of calling the `SpiceMix.py` script:
```{bash}
python SpiceMix.py -m ~/scratch/SpaceHack2/data/LIBD_DLPFC/Br8100_151673/counts.mtx -c ~/scratch/SpaceHack2/data/LIBD_DLPFC/Br8100_151673/coordinates.tsv -o ~/scratch/SpaceHack2/data/LIBD_DLPFC/Br8100_151673/observations.tsv -d ./output_test --n_clusters 7 --seed 0  --config config/config_1.json
```

The config file should look something like this:
```{json}

{
    "K": 15,
    "lambda_Sigma_x_inv": 1e-4,
    "device": "cuda:0",
    "dtype": "float64",
    "num_preiterations": 5,
    "num_iterations": 200
}
```
If you want preprocessing to be done within the script (such as log normalization, HVG selection, neighborhood graph construction), specify the `preprocess` parameter:
```{json}

{
    "K": 15,
    "lambda_Sigma_x_inv": 1e-4,
    "device": "cuda:0",
    "dtype": "float64",
    "num_preiterations": 5,
    "num_iterations": 200,
    "preprocess": {
        "hvgs": 3500
    }
}
```
