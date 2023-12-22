#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Jieran Sun; created script

import argparse

parser = argparse.ArgumentParser(description="Calculate Sørensen–Dice coefficient(SDC)")

parser.add_argument(
    "-l", "--labels", help="Labels from domain clustering.", required=True
)
parser.add_argument("-g", "--ground_truth", help="Groundtruth labels.", required=False)
parser.add_argument(
    "-e",
    "--embedding",
    help="Embedding of points in latent space. Potential usage for metrics without groundtruth.",
    required=False,
)
parser.add_argument(
    "-c",
    "--config",
    help="Optional config file used to pass additional parameters.",
    required=False,
)
parser.add_argument("-o", "--out_file", help="Output file.", required=True)
parser.add_argument(
    "--matched_labels",
    help="Flag indicating ground-truth and clustering labels have already been matched.",
    action='store_true',
)

args = parser.parse_args()

# Use these filepaths as input
label_file = args.labels

if args.ground_truth is not None:
    groundtruth_file = args.ground_truth
if args.config is not None:
    config_file = args.config

## Your code goes here

## Define SDC metric function, adopted from stlearn.metric.jaccard_score
def dice_coefficient(
    y_true,
    y_pred,
    *,
    labels=None,
    pos_label=1,
    average="binary",
    sample_weight=None,
    zero_division="warn",
):
    import numpy as np
    from sklearn.metrics._classification import _check_set_wise_labels, multilabel_confusion_matrix, _prf_divide
    
    labels = _check_set_wise_labels(y_true, y_pred, average, labels, pos_label)
    samplewise = average == "samples"
    MCM = multilabel_confusion_matrix(
        y_true,
        y_pred,
        sample_weight=sample_weight,
        labels=labels,
        samplewise=samplewise,
    )
    
    # calculating dice score
    numerator = 2*MCM[:, 1, 1]
    denominator = 2*MCM[:, 1, 1] + MCM[:, 0, 1] + MCM[:, 1, 0]

    if average == "micro":
        numerator = np.array([numerator.sum()])
        denominator = np.array([denominator.sum()])

    dice = _prf_divide(
        numerator,
        denominator,
        "dice",
        "true or predicted",
        average,
        ("dice",),
        zero_division=zero_division,
    )
    
    if average is None:
        return dice
    if average == "weighted":
        weights = MCM[:, 1, 0] + MCM[:, 1, 1]
        if not np.any(weights):
            # numerator is 0, and warning should have already been issued
            weights = None
    elif average == "samples" and sample_weight is not None:
        weights = sample_weight
    else:
        weights = None
    return np.average(dice, weights=weights)


if args.ground_truth is None:
    raise Exception("Groundtruth labels needed.")

import pandas as pd
from scipy.optimize import linear_sum_assignment

# Importing domain labels and ground truth label
domains = pd.read_table(label_file, index_col=0)["label"].astype("category").cat.codes
groundtruth = (
    pd.read_table(groundtruth_file, index_col=0)["label"].astype("category").cat.codes
)
# Filter common barcodes
common_index = domains.index.intersection(groundtruth.index)
groundtruth = groundtruth.loc[common_index]
domains = domains.loc[common_index]

# Matching the ground truth labels with the cluster labels using linear_sum_assignment
if not args.matched_labels:
    contingency_table = pd.crosstab(domains, groundtruth)
    row_ind, col_ind = linear_sum_assignment(contingency_table, maximize=True)
    domains = domains.map(dict(zip(row_ind, col_ind)))

# Using modified dice coefficient for calculation
metric = dice_coefficient(y_true=groundtruth, y_pred=domains, average='weighted')
domain_scores = dice_coefficient(y_true=groundtruth, y_pred=domains, average=None)

# Generating output dataframe
domains_df = pd.DataFrame({
    "cluster": ["all", *sorted(groundtruth.unique())],
    "dice_coefficient": [metric, *domain_scores]
})

## Write output
from pathlib import Path

Path(args.out_file).parent.mkdir(parents=True, exist_ok=True)

with open(args.out_file, "w") as file:
    file.write(domains_df.to_json(orient='records'))