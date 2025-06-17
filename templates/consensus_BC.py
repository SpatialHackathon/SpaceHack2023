#!/usr/bin/env python

# Author_and_contribution: Jieran Sun & Mark Robinson; implmented method
# Author_and_contribution: Peiying Cai; created template
# Author_and_contribution: ENTER YOUR NAME AND CONTRIBUTION HERE

import argparse
import pandas as pd
from pathlib import Path

# TODO adjust description
description = "... to select base clusterings"

parser = argparse.ArgumentParser(description=description)

parser.add_argument(
    "-i", "--input_file", required=True,
    help="Input containing the aggregated labels."
)
parser.add_argument(
    "-o", "--output_file", required=True,
    help="Desired output file."
)
# TODO add additional arguments

args = parser.parse_args()

# Load input file
label_df = pd.read_csv(args.input_file, sep="\t", index_col=0)

## Your code goes here
# TODO
# output_df: DataFrame with number of clusters as columns, clustering label names as values
# Example:
#     7                          8
# 0   method1_default_7_label    method1_default_8_label
# 1   method2_default_7_label    method2_default_8_label
# 2   method3_default_7_label    method3_default_8_label

## Write output
output_path = Path(args.output_file)
output_path.parent.mkdir(parents=True, exist_ok=True)

# Save the results
output_df.to_csv(output_path, sep="\t", index=False)
