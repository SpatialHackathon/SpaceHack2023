#!/usr/bin/env python

# Author_and_contribution: Jieran Sun; Implemented visualization

import argparse
import os
from pathlib import Path
import shutil


# TODO adjust description
parser = argparse.ArgumentParser(description="Merging all visualization pdfs into one big pdf")

parser.add_argument("-d", "--out_dir", help="Output directory.", required=True)
parser.add_argument("-p", "--pdfs", nargs="+", help="All input pdf datasets", required=True)

from PyPDF2 import PdfMerger

args = parser.parse_args()
out_dir = Path(args.out_dir)

# Create a merger object
merger = PdfMerger()

for pdf_file in args.pdfs:
    merger.append(pdf_file)
    pdf_file = Path(pdf_file)
    shutil.rmtree(pdf_file.parent, ignore_errors=True)

# write the output
out_dir.mkdir(parents=True, exist_ok=True)
merger.write(out_dir / "pp_report.pdf")