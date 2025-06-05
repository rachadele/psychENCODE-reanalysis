import os
import argparse
import glob
import pandas as pd
from functools import reduce
import scanpy as sc
import pandas as pd
import numpy as np
import sys

def parse_arguments():
  parser = argparse.ArgumentParser(description="aggreate pseudobulk matrices by cell type from Gemma data")
  parser.add_argument("--author_degs", type=str, default = "/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/source_data/ASD_DEGcombined.csv")
  parser.add_argument("--contrast", type=str, default = "ASD")
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args

def main():
  args = parse_arguments()
  author_degs = pd.read_csv(args.author_degs, index_col=0)
  contrast = args.contrast
  # split by "cell_type" and save each subset to its own tsv
  # name new file after cell type and contrast
  for cell_type, ct_subset in author_degs.groupby("cell_type"):
    # check if "gene" column has duplicates
    if ct_subset["gene"].duplicated().any():
      # send warning to stderr
      sys.stderr.write(f"Warning: Duplicates found in {cell_type} for {contrast}.\n")
      ct_subset = ct_subset.drop_duplicates(subset=["gene"])
    ct_subset.to_csv(f"{contrast}_{cell_type}_degs.tsv", sep="\t", index=False)
  
if __name__ == "__main__":
  main()
