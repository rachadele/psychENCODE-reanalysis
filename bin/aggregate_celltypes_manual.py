import os
import argparse
import glob
import pandas as pd
from functools import reduce
import scanpy as sc
import pandas as pd
import numpy as np

def parse_arguments():
  parser = argparse.ArgumentParser(description="aggreate pseudobulk matrices by cell type from Gemma data")
  parser.add_argument("--h5ad_files", type=str, nargs="+", default = ["/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/results/aggregated/manual/DevBrain/DevBrain_pseudobulk.h5ad",
                                                                      "/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/results/aggregated/manual/MultiomeBrain/MultiomeBrain_pseudobulk.h5ad"])
  
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args

def main():
  args = parse_arguments()
  h5ad_files =args.h5ad_files
    
  ct_pseudobulks = {}
  
  for file in h5ad_files:
    adata = sc.read_h5ad(file)
    # print unique samples
    
    cohort = os.path.basename(file).replace("_pseudobulk.h5ad", "")
    celltypes = adata.obs["cell_type"].unique()
    for cell_type in celltypes:
      ct_subset = adata[adata.obs["cell_type"] == cell_type].copy()
      # extract the pseudobulk matrix
      matrix = ct_subset.layers["sum"].copy()
      # get hgnc
      var_names = ct_subset.var["feature_name"].tolist()
      new_pseudobulk = pd.DataFrame(matrix.T, index=var_names, columns = ct_subset.obs["sample_id"].tolist())
      # transpose matrix and make rownames var_names
      if cell_type not in ct_pseudobulks:
        ct_pseudobulks[cell_type] = {}
        # combine by row names
      ct_pseudobulks[cell_type][cohort] = new_pseudobulk
      
      # print the shape of the new pseudobulk matrix
      print(f"Processed {cohort} for cell type {cell_type}, shape: {new_pseudobulk.shape}")
      # print sample names 
      print(f"Sample names: {new_pseudobulk.columns.tolist()}")
     
      
  for cell_type, cohort in ct_pseudobulks.items():
    newname = cell_type.replace(" ", "").replace("/", "_")
    os.makedirs(newname, exist_ok=True)
    combined = reduce(lambda left, right: pd.merge(left, right, left_index=True, right_index=True, how="outer"),
                      ct_pseudobulks[cell_type].values())
    combined.index.name = "feature_name"
    print(combined.shape)
    combined.to_csv(os.path.join(newname,f"{newname}_pseudobulk_matrix.tsv.gz"), sep="\t", index=True, # fill NA with 0
                  na_rep="0", compression="gzip")
      
      # append columns to the pseudobulk matrix row wise
      

    
    
    
    
if __name__== "__main__":
    main()