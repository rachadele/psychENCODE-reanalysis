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
  parser.add_argument("--h5ad_files", type=str, nargs="+", default = ["/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/results/experiment_pseudobulks/manual/DevBrain/DevBrain_pseudobulk.h5ad",
                                                                      "/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/results/experiment_pseudobulks/manual/MultiomeBrain/MultiomeBrain_pseudobulk.h5ad"])
  
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args
  
def filter_samples(adata, min_cells=50):
    cell_counts_per_samples = adata.obs["sample_id"].value_counts().reset_index()
    samples_to_keep = cell_counts_per_samples[cell_counts_per_samples >= min_cells].index
    adata_filtered = adata[adata.obs["sample_id"].isin(samples_to_keep)].copy()
    return adata_filtered

    
def main():
  args = parse_arguments()
  h5ad_files =args.h5ad_files
    
  ct_pseudobulks = {}
  meta_pseudobulks = {}
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
      meta = ct_subset.obs.copy()
      new_pseudobulk = pd.DataFrame(matrix.T, index=var_names, columns = ct_subset.obs["sample_id"].tolist())
      # transpose matrix and make rownames var_names
      if cell_type not in ct_pseudobulks:
        ct_pseudobulks[cell_type] = {}
      if cell_type not in meta_pseudobulks:
        meta_pseudobulks[cell_type] = {}
        # combine by row names
      ct_pseudobulks[cell_type][cohort] = new_pseudobulk
      meta_pseudobulks[cell_type][cohort] = meta
      
      # print the shape of the new pseudobulk matrix
      print(f"Processed {cohort} for cell type {cell_type}, shape: {new_pseudobulk.shape}")
      # print sample names 
      #print(f"Sample names: {new_pseudobulk.columns.tolist()}")
     
      
  for cell_type, cohorts in ct_pseudobulks.items():
    print(f"Processing cell type: {cell_type}")
    print(f"Number of cohorts: {len(cohorts)}")
    newname = cell_type.replace(" ", "").replace("/", "_")
    os.makedirs(newname, exist_ok=True)
    combined = reduce(lambda left, right: pd.merge(left, right, left_index=True, right_index=True, how="outer"),
                      ct_pseudobulks[cell_type].values())
    combined.index.name = "feature_name"
    print(combined.shape)
    # if less than 16 samples, skip saving
    if combined.shape[1] < 16:
        print(f"Skipping {newname} as it has less than 16 samples.")
        continue
    combined.to_csv(os.path.join(newname,f"{newname}_pseudobulk_matrix.tsv.gz"), sep="\t", index=True, # fill NA with 0
                  na_rep="0", compression="gzip")
     
     # concatenate metadata row wise
    meta_combined = pd.concat(meta_pseudobulks[cell_type].values(), axis=0)
    meta_combined.index.name = "index" 
    meta_combined.to_csv(os.path.join(newname, f"{newname}_pseudobulk_metadata.tsv"), sep="\t", index=False)
    # print metadata shape
    print(f"Metadata shape: {meta_combined.shape}")

    
    
    
    
if __name__== "__main__":
    main()