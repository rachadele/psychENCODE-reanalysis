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
  parser.add_argument("--h5ad_file", type=str, default = "/space/grp/rschwartz/rschwartz/get_gemma_data.nf/psychEncode_author_false_sample_split_false/homo_sapiens/Velmeshev_et_al.1.h5ad")
  parser.add_argument("--cell_type_column", type=str, default="cell_type")
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args

def get_avg_umi(adata, cell_type_column="cell_type"):
# compute average log umi per cell per sample 
  sc.pp.calculate_qc_metrics(adata, inplace=True)
  adata.obs["avg_UMI_ct"] = adata.obs.groupby(["sample_id", cell_type_column])["total_counts"].transform("mean")
  adata.obs["avg_UMI_sample"] = adata.obs.groupby("sample_id")["total_counts"].transform("mean")
  # return a vector mapping sample id to average umi
  umi_mapping = adata.obs[["sample_id", cell_type_column, "avg_UMI_sample", "avg_UMI_ct"]].value_counts().reset_index()
  return umi_mapping

# donors with <50 cells detected. 
# After filtering, cell types with <16 samples also removed

def filter_samples(pseudobulk, cell_counts, cell_type_column="cell_type", min_cells=50):
    subsets_to_keep = cell_counts[cell_counts["count"] >= min_cells][["sample_id", cell_type_column]]
    subsets_to_keep = subsets_to_keep.set_index(["sample_id", cell_type_column], drop=False)

    pseudobulk.obs = pseudobulk.obs.set_index(["sample_id", cell_type_column], drop=False)
    pseudobulk_filtered = pseudobulk[pseudobulk.obs.index.isin(subsets_to_keep.index)].copy()
    return pseudobulk_filtered

def aggregate_data(adata, cell_type_column="cell_type"):
    # generate pseudobulk matrix
    aggregated =  sc.get.aggregate(adata, by=["sample_id", cell_type_column], func=["sum", "count_nonzero", "mean"])
    return aggregated
  
  
def main():
  args = parse_arguments()
  h5ad_file = args.h5ad_file
  cell_type_column = args.cell_type_column
  adata = sc.read_h5ad(h5ad_file)
  cohort = os.path.basename(h5ad_file).replace(".h5ad", "")
  
  # drop genes with missing feature_name
  
  new_var = adata.var.dropna(subset=["feature_name"])
  adata = adata[:, new_var.index].copy()
  
  # get original metadata for sample id characteristics
  characteristics = [
      "sample_id",
      "sample_name",
      "1000G_ancestry",
      "Age_death",
      "Biological_Sex",
      "Cohort",
      "Disorder",
      "Genotype_data",
      "Individual_ID",
      "PMI",
      "RIN",
      "pH"
  ]

  # get sample id characteristics value counts
  sample_characteristics = adata.obs[characteristics].value_counts(dropna=False).reset_index()

  
  umi_mapping = get_avg_umi(adata, cell_type_column=cell_type_column)
  cell_counts = adata.obs[["sample_id", cell_type_column]].value_counts().reset_index()
  pseudobulk = aggregate_data(adata, cell_type_column=cell_type_column)
  pseudobulk = filter_samples(pseudobulk, cell_counts, cell_type_column=cell_type_column, min_cells=50)
  # add avg_UMI to obs
  pseudobulk.obs = pseudobulk.obs.merge(umi_mapping, on=["sample_id", cell_type_column], how="left")
  pseudobulk.obs = pseudobulk.obs.merge(sample_characteristics, on=["sample_id"], how="left") 
  pseudobulk.write_h5ad(os.path.join(f"{cohort}_pseudobulk.h5ad"))	


if __name__ == "__main__":
  main()
