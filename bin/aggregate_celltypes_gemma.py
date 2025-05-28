import os
import argparse
import glob
import pandas as pd
from functools import reduce

def parse_arguments():
  parser = argparse.ArgumentParser(description="aggreate pseudobulk matrices by cell type from Gemma data")
  parser.add_argument("--pseudobulk_matrix_dir", type=str, default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/results/aggregated",
                      help="Path to the pseudobulk matrix files directory")
  parser.add_argument("--metadata_files", type=str, default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/gemma/metadata",
                      help="Path to the metadata files from Gemma")
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args

def main():
  args = parse_arguments()

  matrix_files = glob.glob(os.path.join(args.pseudobulk_matrix_dir, "*.tsv.gz"))
  pseudobulks = {}
  for file in matrix_files:
      df = pd.read_csv(file, sep="\t", comment="#", dtype=str, compression="gzip")
      key = os.path.basename(file).replace("_expmat.unfilt.aggregated.tsv.gz", "")
      pseudobulks[key] = df

  metadata_files = glob.glob(os.path.join(args.metadata_files, "*_sample_meta.tsv"))
  metadata = {}
  for file in metadata_files:
      key = os.path.basename(file).replace("_sample_meta.tsv", "")
      metadata[key] = pd.read_csv(file, sep="\t", dtype=str)

  matched_names = {}
  for gemma_key in pseudobulks:
      matched = [meta_key for meta_key in metadata if meta_key in gemma_key]
      matched_names[gemma_key] = matched[0] if len(matched) == 1 else None

  cell_type_pseudobulks = {}

  for gemma_key, meta_key in matched_names.items():
      if meta_key is None:
          continue
      mat = pseudobulks[gemma_key]
      #mat.columns.values[0:2] = ["ncbi_id", "gene_description"]
      mat["feature_name"] = mat["Sequence"].str.split(" ").str[0]
      mat.index = mat["feature_name"]
      new_mat = mat.drop(columns=["Probe", "Sequence", "feature_name"])

      cols = new_mat.columns
      mapping = pd.DataFrame({"original": cols})
      mapping["after_delim"] = mapping["original"].str.extract(r"___(.+)")
      mapping[["sample", "cell_type"]] = mapping["after_delim"].str.extract(r"([^\.]+)\.(.+)")

      for cell_type in mapping["cell_type"].dropna().unique():
          sample_ids = mapping.loc[mapping["cell_type"] == cell_type, "original"]
          mat_subset = new_mat[sample_ids]
          mat_subset.columns = mapping.loc[mapping["cell_type"] == cell_type, "sample"].values
          # Set cell type pseudobulks to dictionary for this cell type and fill with subset for this experiment
          cell_type_pseudobulks.setdefault(cell_type, {})[meta_key] = mat_subset

  for cell_type, dataset_mats in cell_type_pseudobulks.items():
    print(f"Processing cell type: {cell_type} with {len(dataset_mats)} datasets")
    combined = reduce(lambda left, right: pd.merge(left, right, left_index=True, right_index=True, how="outer"),
                      dataset_mats.values())
    outdir = cell_type.strip().replace("/", "_")
    os.makedirs(outdir, exist_ok=True)
    output_path = cell_type.strip().replace("/", "_") + "_pseudobulk_matrix.tsv.gz"
    combined.to_csv(os.path.join(outdir,output_path), sep="\t", index=True, na_rep=0, compression="gzip")
    
    # make some kind of summary plot for this pseudobulk matrix
    # PCA?
    

if __name__ == "__main__":
  main()
