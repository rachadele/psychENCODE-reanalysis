#!/bin/python

import numpy as np
import pandas as pd
import os
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import matplotlib.ticker as ticker

def parse_args():
  parser = argparse.ArgumentParser(description="Explore pseudobulk data")
  parser.add_argument("--pseudobulk_dir", type=str,  help="Directory containing pseudobulk data files", default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/all_pseudobulk_sc_unfiltered")
  parser.add_argument("--output_dir", type=str, help="Directory to save the output plots", default="/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis/pseudobulk_exploration/")
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args


def plot_mean_variance(sc_df, output_dir, shortname):
  
  fig, axes = plt.subplots(2, 1, figsize=(10, 12))  # One figure, two subplots

  # Plot mean expression
  sns.histplot(sc_df.T.mean(), bins=100, ax=axes[0])
  axes[0].set_title("Mean Expression Distribution")
  axes[0].set_xlabel("Mean Expression")
  axes[0].set_ylabel("Frequency")

  # Plot variance
  sns.histplot(np.log10(sc_df.T.var()), bins=100, ax=axes[1])
  axes[1].axvline(np.log10(0.05), color='red', linestyle='--')
  axes[1].set_title("Variance Distribution")
  axes[1].set_xlabel("Log10 Variance")
  axes[1].set_ylabel("Frequency")

  # Optional overall title
  fig.suptitle(f"Expression Summary for {shortname}", fontsize=16)

  plt.tight_layout()  # Reserve space for suptitle
  plt.savefig(os.path.join(output_dir, f"{shortname}_mean_variance.png"), dpi=300)
  plt.show()
  plt.close(fig)  # Close the figure to free memory


def main():
  args = parse_args()
  pseudobulk_dir = args.pseudobulk_dir
  parent = args.output_dir
  os.makedirs(parent, exist_ok=True)
  
  for root, dirs, files in os.walk(args.pseudobulk_dir):
    for file in files:
      if file.endswith(".unfilt.data.txt.gz"):
        pseudobulk_file = os.path.join(root, file)
        
        shortname = os.path.basename(pseudobulk_file).split('_')[1]
        output_dir = os.path.join(args.output_dir, shortname)

        print(f"Processing {shortname}...")
        
        # Try to read the file
        try:
          sc_df = pd.read_table(pseudobulk_file, compression='gzip', comment='#')
        except Exception as e:
          print(f"Error reading {pseudobulk_file}: {e}")
          # write error to a file
          with open(os.path.join(parent, "error_files.txt"), "a") as error_file:
            error_file.write(f"{pseudobulk_file}\n")
          continue
        os.makedirs(os.path.join(parent,output_dir), exist_ok=True)
        sc_df = sc_df.set_index(['Probe', 'Sequence', 'GeneSymbol', 'GeneName', 'GemmaId', 'NCBIid'])

        plot_mean_variance(sc_df, output_dir, shortname)
       
 
if __name__ == "__main__":
  main()
 