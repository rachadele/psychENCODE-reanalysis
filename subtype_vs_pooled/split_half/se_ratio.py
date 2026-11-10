"""Part 1: precision of existing limma-voom fits, author subtype vs pooled (merged) author type.

For each pooled group (IT, Lamp5, Sst), each of its author subtypes and each of the eight
Figure 4 contrasts, we read the existing full-data limma tables (results_v2.1.2_limma) and, over
genes present in both tables, computed the median of lfcSE_subtype / lfcSE_pooled (lfcSE = logFC/t,
the moderated standard error). We also report median |t| and the number of genes at BH FDR < 0.05
per table, for the subtype, the merged author type and our pipeline label.
Outputs: se_ratio_per_contrast.tsv, se_ratio_summary.tsv (mean and SD over contrasts).
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import CONTRASTS, DE, GROUPS, OUT, de_path  # noqa: E402


def load(scheme, ct, contrast):
    p = de_path(scheme, ct, contrast)
    if not os.path.exists(p):
        return None
    d = pd.read_csv(p, sep="\t", usecols=["gene", "logFC", "t", "adj.P.Val"]).set_index("gene")
    d["se"] = d["logFC"] / d["t"]
    return d


def main():
    rows = []
    for g, (label, subs) in GROUPS.items():
        for contrast, short in CONTRASTS.items():
            pooled = load("merged", label, contrast)
            ours = load("ours", label, contrast)
            for scheme, ct, d in [("merged", label, pooled), ("ours", label, ours)] + \
                    [("author", s, load("author", s, contrast)) for s in subs]:
                if d is None:
                    continue
                r = {"group": g, "contrast": short, "scheme": scheme, "celltype": ct,
                     "n_genes": len(d), "median_abs_t": float(np.median(np.abs(d["t"]))),
                     "n_fdr05": int((d["adj.P.Val"] < 0.05).sum())}
                if scheme != "merged":
                    sh = d.index.intersection(pooled.index)
                    ratio = d.loc[sh, "se"] / pooled.loc[sh, "se"]
                    r.update(n_shared=len(sh), median_se_ratio_vs_pooled=float(np.median(ratio)))
                rows.append(r)
    df = pd.DataFrame(rows)
    df.to_csv(f"{OUT}/se_ratio_per_contrast.tsv", sep="\t", index=False, float_format="%.4g")
    s = df.groupby(["group", "scheme", "celltype"], sort=False).agg(
        mean_se_ratio=("median_se_ratio_vs_pooled", "mean"), sd_se_ratio=("median_se_ratio_vs_pooled", "std"),
        mean_median_abs_t=("median_abs_t", "mean"), mean_n_fdr05=("n_fdr05", "mean"),
        n_contrasts=("contrast", "nunique")).reset_index()
    s.to_csv(f"{OUT}/se_ratio_summary.tsv", sep="\t", index=False, float_format="%.4g")
    print(s.to_string())


if __name__ == "__main__":
    main()
