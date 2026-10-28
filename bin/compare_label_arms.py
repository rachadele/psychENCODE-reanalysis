#!/usr/bin/env python3
"""Spearman correlation of DE results between the pipeline-label arm and the author-label arm,
for each author cell type and the pipeline cell type it maps to (label map JSON), per contrast.

Correlates log2FoldChange, the test statistic and pvalue over the genes tested in both arms.
The statistic is the `stat` column: the moderated t for limma-voom (bin/limma_voom_analysis.R)
and the Wald statistic log2FoldChange / lfcSE (signed) for DESeq2; when `stat` is missing it is
computed as log2FoldChange / lfcSE.

Input layout (both arms): <dir>/<cell_type>/<contrast>/<contrast>_<cell_type>_results.tsv
"""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

METRICS = ["log2FoldChange", "stat", "pvalue"]


def load(path):
    df = pd.read_csv(path, sep="\t")
    if "stat" not in df.columns:
        df["stat"] = df["log2FoldChange"] / df["lfcSE"]
    df = df.dropna(subset=["gene"]).drop_duplicates(subset=["gene"]).set_index("gene")
    return df[METRICS].replace([np.inf, -np.inf], np.nan)


def results(d):
    """{(cell_type, contrast): path}"""
    out = {}
    for f in glob.glob(os.path.join(d, "*", "*", "*_results.tsv")):
        contrast = os.path.basename(os.path.dirname(f))
        ct = os.path.basename(os.path.dirname(os.path.dirname(f)))
        out[(ct, contrast)] = f
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pavlab_dir", required=True, help="pipeline-label DE dir")
    p.add_argument("--author_dir", required=True, help="author-label DE dir")
    p.add_argument("--label_map", required=True,
                   help="author -> pipeline label JSON (or a params JSON with gemma_to_gemma_map)")
    p.add_argument("--label", default="", help="value of the 'run' column")
    p.add_argument("--min_genes", type=int, default=10)
    p.add_argument("--out", required=True, help="output TSV, one row per pair and contrast")
    a = p.parse_args()

    m = json.load(open(a.label_map))
    m = m.get("gemma_to_gemma_map", m)
    pav, aut = results(a.pavlab_dir), results(a.author_dir)
    rows = []
    for (act, contrast), apath in sorted(aut.items()):
        pct = m.get(act)
        ppath = pav.get((pct, contrast)) if pct else None
        row = dict(run=a.label, contrast=contrast, ct_author=act, ct_sc_pipeline=pct)
        if ppath is None:
            rows.append({**row, "n_genes": 0, "status": "no_pipeline_result"})
            continue
        merged = load(apath).join(load(ppath), how="inner", lsuffix="_author", rsuffix="_pipeline")
        row["n_genes"] = len(merged)
        for metric in METRICS:
            pair = merged[[f"{metric}_author", f"{metric}_pipeline"]].dropna()
            row[f"spearman_{metric}"] = (spearmanr(pair.iloc[:, 0], pair.iloc[:, 1])[0]
                                         if len(pair) >= a.min_genes else np.nan)
        row["status"] = "ok" if len(merged) >= a.min_genes else "too_few_genes"
        rows.append(row)
    pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)
    print(f"wrote {a.out}: {len(rows)} rows")


if __name__ == "__main__":
    main()
