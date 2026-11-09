#!/usr/bin/env python3
"""Per-sample subclass precision, recall and F1 of the deployed pipeline labels against
author labels, for one cohort, then the mean over samples per cell type.

Cells: QC-masked cells and cells without a pipeline label are dropped (as in
de_agreement_vs_annotation.py and build_local_pseudobulks.py). Author labels come from
cell_id (<label>[.<n>]_<sample_name>). An author type t maps to pipeline label L = map[t].
Within each sample:
  recall(t)    = cells of t labelled L / cells of t            (NaN if the sample has no t cells)
  precision(t) = cells labelled L whose author type maps to L / cells labelled L
                 (NaN if no cell is labelled L; types pooled into L share it)
  F1(t)        = 2PR / (P + R); 0 when recall is 0; NaN when the sample has no t cells.
The per-cohort value is the unweighted mean over samples with at least one t cell (each
sample counts once; cells are never pooled across samples).
Writes <outdir>/<study>_per_sample_f1.tsv and <outdir>/<study>_per_cohort_f1.tsv.
"""
import argparse
import json
import os
import re

import anndata as ad
import numpy as np
import pandas as pd


def rname(s):
    return re.sub(r"[^A-Za-z0-9_.]", ".", str(s))


def load_cells(h5ad, cta, mask):
    obs = ad.read_h5ad(h5ad, backed="r").obs
    df = pd.DataFrame({"sample_id": obs["sample_id"].astype(str).values,
                       "row_id": obs.index.astype(str),
                       "cell_id": obs["cell_id"].astype(str).values,
                       "sample_name": obs["sample_name"].astype(str).values})
    df["author"] = [rname(re.sub(r"\.\d+$", "", c[:-len(s) - 1])) for c, s in zip(df["cell_id"], df["sample_name"])]
    c = pd.read_csv(cta, sep="\t", dtype=str)
    c["pipeline"] = c["cell_type"].map(rname)
    df = df.merge(c[["sample_id", "cell_id", "pipeline"]].rename(columns={"cell_id": "row_id"}),
                  on=["sample_id", "row_id"], how="left")
    m = pd.read_csv(mask, sep="\t", dtype=str)
    m = m[m["value"].str.lower() == "true"]
    masked = pd.MultiIndex.from_frame(m[["sample_id", "cell_id"]])
    df["masked"] = pd.MultiIndex.from_frame(df[["sample_id", "row_id"]]).isin(masked)
    n_all = len(df)
    df = df[~df["masked"] & df["pipeline"].notna()].copy()
    return df, {"cells_total": n_all, "cells_excluded_masked_or_unlabelled": int(n_all - len(df)),
                "cells_used": int(len(df))}


def per_sample_f1(df, label_map):
    df = df.assign(author_mapped=df["author"].map(label_map))
    rows = []
    for s, g in df.groupby("sample_name", sort=True):
        n_author = g["author"].value_counts()
        n_pred = g["pipeline"].value_counts()
        for t, lab in label_map.items():
            nt = int(n_author.get(t, 0))
            if nt == 0:
                continue
            npred = int(n_pred.get(lab, 0))
            tp_t = int(((g["author"] == t) & (g["pipeline"] == lab)).sum())
            recall = tp_t / nt
            if npred:
                precision = int(((g["author_mapped"] == lab) & (g["pipeline"] == lab)).sum()) / npred
            else:
                precision = np.nan
            if recall == 0:
                f1 = 0.0
            else:
                f1 = 2 * precision * recall / (precision + recall)
            rows.append({"sample": s, "ct_author": t, "pipeline_label": lab, "n_author_cells": nt,
                         "n_pipeline_cells": npred, "precision": precision, "recall": recall, "f1": f1})
    return pd.DataFrame(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--h5ad", required=True)
    p.add_argument("--cta", required=True)
    p.add_argument("--mask", required=True)
    p.add_argument("--label_map", required=True)
    p.add_argument("--study", required=True)
    p.add_argument("--outdir", required=True)
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    label_map = json.load(open(a.label_map))
    label_map = label_map.get("gemma_to_gemma_map", label_map)
    df, info = load_cells(a.h5ad, a.cta, a.mask)
    unmapped = sorted(set(df["author"]) - set(label_map))
    ps = per_sample_f1(df, label_map)
    ps.insert(0, "study", a.study)
    ps.to_csv(os.path.join(a.outdir, f"{a.study}_per_sample_f1.tsv"), sep="\t", index=False, float_format="%.6g")
    pc = (ps.groupby(["study", "ct_author", "pipeline_label"])
          .agg(precision=("precision", "mean"), recall=("recall", "mean"), f1=("f1", "mean"),
               n_samples=("f1", "size"), n_samples_precision=("precision", "count"),
               median_author_cells_in_samples=("n_author_cells", "median"))
          .reset_index())
    pc.to_csv(os.path.join(a.outdir, f"{a.study}_per_cohort_f1.tsv"), sep="\t", index=False, float_format="%.6g")
    info.update({"study": a.study, "n_samples": int(df["sample_name"].nunique()),
                 "author_labels_not_in_map": unmapped,
                 "cells_with_unmapped_author_label": int(df["author"].isin(unmapped).sum())})
    json.dump(info, open(os.path.join(a.outdir, f"{a.study}_cells.json"), "w"), indent=2)
    print(json.dumps(info))


if __name__ == "__main__":
    main()
