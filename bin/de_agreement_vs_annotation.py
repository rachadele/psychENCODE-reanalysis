#!/usr/bin/env python3
"""DE agreement (author vs pipeline labels) against annotation quality and support.

Per author subclass of one study (default CMC):
  * support: cells per donor in the QC-masked pseudobulks (median over all donors of the
    study, zeros included, and total), for the author-label pseudobulk of the type and for
    the pipeline-label pseudobulk of the label it maps to;
  * CMC F1: F1 of the pipeline labels against the author labels on the cells used in the
    DE (QC-masked cells excluded). An author type t maps to pipeline label L = map[t];
    recall = share of t cells labelled L; precision = share of L cells whose author type maps
    to L (so author types that share L share its precision);
  * mean Spearman rho of log2 fold changes across contrasts (as in the heatmap);
  * benchmark F1 of L (extract_benchmark_f1.py).
Then Spearman correlations across cell types of mean rho with benchmark F1, CMC F1 and
log10 support, and partial Spearman of rho with each F1 controlling for log10 support.
Writes summary TSVs, the combined figure and a provenance JSON to --outdir.
"""
import argparse
import datetime
import hashlib
import json
import os
import re

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

from utils import (author_display, contrast_pivot, load_benchmark_f1, load_label_map,
                   own_label_types, plot_correlation_heatmap)


def rname(s):
    return re.sub(r"[^A-Za-z0-9_.]", ".", str(s))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = ~np.isnan(x) & ~np.isnan(y)
    r, p = stats.spearmanr(x[ok], y[ok])
    return r, p, int(ok.sum())


def partial_spearman(x, y, z):
    """Spearman of x and y controlling for z: Pearson of rank residuals, t test on n-3 df."""
    x, y, z = (np.asarray(v, float) for v in (x, y, z))
    ok = ~np.isnan(x) & ~np.isnan(y) & ~np.isnan(z)
    rx, ry, rz = (stats.rankdata(v[ok]) for v in (x, y, z))
    res = lambda a: a - np.polyval(np.polyfit(rz, a, 1), rz)
    r = np.corrcoef(res(rx), res(ry))[0, 1]
    n = int(ok.sum())
    t = r * np.sqrt((n - 3) / (1 - r ** 2))
    return r, 2 * stats.t.sf(abs(t), n - 3), n


def cmc_f1(h5ad, cta, mask, label_map):
    obs = ad.read_h5ad(h5ad, backed="r").obs
    df = pd.DataFrame({"sample_id": obs["sample_id"].astype(str).values,
                       "row_id": obs.index.astype(str),
                       "cell_id": obs["cell_id"].astype(str).values,
                       "sample_name": obs["sample_name"].astype(str).values})
    # author label from cell_id (<label>[.<n>]_<sample_name>), as build_local_pseudobulks.py
    df["author"] = [rname(re.sub(r"\.\d+$", "", c[:-len(s) - 1])) for c, s in zip(df["cell_id"], df["sample_name"])]
    c = pd.read_csv(cta, sep="\t", dtype=str)
    c["pipeline"] = c["cell_type"].map(rname)
    df = df.merge(c[["sample_id", "cell_id", "pipeline"]].rename(columns={"cell_id": "row_id"}),
                  on=["sample_id", "row_id"], how="left")
    m = pd.read_csv(mask, sep="\t", dtype=str)
    m = m[m["value"].str.lower() == "true"]
    masked = set(zip(m["sample_id"], m["cell_id"]))
    df["masked"] = [(s, r) in masked for s, r in zip(df["sample_id"], df["row_id"])]
    n_all = len(df)
    df = df[~df["masked"] & df["pipeline"].notna()]
    df["author_mapped"] = df["author"].map(label_map)
    rows = []
    for t, lab in label_map.items():
        is_t = df["author"] == t
        if is_t.sum() == 0:
            continue
        pred_l = df["pipeline"] == lab
        recall = (is_t & pred_l).sum() / is_t.sum()
        precision = ((df["author_mapped"] == lab) & pred_l).sum() / max(pred_l.sum(), 1)
        f1 = 0 if precision + recall == 0 else 2 * precision * recall / (precision + recall)
        rows.append({"ct_author": t, "pipeline_label": lab, "n_author_cells": int(is_t.sum()),
                     "n_pipeline_cells": int(pred_l.sum()), "cmc_precision": precision,
                     "cmc_recall": recall, "cmc_f1": f1})
    return pd.DataFrame(rows), {"cells_total": n_all, "cells_excluded_masked_or_unlabelled": int(n_all - len(df)),
                                "cells_used": int(len(df))}


def support(cell_counts, study, key):
    cc = pd.read_csv(cell_counts, sep="\t")
    cc = cc[cc["study"] == study]
    wide = cc.pivot_table(index="sample", columns="ct", values="n_kept", aggfunc="sum").fillna(0)
    return pd.DataFrame({f"{key}_median_cells_per_donor": wide.median(),
                         f"{key}_total_cells": wide.sum(),
                         f"{key}_donors_with_cells": (wide > 0).sum(),
                         f"{key}_n_donors": len(wide)})


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--combined_corr", required=True, help="subclass combined_filtered_corr.tsv")
    p.add_argument("--params", required=True)
    p.add_argument("--label_map", required=True, help="author -> pipeline subclass map JSON")
    p.add_argument("--benchmark_f1", required=True)
    p.add_argument("--h5ad", required=True)
    p.add_argument("--cta", required=True, help="<study>_subclass_cell_type.tsv (pipeline v2.1.2)")
    p.add_argument("--mask", required=True, help="<study>_mask.tsv")
    p.add_argument("--author_cell_counts", required=True, help="pseudobulks/author/ct_pseudobulks/cell_counts.tsv")
    p.add_argument("--pipeline_cell_counts", required=True, help="pseudobulks/pavlab_subclass/ct_pseudobulks/cell_counts.tsv")
    p.add_argument("--study", default="CMC")
    p.add_argument("--level", default="subclass")
    p.add_argument("--metric", default="spearman_log2FoldChange")
    p.add_argument("--outdir", required=True)
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    label_map = load_label_map(a.params, a.label_map)
    combined = pd.read_csv(a.combined_corr, sep="\t")
    pivot = contrast_pivot(combined, label_map, a.metric)
    tab = pd.DataFrame({"ct_author": pivot.columns, "display_name": [author_display(c) for c in pivot.columns],
                        "pipeline_label": [label_map[c] for c in pivot.columns],
                        "mean_rho": pivot.mean().values, "n_contrasts": pivot.notna().sum().values})
    a2l = dict(zip(tab["ct_author"], tab["pipeline_label"]))
    bf1 = load_benchmark_f1(a.benchmark_f1, a.level, a2l)
    own = own_label_types(a2l, a.level)
    tab["benchmark_f1"] = tab["ct_author"].map(bf1)
    tab["benchmark_own_label"] = tab["ct_author"].isin(own)

    f1_df, cell_info = cmc_f1(a.h5ad, a.cta, a.mask, label_map)
    tab = tab.merge(f1_df.drop(columns="pipeline_label"), on="ct_author", how="left")
    sa = support(a.author_cell_counts, a.study, "author")
    sp = support(a.pipeline_cell_counts, a.study, "pipeline")
    tab = tab.merge(sa, left_on="ct_author", right_index=True, how="left")
    tab = tab.merge(sp, left_on="pipeline_label", right_index=True, how="left")
    tab["log10_author_median_cells"] = np.log10(tab["author_median_cells_per_donor"] + 1)
    tab["log10_pipeline_median_cells"] = np.log10(tab["pipeline_median_cells_per_donor"] + 1)
    tab = tab.sort_values("mean_rho", ascending=False)
    tab.to_csv(os.path.join(a.outdir, "celltype_summary.tsv"), sep="\t", index=False)

    y = tab["mean_rho"]
    own_mask = tab["benchmark_own_label"]
    logs = tab["log10_author_median_cells"]
    tests = []
    def add(name, subset, fn, *cols):
        r, pv, n = fn(*[c[subset] for c in cols])
        tests.append({"test": name, "subset": "own benchmark label" if subset is own_mask else "all",
                      "spearman_rho": r, "p": pv, "n": n})
    allm = pd.Series(True, index=tab.index)
    add("rho ~ benchmark F1", own_mask, spearman, tab["benchmark_f1"], y)
    add("rho ~ benchmark F1", allm, spearman, tab["benchmark_f1"], y)
    add("rho ~ CMC F1", allm, spearman, tab["cmc_f1"], y)
    add("rho ~ log10 median author cells per donor", allm, spearman, logs, y)
    add("rho ~ log10 median pipeline-label cells per donor", allm, spearman, tab["log10_pipeline_median_cells"], y)
    add("benchmark F1 ~ log10 median author cells", own_mask, spearman, tab["benchmark_f1"], logs)
    add("CMC F1 ~ log10 median author cells", allm, spearman, tab["cmc_f1"], logs)
    add("benchmark F1 ~ CMC F1", own_mask, spearman, tab["benchmark_f1"], tab["cmc_f1"])
    add("rho ~ benchmark F1 | log10 support (partial)", own_mask, partial_spearman, tab["benchmark_f1"], y, logs)
    add("rho ~ benchmark F1 | log10 support (partial)", allm, partial_spearman, tab["benchmark_f1"], y, logs)
    add("rho ~ CMC F1 | log10 support (partial)", allm, partial_spearman, tab["cmc_f1"], y, logs)
    pd.DataFrame(tests).to_csv(os.path.join(a.outdir, "correlations.tsv"), sep="\t", index=False)

    cmc = dict(zip(tab["ct_author"], tab["cmc_f1"]))
    logd = dict(zip(tab["ct_author"], logs))
    lo = min(0.4, np.floor(np.nanmin(tab["cmc_f1"]) * 10) / 10)
    plot_correlation_heatmap(
        combined, label_map, metric=a.metric, f1_path=a.benchmark_f1, annotation_level=a.level,
        output_path=os.path.join(a.outdir, f"fig_de_agreement_{a.level}.png"),
        extra_bars=[(f"{a.study} F1", cmc)], f1_range=(lo, 1.0),
        scatter_panels=[
            dict(x=cmc, xlabel=f"Annotation F1 ({a.study}, pipeline vs author labels)", label="outliers",
                 xlim=(lo - 0.03, 1.03)),
            dict(x=logd, xlabel=f"log10(median {a.study} cells per donor + 1), author label", label="outliers"),
        ])

    prov = {"created": datetime.datetime.now().astimezone().isoformat(), "study": a.study, "level": a.level,
            "cells": cell_info, "inputs": {}}
    for k in ("combined_corr", "label_map", "benchmark_f1", "cta", "mask", "author_cell_counts", "pipeline_cell_counts"):
        f = getattr(a, k)
        prov["inputs"][k] = {"path": os.path.abspath(f), "sha256": sha256(f),
                             "mtime": datetime.datetime.fromtimestamp(os.path.getmtime(f)).astimezone().isoformat()}
    prov["inputs"]["h5ad"] = {"path": os.path.abspath(a.h5ad),
                              "mtime": datetime.datetime.fromtimestamp(os.path.getmtime(a.h5ad)).astimezone().isoformat()}
    bj = os.path.splitext(a.benchmark_f1)[0] + ".json"
    if os.path.exists(bj):
        prov["benchmark_f1_provenance"] = json.load(open(bj))
    prov["data_version"] = {"pipeline_labels": "sc-annotation-pipeline 2.1.2", "de": "results_v2.1.2_limma"}
    json.dump(prov, open(os.path.join(a.outdir, "manifest.json"), "w"), indent=2)
    print(pd.DataFrame(tests).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
