#!/usr/bin/env python3
"""Build a per-study pseudobulk matrix from a local h5ad and sc-annotation-pipeline
cell type assignments, in the same layout that
`gemma-cli getSingleCellDataMatrix --aggregate-by-assay --aggregate-by-cell-type-assignment`
writes, so that bin/aggregate_celltypes_gemma.py and bin/DESeq2_analysis.R (gemma mode)
consume it unchanged.

What it does:
  * sums raw counts (h5ad X) over cells per (sample, cell type) -- Gemma's default SUM;
  * drops cells with no assignment (Gemma's default, aggregateUnknownCharacteristics off);
  * optionally drops cells flagged by the pipeline QC mask (--mask). The Gemma CLI path used
    by GET_GEMMA_PSEUDOBULKS applies no mask, so the default (no --mask) matches it;
  * writes one column per sample x cell type, named `<sample>...<ct>___<sample>...<ct>`
    with R-friendly names (non [A-Za-z0-9_.] characters -> '.'), as Gemma does;
  * rows: Probe = NCBI gene id and Sequence = "SYMBOL [NCBI ID=..] ..." taken by position
    from the study's MEX features.tsv.gz (same gene order as the h5ad) and, when given,
    the symbol strings of an earlier Gemma pseudobulk file for the same study.
Also writes <study>_cell_counts.tsv (cells per sample x cell type, kept and masked).
"""
import argparse
import gzip
import re

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def rname(s):
    return re.sub(r"[^A-Za-z0-9_.]", ".", str(s))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True)
    p.add_argument("--cta", required=True, help="pipeline <study>_<level>_cell_type.tsv")
    p.add_argument("--mex_features", required=True, help="features.tsv.gz of one MEX sample")
    p.add_argument("--gemma_rows", default=None, help="earlier Gemma pseudobulk tsv.gz (Probe/Sequence)")
    p.add_argument("--mask", default=None, help="pipeline <study>_mask.tsv; exclude cells with value true")
    p.add_argument("--study_name", required=True)
    p.add_argument("--chunk", type=int, default=50000)
    a = p.parse_args()

    adata = ad.read_h5ad(a.h5ad, backed="r")
    obs = adata.obs[["sample_id", "sample_name"]].copy()
    obs["sample_id"] = obs["sample_id"].astype(str)
    obs["row_id"] = adata.obs_names.astype(str)

    cta = pd.read_csv(a.cta, sep="\t", dtype=str)
    # process_query.py sets cell_id to the h5ad obs_names
    if not set(cta["cell_id"]).issubset(set(obs["row_id"])):
        raise SystemExit("CTA cell_id values do not match h5ad obs_names")
    obs = obs.merge(cta[["sample_id", "cell_id", "cell_type"]], how="left",
                    left_on=["sample_id", "row_id"], right_on=["sample_id", "cell_id"])
    if len(obs) != adata.n_obs:
        raise SystemExit("join changed cell count")
    obs["masked"] = False
    if a.mask:
        m = pd.read_csv(a.mask, sep="\t", dtype=str)
        m = m[m["value"].str.lower() == "true"]
        key = set(zip(m["sample_id"], m["cell_id"]))
        obs["masked"] = [(s, c) in key for s, c in zip(obs["sample_id"], obs["row_id"])]

    obs["sample"] = obs["sample_name"].astype(str).map(rname)
    obs["ct"] = obs["cell_type"].map(lambda x: rname(x) if isinstance(x, str) else np.nan)

    counts = (obs.groupby(["sample", "ct"], dropna=False)
                 .agg(n_cells=("row_id", "size"), n_masked=("masked", "sum")).reset_index())
    counts.insert(0, "study", a.study_name)
    counts.to_csv(f"{a.study_name}_cell_counts.tsv", sep="\t", index=False)

    keep = obs["ct"].notna() & ~obs["masked"]
    samples = sorted(obs["sample"].unique())
    cts = sorted(obs.loc[obs["ct"].notna(), "ct"].unique())
    groups = [(s, c) for s in samples for c in cts]
    gidx = {g: i for i, g in enumerate(groups)}
    cell_group = np.array([gidx[(s, c)] if k else -1
                           for s, c, k in zip(obs["sample"], obs["ct"], keep)])

    out = np.zeros((len(groups), adata.n_vars))
    for start in range(0, adata.n_obs, a.chunk):
        stop = min(start + a.chunk, adata.n_obs)
        X = sp.csr_matrix(adata.X[start:stop])
        g = cell_group[start:stop]
        sel = g >= 0
        ind = sp.csr_matrix((np.ones(sel.sum()), (g[sel], np.nonzero(sel)[0])),
                            shape=(len(groups), stop - start))
        out += (ind @ X).toarray()

    feats = pd.read_csv(a.mex_features, sep="\t", header=None, dtype=str)
    if len(feats) != adata.n_vars:
        raise SystemExit("MEX features and h5ad var differ in length")
    probe = feats[0].values
    sym_h5ad = adata.var["feature_name"].astype(object).values
    seq = {}
    if a.gemma_rows:
        g = pd.read_csv(a.gemma_rows, sep="\t", comment="#", usecols=[0, 1], dtype=str)
        seq = dict(zip(g["Probe"], g["Sequence"]))
    sequence = [seq.get(pr) or f"{s if isinstance(s, str) else pr} [NCBI ID={pr}]"
                for pr, s in zip(probe, sym_h5ad)]

    cols = [f"{s}...{c}___{s}...{c}" for s, c in groups]
    df = pd.DataFrame(out.T, columns=cols)
    df.insert(0, "Sequence", sequence)
    df.insert(0, "Probe", probe)
    with gzip.open(f"{a.study_name}_aggregated_gemma.tsv.gz", "wt") as fh:
        fh.write(f"# Local pseudobulk for {a.study_name}: SUM over cells by sample and cell type\n")
        fh.write(f"# h5ad: {a.h5ad}\n# cta: {a.cta}\n# mask: {a.mask}\n")
        df.to_csv(fh, sep="\t", index=False)
    print(f"{a.study_name}: {adata.n_obs} cells, {keep.sum()} kept, {len(samples)} samples, "
          f"{len(cts)} cell types, {int((obs['ct'].isna()).sum())} unassigned, "
          f"{int(obs['masked'].sum())} masked")


if __name__ == "__main__":
    main()
