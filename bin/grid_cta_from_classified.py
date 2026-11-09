#!/usr/bin/env python3
"""Cell type assignment (CTA) table from one nf-annotation-grid classified h5ad, in the layout of
the sc-annotation-pipeline <study>_<level>_cell_type.tsv that build_local_pseudobulks.py and
per_cohort_f1.py read (sample_id, cell_id = source h5ad obs_name, cell_type, cell_type_uri).

The grid was run on the full source h5ad (annotate_query.py resets obs_names to 0..n-1, which
equal the source obs_names). The script checks this: same cell count, and sample_id and the
author cell_id column identical row by row. Cells labelled 'unknown' are left out (no
assignment, as Gemma drops unassigned cells).
"""
import argparse
import json

import anndata as ad
import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--classified", required=True, help="grid *_annotated.h5ad")
    p.add_argument("--source_h5ad", required=True, help="h5ad the grid annotated (psychENCODE source)")
    p.add_argument("--level", default="subclass")
    p.add_argument("--out", required=True)
    a = p.parse_args()

    cl = ad.read_h5ad(a.classified, backed="r").obs
    src = ad.read_h5ad(a.source_h5ad, backed="r")
    so = src.obs
    if len(cl) != len(so):
        raise SystemExit(f"cell count differs: {len(cl)} vs {len(so)}")
    if not np.array_equal(so.index.astype(str), np.arange(len(so)).astype(str)):
        raise SystemExit("source obs_names are not 0..n-1")
    for col in ("sample_id", "cell_id"):
        if not np.array_equal(cl[col].astype(str).values, so[col].astype(str).values):
            raise SystemExit(f"{col} differs between classified and source h5ad")
    lab = cl[f"predicted_{a.level}"].astype(str).values
    cta = pd.DataFrame({"sample_id": so["sample_id"].astype(str).values,
                        "cell_id": so.index.astype(str), "cell_type": lab, "cell_type_uri": ""})
    unk = cta["cell_type"].str.lower().isin(["unknown", "nan", ""])
    cta[~unk].to_csv(a.out, sep="\t", index=False)
    print(json.dumps({"cells": int(len(cta)), "cells_unknown_dropped": int(unk.sum()),
                      "config_id": str(cl["annotation_config_id"].iloc[0]),
                      "label_counts": cta.loc[~unk, "cell_type"].value_counts().to_dict()}))


if __name__ == "__main__":
    main()
