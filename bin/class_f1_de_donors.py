#!/usr/bin/env python3
"""Class-level F1 of the pipeline labels per author cell type, under the de_donors rule
(docs/DE_AGREEMENT_DECISIONS.md): per-sample F1 over every donor in the DE of the author type,
then the plain mean (cells never pooled, no weighting).

Same code path as celltype_f1_cells.tsv: per_cohort_f1.load_cells / per_sample_f1 per cohort,
then de_agreement_vs_annotation.de_donor_table. The author type t maps to its class through
assets/author_to_pipeline_map.class.json; pipeline labels are the pipeline's class-level labels
(<study>_class_cell_type.tsv). Author types that share a class share its precision.

Writes <outdir>/per_cohort_class/<study>_per_sample_f1.tsv, <outdir>/celltype_f1_cells_class.tsv
and <outdir>/manifest_class_f1.json.
"""
import argparse
import datetime
import hashlib
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from per_cohort_f1 import load_cells, per_sample_f1  # noqa: E402
from de_agreement_vs_annotation import de_donor_table  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STUDIES = ["CMC", "DevBrain", "MultiomeBrain", "PTSDBrainomics", "SZBDMulti-Seq", "UCLA-ASD",
           "Velmeshev-2019.1", "Velmeshev-2019.2"]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--h5ad_dir", default="/space/grp/rschwartz/rschwartz/get_gemma_data.nf/"
                   "psychencode.txt_author_false_process_samples_false/h5ad")
    p.add_argument("--annotation_dir", default=os.path.join(ROOT, "results_v2.1.2", "annotation"))
    p.add_argument("--label_map", default=os.path.join(ROOT, "assets", "author_to_pipeline_map.class.json"))
    p.add_argument("--author_cell_counts", default=os.path.join(
        ROOT, "results_v2.1.2_limma", "pseudobulks", "author", "ct_pseudobulks", "cell_counts.tsv"))
    p.add_argument("--meta_dir", default=os.path.join(ROOT, "gemma", "metadata"))
    p.add_argument("--outdir", default=os.path.join(ROOT, "results_v2.1.2_limma", "de_agreement_vs_annotation"))
    a = p.parse_args()

    label_map = json.load(open(a.label_map))
    pc_dir = os.path.join(a.outdir, "per_cohort_class")
    os.makedirs(pc_dir, exist_ok=True)
    inputs, cells = [a.label_map, a.author_cell_counts], {}
    for s in STUDIES:
        h5ad = os.path.join(a.h5ad_dir, s, f"{s}.h5ad")
        cta = os.path.join(a.annotation_dir, s, "predicted_celltypes", f"{s}_class_cell_type.tsv")
        mask = os.path.join(a.annotation_dir, s, "qc", f"{s}_mask.tsv")
        df, info = load_cells(h5ad, cta, mask)
        unmapped = sorted(set(df["author"]) - set(label_map))
        info["author_labels_not_in_map"] = unmapped
        cells[s] = info
        ps = per_sample_f1(df, label_map)
        ps.insert(0, "study", s)
        ps.to_csv(os.path.join(pc_dir, f"{s}_per_sample_f1.tsv"), sep="\t", index=False, float_format="%.6g")
        inputs += [cta, mask]
        print(s, info, file=sys.stderr)

    tab = de_donor_table(pc_dir, a.meta_dir, {"author": a.author_cell_counts}, label_map)
    tab.insert(0, "level", "class")
    tab["class_label"] = tab["celltype"].map(label_map)
    out = os.path.join(a.outdir, "celltype_f1_cells_class.tsv")
    tab.to_csv(out, sep="\t", index=False, float_format="%.6g")

    now = datetime.datetime.now().isoformat(timespec="seconds")
    man = {"created": now, "producer": os.path.abspath(__file__),
           "rule": "de_donors: plain mean of per-sample F1 over the donors in the DE of each author type",
           "data_version": "pipeline labels 2.1.2 (class level), DE results_v2.1.2_limma",
           "cells": cells,
           "inputs": {f: {"sha256": sha256(f), "mtime": datetime.datetime.fromtimestamp(os.path.getmtime(f)).isoformat(timespec="seconds")}
                      for f in inputs},
           "inputs_not_hashed": {"h5ad_dir (obs only; size)": a.h5ad_dir},
           "outputs": {out: {"sha256": sha256(out), "written": now}}}
    with open(os.path.join(a.outdir, "manifest_class_f1.json"), "w") as f:
        json.dump(man, f, indent=2)
    print(out)


if __name__ == "__main__":
    main()
