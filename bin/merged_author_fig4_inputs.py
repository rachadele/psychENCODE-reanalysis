#!/usr/bin/env python3
"""Inputs of Figure 4 (merged-author arm) that fig_de_agreement.merged_figure() draws, written to
results_v2.1.2_limma/merged_author_2026-11-09/fig4_inputs/ (each step skipped when its output exists):
  1. celltype_f1_cells_class_merged.tsv: class F1 per merged type under the de_donors rule
     (class_f1_de_donors.py logic, author labels first mapped to merged pipeline subclass labels;
     every merged type has one class label);
  2. null/: benchmark-post-hoc/de_permutation_null de_permutation_null.run_test ("labels" test only,
     same seed and permutation counts) on the 168 merged pairs, and matched_effect_adjusted.main, both
     imported unchanged with their input and output paths pointed at the merged arm;
     celltype_f1_cells_merged.tsv (merged author rows plus the current pipeline rows) feeds the latter.
"""
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

REPO = "/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis"
R = f"{REPO}/results_v2.1.2_limma"
S = f"{R}/merged_author_2026-11-09"
NULL_REPO = "/space/grp/rschwartz/rschwartz/benchmark-post-hoc/de_permutation_null"
D = f"{S}/fig4_inputs"
sys.path.insert(0, f"{REPO}/bin")

from de_agreement_vs_annotation import de_donor_table  # noqa: E402

SUBMAP = json.load(open(f"{REPO}/assets/author_to_pipeline_map.subclass.json"))
CLASSMAP = json.load(open(f"{REPO}/assets/author_to_pipeline_map.class.json"))
MEMBERS = {}
for t, lab in SUBMAP.items():
    MEMBERS.setdefault(lab, []).append(t)
MERGED_CC = f"{S}/pseudobulks/author_merged/ct_pseudobulks/cell_counts.tsv"
STUDIES = open(f"{REPO}/shortnames.txt").read().split()


def class_f1():
    out = f"{D}/celltype_f1_cells_class_merged.tsv"
    if os.path.exists(out):
        return out
    from per_cohort_f1 import load_cells, per_sample_f1
    cmap = {}
    for lab, ms in MEMBERS.items():
        cs = {CLASSMAP[m] for m in ms}
        assert len(cs) == 1, (lab, cs)
        cmap[lab] = cs.pop()
    pc = f"{D}/per_cohort_class"
    os.makedirs(pc, exist_ok=True)
    h5 = "/space/grp/rschwartz/rschwartz/get_gemma_data.nf/psychencode.txt_author_false_process_samples_false/h5ad"
    ann = f"{REPO}/results_v2.1.2/annotation"
    for s in STUDIES:
        df, _ = load_cells(f"{h5}/{s}/{s}.h5ad", f"{ann}/{s}/predicted_celltypes/{s}_class_cell_type.tsv",
                           f"{ann}/{s}/qc/{s}_mask.tsv")
        df["author"] = df["author"].map(SUBMAP)
        assert df["author"].notna().all()
        ps = per_sample_f1(df, cmap)
        ps.insert(0, "study", s)
        ps.to_csv(f"{pc}/{s}_per_sample_f1.tsv", sep="\t", index=False, float_format="%.6g")
        print("class F1", s, file=sys.stderr)
    tab = de_donor_table(pc, f"{REPO}/gemma/metadata", {"author": MERGED_CC}, cmap)
    tab.insert(0, "level", "class")
    tab["class_label"] = tab["celltype"].map(cmap)
    tab.to_csv(out, sep="\t", index=False, float_format="%.6g")
    return out


def f1_cells_table():
    """celltype_f1_cells-style table: merged author rows + current pipeline rows (matched_effect_adjusted input)."""
    out = f"{D}/celltype_f1_cells_merged.tsv"
    if not os.path.exists(out):
        ident = {lab: lab for lab in MEMBERS}
        per_sample = f"{S}/per_cohort_f1"
        au = de_donor_table(per_sample, f"{REPO}/gemma/metadata", {"author": MERGED_CC}, ident)
        cur = pd.read_csv(f"{R}/de_agreement_vs_annotation/celltype_f1_cells.tsv", sep="\t")
        pd.concat([au[au["scheme"] == "author"], cur[cur["scheme"] != "author"]], ignore_index=True) \
          .to_csv(out, sep="\t", index=False, float_format="%.6g")
    return out


def null():
    nd = f"{D}/null"
    os.makedirs(nd, exist_ok=True)
    sys.path.insert(0, NULL_REPO)
    import de_permutation_null as dpn
    if not os.path.exists(f"{nd}/null_summary.tsv"):
        adir = f"{S}/author_submitted_true/subclass/limma_voom/gemma"
        pdir = f"{R}/author_submitted_false/subclass/limma_voom/gemma"

        def label_tables():
            out = {}
            for c in dpn.LABEL_CONTRASTS:
                auth, pipe, match = {}, {}, {}
                for lab in MEMBERS:
                    ap = f"{adir}/{lab}/{c}/{c}_{lab}_results.tsv"
                    pp = f"{pdir}/{lab}/{c}/{c}_{lab}_results.tsv"
                    if not (os.path.exists(ap) and os.path.exists(pp)):
                        continue
                    auth[lab] = dpn.read_our(ap)[["gene", "log2FoldChange"]]
                    pipe[lab] = dpn.read_our(pp)[["gene", "log2FoldChange"]]
                    match[lab] = lab
                out[c] = (auth, pipe, match)
            return out

        rng = np.random.default_rng(dpn.SEED)
        pairs, summ, nulls = dpn.run_test("labels", "subclass", label_tables(), dpn.pair_labels, rng)
        # matched rho must equal the merged arm's comparison table
        src = f"{S}/comparison/subclass/average_corr/combined_filtered_corr.tsv"
        ref = pd.read_csv(src, sep="\t")
        ref = ref[ref["ct_author"] == ref["ct_sc_pipeline"]].set_index(["contrast", "ct_author"])["spearman_log2FoldChange"]
        got = pairs[pairs["matched"]].set_index(["contrast", "celltype"])["rho"]
        j = pd.concat([got.rename("g"), ref.rename("r")], axis=1, join="inner")
        assert len(j) == len(got) and np.allclose(j["g"], j["r"], atol=1e-8)
        pairs.to_csv(f"{nd}/pairs.tsv", sep="\t", index=False, float_format="%.6g")
        pd.DataFrame(summ).to_csv(f"{nd}/null_summary.tsv", sep="\t", index=False, float_format="%.6g")
        nulls.to_csv(f"{nd}/null_means.tsv.gz", sep="\t", index=False, float_format="%.6g")
        json.dump({"script": os.path.abspath(__file__), "null_module_sha256": dpn.sha(dpn.__file__),
                   "seed": dpn.SEED, "n_perm_celltype": dpn.N_PERM_CT, "n_perm_gene": dpn.N_PERM_GENE,
                   "inputs": {p: dpn.sha(p) for p in sorted(dpn.INPUTS)},
                   "outputs": {o: dpn.sha(f"{nd}/{o}") for o in ("pairs.tsv", "null_summary.tsv", "null_means.tsv.gz")}},
                  open(f"{nd}/manifest.json", "w"), indent=1)
    if not os.path.exists(f"{nd}/matched_effect_adjusted.tsv"):
        import matched_effect_adjusted as mea
        mea.HERE, mea.PAIRS, mea.TABLE = nd, f"{nd}/pairs.tsv", f1_cells_table()
        mea.COUNTS = dict(mea.COUNTS, author=MERGED_CC)
        mea.main()
    return nd


def main():
    os.makedirs(D, exist_ok=True)
    class_f1()
    null()
    print(pd.read_csv(f"{D}/null/null_summary.tsv", sep="\t").round(4).to_string(index=False))


if __name__ == "__main__":
    main()
