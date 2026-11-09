#!/usr/bin/env python3
"""manifest.json for the merged-author arm: code commit, sha256 of inputs, pseudobulks, DE results,
comparison and agreement outputs, with modification times."""
import datetime, glob, hashlib, json, os, subprocess
REPO = "/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis"
R = f"{REPO}/results_v2.1.2_limma"
S = f"{R}/merged_author_2026-11-09"
def sha(f):
    h = hashlib.sha256()
    with open(f, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()
ts = lambda f: datetime.datetime.fromtimestamp(os.path.getmtime(f)).astimezone().isoformat(timespec="seconds")
rec = lambda fs: {os.path.relpath(f, REPO): {"sha256": sha(f), "mtime": ts(f)} for f in sorted(fs)}
git = lambda *a: subprocess.check_output(["git", "-C", REPO, *a], text=True).strip()
studies = open(f"{REPO}/shortnames.txt").read().split()
ann = f"{REPO}/results_v2.1.2/annotation"
m = {
 "created": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
 "data_version": "merged_author_2026-11-09: author labels mapped to pipeline subclass labels before pseudobulking; "
                 "pipeline labels sc-annotation-pipeline 2.1.2; pipeline-subclass DE from results_v2.1.2_limma",
 "status": "pending Rachel's decision (docs/DE_AGREEMENT_DECISIONS.md, 2026-11-09)",
 "repo": {"branch": git("rev-parse", "--abbrev-ref", "HEAD"), "commit": git("rev-parse", "HEAD")},
 "tmux_session": "merged_author_de",
 "nextflow_work": "/space/scratch/rschwartz_psyrerun_work/merged_author_2026-11-09",
 "de": {"design": "~ Disorder + Age_death + PMI + Biological_Sex + X1000G_ancestry (one fit per type over all cohorts, no cohort term)",
        "method": "limma_voom, pipeline defaults as results_v2.1.2_limma/scripts/run_de.sh"},
 "inputs": rec([f"{REPO}/assets/author_to_pipeline_map.subclass.json", f"{REPO}/params.subclass.json",
                f"{R}/comparison/subclass/average_corr/combined/combined_filtered_corr.tsv",
                f"{R}/de_agreement_vs_annotation/correlations_per_cohort.tsv",
                f"{R}/de_agreement_vs_annotation/celltype_f1_cells.tsv"]
               + [f"{ann}/{c}/qc/{c}_mask.tsv" for c in studies]
               + [f"{ann}/{c}/predicted_celltypes/{c}_subclass_cell_type.tsv" for c in studies]
               + glob.glob(f"{REPO}/gemma/metadata/*_sample_meta.tsv")
               + glob.glob(f"{R}/author_submitted_false/subclass/limma_voom/gemma/*/*/*_results.tsv")),
 "outputs": rec(glob.glob(f"{S}/pseudobulks/author_merged/ct_pseudobulks/*/*.tsv.gz")
                + [f"{S}/pseudobulks/author_merged/ct_pseudobulks/cell_counts.tsv", f"{S}/identity_map.subclass.json"]
                + glob.glob(f"{S}/author_submitted_true/subclass/limma_voom/gemma/*/*/*_results.tsv")
                + glob.glob(f"{S}/comparison/subclass/all_corr/*/*_pairwise_corr.tsv")
                + glob.glob(f"{S}/comparison/subclass/average_corr/*.tsv")
                + glob.glob(f"{S}/per_cohort_f1/*.tsv") + glob.glob(f"{S}/per_cohort_f1/*.json")
                + [f for f in glob.glob(f"{S}/de_agreement/*") if not f.endswith("manifest.json")]),
 "scripts": rec(glob.glob(f"{S}/scripts/*") + [f"{REPO}/bin/{b}" for b in
               ("build_local_pseudobulks.py", "aggregate_celltypes_gemma.py", "limma_voom_analysis.R", "aggregate_pairwise.py",
                "average_corr.py", "per_cohort_f1.py", "merged_author_de_agreement.py", "de_agreement_vs_annotation.py")]),
}
json.dump(m, open(f"{S}/manifest.json", "w"), indent=1)
print(f"{S}/manifest.json: {len(m['inputs'])} inputs, {len(m['outputs'])} outputs")
