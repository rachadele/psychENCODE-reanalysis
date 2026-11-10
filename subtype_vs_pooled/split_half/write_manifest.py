"""manifest.json: sha256 and mtimes of inputs, scripts and outputs, seeds, data version, timestamps."""
import datetime
import glob
import hashlib
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import CONTRASTS, OUT, PB, REPO, SEEDS, WORK, de_path, jobs  # noqa: E402


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def entry(p):
    return {"sha256": sha(p), "mtime": datetime.datetime.fromtimestamp(os.path.getmtime(p)).astimezone().isoformat()}


def main():
    rel = lambda p: os.path.relpath(p, REPO)
    inputs = [f"{REPO}/bin/limma_voom_analysis.R", f"{REPO}/nextflow.config",
              f"{REPO}/subtype_vs_pooled/cells_per_donor.tsv"] + sorted(glob.glob(f"{REPO}/gemma/metadata/*.tsv"))
    for _, sch, ct in jobs():
        inputs += [f"{PB[sch]}/{ct}/{ct}_pseudobulk_matrix.tsv.gz", f"{PB[sch]}/{ct}/{ct}_gene_ncells.tsv.gz",
                   f"{PB[sch]}/cell_counts.tsv"] + [de_path(sch, ct, c) for c in CONTRASTS]
    inputs = sorted(set(inputs))
    outs = sorted(glob.glob(f"{OUT}/*.tsv"))
    scripts = sorted(glob.glob(f"{OUT}/*.py") + glob.glob(f"{OUT}/*.sh"))
    timing = open(f"{WORK}/logs/timing.log").read().splitlines() if os.path.exists(f"{WORK}/logs/timing.log") else []
    m = {"created": datetime.datetime.now().astimezone().isoformat(),
         "data_version": "pipeline labels sc-annotation-pipeline 2.1.2; DE results_v2.1.2_limma (author, ours) and "
                         "results_v2.1.2_limma/merged_author_2026-11-09 (merged author); limma-voom, "
                         "~ Disorder + Age_death + PMI + Biological_Sex + X1000G_ancestry, no cohort term",
         "repo": {"branch": subprocess.run(["git", "-C", REPO, "branch", "--show-current"], capture_output=True, text=True).stdout.strip(),
                  "commit": subprocess.run(["git", "-C", REPO, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()},
         "seeds": SEEDS, "rng": "numpy.random.default_rng(seed)",
         "stratification": "cohort x Disorder x Sex strata, donors sorted by age and paired, one of each pair per half",
         "limma_params": {"min_cells_per_sample": 0, "min_cells_per_gene": 3, "min_variance": 0.01},
         "tmux_session": "split_half_de", "fit_store": f"{WORK}/fits", "timing": timing,
         "check_full_data_reproduction": "fit_half.py with all donors reproduced existing t exactly (max abs diff 0) "
                                         "for author L6_IT_Car3 and merged Lamp5, all eight contrasts",
         "inputs": {rel(p): entry(p) for p in inputs},
         "scripts": {rel(p): entry(p) for p in scripts},
         "outputs": {rel(p): entry(p) for p in outs}}
    json.dump(m, open(f"{OUT}/manifest.json", "w"), indent=1)


if __name__ == "__main__":
    main()
