#!/usr/bin/env python3
"""manifest.json for grid_labels_full_2026-11-09: grid code state, inputs and outputs (sha256), timestamps."""
import datetime, glob, hashlib, json, os, subprocess
REPO = "/space/grp/rschwartz/rschwartz/psychENCODE-reanalysis"
S = f"{REPO}/results_v2.1.2_limma/grid_labels_full_2026-11-09"
G = "/cosmos/data/aim2-annotation-grid-results/full_cells_2026-11-09"
GRID = "/space/grp/rschwartz/rschwartz/aim2-reannotation/nf-annotation-grid"
CACHE = "/cosmos/data/aim2-annotation-grid-cache"


def sha(f):
    h = hashlib.sha256()
    with open(f, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def entry(f):
    return {"sha256": sha(f), "mtime": datetime.datetime.fromtimestamp(os.path.getmtime(f)).isoformat(timespec="seconds")}


git = lambda d, *a: subprocess.check_output(["git", "-C", d, *a], text=True).strip()
prev = json.load(open(f"{REPO}/results_v2.1.2/manifest/manifest.json"))
m = {
    "created": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
    "data_version": "grid_labels_full_2026-11-09: two nf-annotation-grid scVI configs (cutoff 0, 100 ref cells per type) "
                    "on all cells of CMC, MultiomeBrain, SZBDMulti-Seq; QC-masked local pseudobulks; limma-voom subclass DE",
    "grid_code": {"repo": GRID, "commit": git(GRID, "rev-parse", "HEAD"),
                  "uncommitted": git(GRID, "status", "--short").splitlines(),
                  "scripts_used_sha256": {f: sha(f) for f in [f"{GRID}/bin/{s}" for s in
                                          ("annotate_query.py", "predict_scvi.py", "classify_query.py", "utils.py")]}},
    "reanalysis_repo": {"branch": git(REPO, "rev-parse", "--abbrev-ref", "HEAD"), "commit": git(REPO, "rev-parse", "HEAD")},
    "scvi_model": entry(f"{CACHE}/scvi_model/homo_sapiens/2024-07-01/scvi-homo_sapiens-2024-07-01/model.pt"),
    "references": {f: entry(f) for f in [f"{CACHE}/refs/homo_sapiens/2024-07-01/brain/dataset_id/sub_100/refs/{r}.h5ad"
                                         for r in ("Human_Multiple_Cortical_Areas_SMART-seq", "aggregated")]},
    "source_h5ad_sha256": {k: v for k, v in prev["input_h5ad_sha256_after"].items()
                           if any(c in k for c in ("CMC", "MultiomeBrain", "SZBDMulti-Seq"))},
    "grid_probs": {f: entry(f) for f in sorted(glob.glob(f"{G}/work/*/predict_*/probs/*.tsv.gz"))},
    "grid_cta": {f: entry(f) for f in sorted(glob.glob(f"{G}/cta/*/*.tsv"))},
    "cell_counts": {f: entry(f) for f in sorted(glob.glob(f"{S}/pseudobulks/*/ct_pseudobulks/cell_counts.tsv"))},
    "de_results_per_arm": {a: len(glob.glob(f"{S}/de/{a}/*/subclass/limma_voom/gemma/*/*/*_results.tsv"))
                           for a in sorted(os.listdir(f"{S}/de"))},
    "comparison": {f: entry(f) for f in sorted(glob.glob(f"{S}/comparison/*.tsv"))},
    "agreement_outputs": {f: entry(f) for f in sorted(glob.glob(f"{S}/de_agreement/*"))},
    "scripts": sorted(glob.glob(f"{S}/scripts/*")),
    "tmux_sessions": ["grid_full_cells", "grid_de", "grid_agreement"],
}
json.dump(m, open(f"{S}/manifest.json", "w"), indent=1)
print("wrote", f"{S}/manifest.json")
