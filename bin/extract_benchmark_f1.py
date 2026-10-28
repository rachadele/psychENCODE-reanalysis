#!/usr/bin/env python3
"""Extract per-label benchmark F1 for one classifier configuration.

Reads label_results_summary.tsv.gz from an evaluation_summary.nf run and keeps the
rows for one method / reference / reference size / cutoff, at every taxonomy level.
Writes a TSV (key, label, mean_f1_score, sd_f1_score, n_samples) that
average_corr.py --f1_path reads, plus a JSON provenance sidecar (sha256 of input and
output, data_version, timestamps).
"""
import argparse
import datetime
import hashlib
import json
import os

import pandas as pd


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--label_summary", required=True, help="label_results_summary.tsv.gz")
    p.add_argument("--run_manifest", default=None, help="manifest.json of the evaluation_summary.nf run")
    p.add_argument("--method", default="scvi_rf")
    p.add_argument("--reference", default="whole cortex")
    p.add_argument("--subsample_ref", default="500")
    p.add_argument("--cutoff", type=float, default=0.0)
    p.add_argument("--out", required=True)
    a = p.parse_args()

    df = pd.read_csv(a.label_summary, sep="\t")
    sel = df[(df["method"] == a.method)
             & (df["reference"] == a.reference)
             & (df["subsample_ref"].astype(str) == str(a.subsample_ref))
             & (df["cutoff"].astype(float) == a.cutoff)]
    if sel.empty:
        raise SystemExit("No rows match the requested configuration")
    out = sel.rename(columns={"f1_score_mean": "mean_f1_score", "f1_score_std": "sd_f1_score",
                              "count": "n_samples"})
    out = out[["key", "label", "mean_f1_score", "sd_f1_score", "n_samples"]].sort_values(["key", "label"])
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    out.to_csv(a.out, sep="\t", index=False)

    run = {}
    if a.run_manifest and os.path.exists(a.run_manifest):
        with open(a.run_manifest) as f:
            run = json.load(f)
    side = {
        "output": os.path.abspath(a.out),
        "output_sha256": sha256(a.out),
        "input": os.path.abspath(a.label_summary),
        "input_sha256": sha256(a.label_summary),
        "input_mtime": datetime.datetime.fromtimestamp(os.path.getmtime(a.label_summary)).astimezone().isoformat(),
        "data_version": {
            "census_version": run.get("params", {}).get("census_version"),
            "evaluation_summary_commit": run.get("git_commit"),
            "evaluation_summary_run": run.get("timestamp"),
            "removed_studies": run.get("params", {}).get("remove_outliers"),
        },
        "filter": {"method": a.method, "reference": a.reference,
                   "subsample_ref": a.subsample_ref, "cutoff": a.cutoff},
        "created": datetime.datetime.now().astimezone().isoformat(),
    }
    with open(os.path.splitext(a.out)[0] + ".json", "w") as f:
        json.dump(side, f, indent=2)
    print(f"Wrote {len(out)} rows to {a.out}")


if __name__ == "__main__":
    main()
