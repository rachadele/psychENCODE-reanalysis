import argparse
import pandas as pd


def parse_arguments():
    parser = argparse.ArgumentParser(description="Aggregate per-query F1 scores to per-label means for heatmap annotation")
    parser.add_argument("--input", required=True, help="Path to label_results.tsv.gz")
    parser.add_argument("--output", required=True, help="Output path for label_metrics_stats_per_label.tsv")
    parser.add_argument("--method", default="scvi", help="Annotation method to filter on")
    parser.add_argument("--cutoff", type=float, default=0.05, help="Confidence cutoff to filter on")
    parser.add_argument("--subsample_ref", type=int, default=100, help="Reference subsample size to filter on")
    parser.add_argument("--keys", nargs="+", default=["class", "subclass"], help="Annotation levels to include")
    return parser.parse_args()


def main():
    args = parse_arguments()

    df = pd.read_csv(args.input, sep="\t")
    df = df[
        (df["method"] == args.method) &
        (df["cutoff"] == args.cutoff) &
        (df["subsample_ref"] == args.subsample_ref) &
        (df["key"].isin(args.keys))
    ]

    result = (
        df.groupby(["key", "label"])["f1_score"]
        .mean()
        .reset_index()
        .rename(columns={"f1_score": "mean"})
    )
    result["metric"] = "f1_score"
    result[["key", "metric", "label", "mean"]].to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(result)} rows to {args.output}")


if __name__ == "__main__":
    main()
