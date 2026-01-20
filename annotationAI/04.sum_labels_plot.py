#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

# -------------------------
# Helpers
#python3 sum_labels_plot.py \
# --indir out_group_traits \
# --outdir out_sum --top 20 \
# --title "Visualization Annotation Database"
# -------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default=".", help="folder with top*.csv")
    ap.add_argument("--outdir", default="out_sum", help="output folder")
    ap.add_argument("--top", type=int, default=20, help="top N labels")
    ap.add_argument("--title", default="Visualization Annotation Database", help="main title")
    args = ap.parse_args()

    indir = Path(args.indir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for f in sorted(indir.glob("top*.tsv")):
        df = pd.read_csv(f, sep='\t')

        # must contain columns: label, count
        df["label"] = df["label"].astype(str).str.strip()
        subtitle = f.stem.split("_")[-1].title()
        top = (df.groupby("label", as_index=False)["count"].sum()
                 .sort_values("count", ascending=False)
                 .head(args.top))

        # plot
        plt.figure(figsize=(14, max(4, 0.45 * len(top))))
        plt.barh(top["label"], top["count"])
        plt.gca().invert_yaxis()
        plt.xlabel("Count")
        plt.title(f"{args.title} — {subtitle}")
        plt.tight_layout()
        plt.savefig(outdir / f"{subtitle}_selected.png", dpi=200)
        plt.close()

        # save summed csv too
        top.to_csv(outdir / f"{subtitle}_selected.csv", index=False)

    print("Done.")


if __name__ == "__main__":
    main()
