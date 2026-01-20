#!/usr/bin/env python3
import argparse
import json
import os
import re
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# -------------------------
# Helpers
#python3 plot_top_traits_by_group.py   --long out_traits/traits_long_with_groups.tsv   --sep $'\t' \
# --rules rules.json   --outdir out_group_traits   --top 20 \
#   --main-title "Visualization Annotation Database"   --unique-variants   --pdf
# -------------------------



def load_rules(dict_path: str) -> dict[str, str]:
    # We only need group names (keys) here to control the order
    p = Path(dict_path)
    if p.suffix.lower() == ".json":
        with open(p, "r", encoding="utf-8") as f:
            rules = json.load(f)
        if not isinstance(rules, dict):
            raise ValueError("Rules JSON must be a dict: {group_name: regex}")
        return rules

    # TXT format: group<TAB>regex (ignore blank and # comments)
    rules = {}
    with open(p, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t", 1)
            if len(parts) != 2:
                raise ValueError(f"Bad rules line (need group<TAB>regex): {line}")
            k, v = parts[0].strip(), parts[1].strip()
            rules[k] = v
    return rules

def short_label(trait: str, row_id) -> str:
    words = str(trait).split()
    if len(words) > 5:
        rid = "" if pd.isna(row_id) else str(int(row_id))
        return f"{rid}: " + " ".join(words[:3])
    return str(trait)

def safe_name(s: str) -> str:
    s = s.strip().lower()
    s = re.sub(r"[^a-z0-9]+", "_", s)
    return s.strip("_") or "group"


def main():
    ap = argparse.ArgumentParser(description="Per-group TOP traits plots from traits_long_with_groups table.")
    ap.add_argument("--long", required=True, help="traits_long_with_groups.tsv(.gz) produced by annot_groups.py")
    ap.add_argument("--sep", default="\t", help="Separator for --long (default: tab)")
    ap.add_argument("--rules", required=True, help="Rules file used (json or txt). Keys define group order.")
    ap.add_argument("--outdir", default="out_group_traits", help="Output directory")
    ap.add_argument("--top", type=int, default=20, help="Top N traits per group")
    ap.add_argument("--main-title", default="Visualization Annotation Database", help="Main title on plots")
    ap.add_argument("--unique-variants", action="store_true",
                    help="Count unique CHR_POS_REF_ALT per trait (recommended). If missing, falls back to row-count.")
    ap.add_argument("--pdf", action="store_true", help="Also write a single multipage PDF with all group plots")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.long, sep=args.sep)
    needed = {"trait", "group"}
    missing = sorted(list(needed - set(df.columns)))
    if missing:
        raise SystemExit(f"[ERROR] Missing required columns in --long: {missing}. Found: {list(df.columns)}")

    rules = load_rules(args.rules)
    group_order = [g for g in rules.keys() if g in set(df["group"].dropna().unique())]

    # If some groups exist in data but not in rules, append them at the end
    extra = [g for g in df["group"].dropna().unique() if g not in rules]
    group_order.extend(sorted(extra))

    pdf_path = Path(args.outdir) / "top_traits_by_group.pdf"
    pdf = PdfPages(pdf_path) if args.pdf else None

    for g in group_order:
        sub = df[df["group"] == g].copy()
        if sub.empty:
            continue

        # Count
        if args.unique_variants and "CHR_POS_REF_ALT" in sub.columns:
            counts = sub.groupby("trait")["CHR_POS_REF_ALT"].nunique()
            count_label = "Unique variants"
        else:
            counts = sub.groupby("trait").size()
            count_label = "Count"
        
        row_id_min = sub.groupby("trait")["row_id"].min() if "row_id" in sub.columns else None

        top_df = (
            counts.sort_values(ascending=False)
                .head(args.top)
                .reset_index(name="count")
        )

        if row_id_min is not None:
            top_df["row_id"] = top_df["trait"].map(row_id_min)
        else:
            top_df["row_id"] = pd.NA

        top_df["label"] = [short_label(t, r) for t, r in zip(top_df["trait"], top_df["row_id"])]

        # Save table
        out_table = Path(args.outdir) / f"top{args.top}_traits_{safe_name(g)}.tsv"
        top_df.to_csv(out_table, sep="\t", index=False)

        # Plot (horizontal bar)
        h = max(4.5, 0.45 * len(top_df))
        fig, ax = plt.subplots(figsize=(13, h))

        ax.barh(top_df["label"], top_df["count"])
        ax.invert_yaxis()
        ax.set_xlabel(count_label)

        # Title = main + subgroup
        fig.suptitle(args.main_title, fontsize=16, y=0.98)
        ax.set_title(f"{g} — Top {args.top} traits", fontsize=13, pad=10)

        # Add numbers on bars
        for i, v in enumerate(top_df["count"]):
            ax.text(v, i, f" {v}", va="center")

        # ---- Centering / layout ----
        max_lab_len = int(top_df["label"].astype(str).str.len().max())
        left = min(0.38, max(0.18, 0.12 + 0.006 * max_lab_len))   # dynamic left margin
        fig.subplots_adjust(left=left, right=0.96, top=0.88, bottom=0.10)

        out_png = Path(args.outdir) / f"top{args.top}_traits_{safe_name(g)}.png"
        fig.savefig(out_png, dpi=200)
        plt.close(fig)

    if pdf is not None:
        pdf.close()
        print(f"[OK] PDF written: {pdf_path}")

    print(f"[OK] Done. Outputs in: {args.outdir}")


if __name__ == "__main__":
    main()
