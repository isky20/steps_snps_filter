#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -------------------------
# Helpers
# python3 annot_groups.py   --in annotated_paragenomicarray_mafmore1.csv   --sep "," \
# --col nhgri_traits   --dict rules.json   --keep-cols CHR_POS_REF_ALT,RSID \
# --outdir out_traits   --top 20   --orient h   --title "Visualization Annotation Database NHGRI MAF>0.01"
# -------------------------
def cell_to_list(v, split_char="|"):
    """Turn a cell into a cleaned list of trait strings."""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return []
    if isinstance(v, np.ndarray):
        v = v.tolist()
    if not isinstance(v, list):
        v = [v]

    out = []
    for item in v:
        for part in str(item).split(split_char):
            part = part.strip()
            if part and part != "---" and part.lower() not in {"nan", "none"}:
                out.append(part)
    return out


def read_table(path: str, sep: str, col: str) -> pd.DataFrame:
    p = str(path)
    if p.endswith(".parquet"):
        df = pd.read_parquet(p)
    else:
        df = pd.read_csv(p, sep=sep, compression="infer", low_memory=False)
    if col not in df.columns:
        raise SystemExit(f"[ERROR] Column not found: {col}. Available: {list(df.columns)[:30]} ...")
    return df


def load_rules(dict_path: str, dict_format: str | None = None) -> dict[str, str]:
    """
    Load group->regex rules from JSON or TXT.

    JSON format:
      { "Group A": "regex", "Group B": "regex" }

    TXT format (one rule per line):
      Group A<TAB>regex
      Group B: regex
    Lines starting with # are ignored.
    """
    dp = Path(dict_path)
    if dict_format is None:
        ext = dp.suffix.lower()
        dict_format = "json" if ext == ".json" else "txt"

    if dict_format == "json":
        with open(dp, "r", encoding="utf-8") as f:
            rules = json.load(f)
        if not isinstance(rules, dict):
            raise SystemExit("[ERROR] JSON must be an object: {group: regex, ...}")
        # ensure str->str
        rules = {str(k): str(v) for k, v in rules.items()}
        return rules

    if dict_format == "txt":
        rules = {}
        with open(dp, "r", encoding="utf-8") as f:
            for ln in f:
                ln = ln.strip()
                if not ln or ln.startswith("#"):
                    continue
                if "\t" in ln:
                    k, v = ln.split("\t", 1)
                elif ":" in ln:
                    k, v = ln.split(":", 1)
                else:
                    raise SystemExit(f"[ERROR] Bad rule line (need TAB or ':'): {ln}")
                k, v = k.strip(), v.strip()
                if not k or not v:
                    continue
                rules[k] = v
        return rules

    raise SystemExit(f"[ERROR] Unsupported dict_format: {dict_format} (use json or txt)")


def compile_rules(rules: dict[str, str], ignore_case: bool = True) -> list[tuple[str, re.Pattern]]:
    flags = re.IGNORECASE if ignore_case else 0
    compiled = []
    for name, pat in rules.items():
        try:
            compiled.append((name, re.compile(pat, flags=flags)))
        except re.error as e:
            raise SystemExit(f"[ERROR] Regex compile failed for '{name}': {e}\nPattern: {pat}")
    return compiled


def assign_groups_to_trait(trait: str, compiled: list[tuple[str, re.Pattern]], first_match: bool) -> list[str]:
    hits = []
    for gname, gpat in compiled:
        if gpat.search(trait):
            hits.append(gname)
            if first_match:
                break
    return hits


# -------------------------
# Main
# -------------------------
def main():
    ap = argparse.ArgumentParser(description="Group NHGRI traits using regex rules (JSON/TXT), output counts + plot.")
    ap.add_argument("--in", dest="in_path", required=True, help="Input table (csv/tsv, optionally .gz; or .parquet)")
    ap.add_argument("--col", default="nhgri_traits", help="Column containing trait strings (default: nhgri_traits)")
    ap.add_argument("--sep", default="\t", help=r"Input separator (default: \t). Use ',' for CSV.")
    ap.add_argument("--dict", dest="dict_path", required=True, help="Rules file path (json or txt)")
    ap.add_argument("--dict-format", choices=["json", "txt"], default=None, help="Force rules format (auto if omitted)")
    ap.add_argument("--outdir", default="out_groups", help="Output directory")
    ap.add_argument("--split-char", default="|", help="Split character inside each cell (default: '|')")
    ap.add_argument("--top", type=int, default=20, help="Top N groups to plot (default: 20)")
    ap.add_argument("--first-match", action="store_true", help="Assign only the first matching group per trait")
    ap.add_argument("--percent-of", choices=["none", "traits", "rows"], default="none",
                    help="Show % instead of counts: traits=all assignments; rows=unique row hits")
    ap.add_argument("--title", default=None, help="Plot title (optional)")
    ap.add_argument("--orient", choices=["h", "v"], default="h", help="Plot orientation: h=horizontal, v=vertical")
    ap.add_argument("--keep-unmatched", action="store_true", help="Keep 'Unmatched' group for non-matching traits")
    ap.add_argument( "--keep-cols", default="", help="Comma-separated extra columns to keep in long output (e.g. CHR_POS_REF_ALT,RSID)"
)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load data + rules
    df = read_table(args.in_path, sep=args.sep, col=args.col)
    rules = load_rules(args.dict_path, dict_format=args.dict_format)
    compiled = compile_rules(rules, ignore_case=True)

    # Build long table: row_id, trait
    keep_cols = [c.strip() for c in args.keep_cols.split(",") if c.strip()]
    missing = [c for c in keep_cols if c not in df.columns]
    if missing:
        print(f"[WARN] keep-cols not found in input and will be skipped: {missing}")
    keep_cols = [c for c in keep_cols if c in df.columns]

    tmp = df[keep_cols + [args.col]].copy()
    tmp["row_id"] = df.index

    tmp["trait_list"] = tmp[args.col].apply(lambda x: cell_to_list(x, split_char=args.split_char))

    long_df = (
        tmp[keep_cols + ["row_id", "trait_list"]]
        .explode("trait_list")
        .rename(columns={"trait_list": "trait"})
    )
    long_df = long_df.dropna()
    long_df["trait"] = long_df["trait"].astype(str).str.strip()
    long_df = long_df[(long_df["trait"] != "") & (long_df["trait"] != "---") & (long_df["trait"].str.lower() != "nan")]

    # Assign groups
    groups_col = []
    for t in long_df["trait"].tolist():
        hits = assign_groups_to_trait(t, compiled, first_match=args.first_match)
        if not hits and args.keep_unmatched:
            hits = ["Unmatched"]
        groups_col.append(hits)

    long_df["group"] = groups_col
    grouped = long_df.explode("group")
    grouped = grouped.dropna(subset=["group"])

    # Count
    if args.percent_of == "rows":
        counts = grouped.groupby("group")["row_id"].nunique().sort_values(ascending=False)
        denom = df.shape[0]
        values = counts / denom * 100.0
        xlab = "% of rows"
    elif args.percent_of == "traits":
        counts = grouped["group"].value_counts()
        denom = counts.sum()
        values = counts / denom * 100.0
        xlab = "% of assignments"
    else:
        counts = grouped["group"].value_counts()
        values = counts
        xlab = "Count"

    # Save counts
    out_counts = pd.DataFrame({"group": values.index, "value": values.values})
    out_counts_path = outdir / "group_counts.tsv"
    out_counts.to_csv(out_counts_path, sep="\t", index=False)

    # Plot top N
    topn = out_counts.head(args.top).copy()
    title = args.title or f"Top {args.top} groups ({Path(args.in_path).name})"

    plt.figure(figsize=(10, max(4, 0.35 * len(topn))))
    if args.orient == "h":
        plt.barh(topn["group"][::-1], topn["value"][::-1])
        plt.xlabel(xlab)
        title = args.title 
        plt.title(title)
    else:
        plt.bar(topn["group"], topn["value"])
        plt.ylabel(xlab)
        title = args.title 
        plt.title(title)
        plt.xticks(rotation=45, ha="right")

    plt.tight_layout()
    plot_path = outdir / "group_counts_top.png"
    plt.savefig(plot_path, dpi=200)
    plt.close()

    # Optional: save long annotated traits
    long_path = outdir / "traits_long_with_groups.tsv.gz"
    grouped.to_csv(long_path, sep="\t", index=False, compression="gzip")

    print(f"[OK] Wrote: {out_counts_path}")
    print(f"[OK] Wrote: {plot_path}")
    print(f"[OK] Wrote: {long_path}")


if __name__ == "__main__":
    main()

