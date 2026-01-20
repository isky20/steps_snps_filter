#!/usr/bin/env python3
import argparse, json, os
import pandas as pd

python3 annot_from_map.py   --in top100_traits_cardiology.tsv   --sep $'\t'   --trait-col trait   --map cardio.json   --out top100_traits_cardiology_annot.tsv



def load_map(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def main():
    ap = argparse.ArgumentParser(description="Annotate traits using a direct mapping dict (JSON).")
    ap.add_argument("--in", dest="in_path", required=True, help="Input TSV")
    ap.add_argument("--sep", default="\t", help="Separator (default: tab)")
    ap.add_argument("--trait-col", default="trait", help="Trait column name (default: trait)")
    ap.add_argument("--map", dest="map_path", required=True, help="JSON file: {trait: new_label}")
    ap.add_argument("--out", required=True, help="Output TSV")
    ap.add_argument("--fallback", default="keep", choices=["keep", "na"], help="If trait not found: keep original or NA")
    args = ap.parse_args()

    df = pd.read_csv(args.in_path, sep=args.sep, dtype=str)
    if args.trait_col not in df.columns:
        raise SystemExit(f"[ERROR] trait column not found: {args.trait_col}. Available: {list(df.columns)}")

    mp = load_map(args.map_path)

    # extract trait column as clean string
    df[args.trait_col] = df[args.trait_col].astype(str).str.strip()

    # build new annotation columns
    df["map_label"] = df[args.trait_col].map(mp)

    if args.fallback == "keep":
        df["map_label"] = df["map_label"].fillna(df[args.trait_col])
    else:
        df["map_label"] = df["map_label"].astype("string")

    # 4 annotated columns output (you can rename as you like)
    # trait, count, row_id, map_label
    keep_cols = []
    for c in [args.trait_col, "count", "row_id", "label"]:
        if c in df.columns:
            keep_cols.append(c)

    out_df = df[keep_cols].copy()
    out_df.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
