#!/usr/bin/env python3
"""List all MGnify biomes with sample counts and write them to biomes.csv.

Run this once (or whenever the MGnify catalogue changes) so `fetch_top_taxa.py`
has an up-to-date lineage table to look up biome short-names.

Output columns: Environment, API_Lineage, Sample_Count, Category.
"""

import argparse
import time
from pathlib import Path

import pandas as pd
import requests


MGNIFY_BIOMES_URL = "https://www.ebi.ac.uk/metagenomics/api/v1/biomes"
DEFAULT_OUT = Path(__file__).parent / "biomes.csv"


def categorize(lineage: str) -> str:
    if "Environmental:Aquatic:Marine" in lineage:
        return "Marine"
    if "Environmental:Aquatic" in lineage:
        return "Freshwater"
    if "Environmental:Terrestrial" in lineage:
        return "Terrestrial"
    if "Host-associated:Human" in lineage:
        return "Human-associated"
    if "Host-associated" in lineage:
        return "Host-associated"
    if "Engineered" in lineage:
        return "Engineered"
    return "Other"


def fetch_all_biomes(session: requests.Session) -> list:
    biomes, page = [], 1
    while True:
        r = session.get(MGNIFY_BIOMES_URL,
                        params={"page": page, "page_size": 100}, timeout=30)
        r.raise_for_status()
        data = r.json().get("data", [])
        if not data:
            break
        biomes.extend(data)
        print(f"  page {page}: +{len(data)} biomes (total {len(biomes)})")
        if not r.json().get("links", {}).get("next"):
            break
        page += 1
        time.sleep(0.1)
    return biomes


def to_dataframe(biomes: list) -> pd.DataFrame:
    rows = []
    for b in biomes:
        attrs = b["attributes"]
        if attrs["samples-count"] <= 0:
            continue
        rows.append({
            "Environment":  attrs["biome-name"],
            "API_Lineage":  attrs["lineage"],
            "Sample_Count": attrs["samples-count"],
            "Category":     categorize(attrs["lineage"]),
        })
    df = pd.DataFrame(rows).sort_values("Sample_Count", ascending=False)
    return df.reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(DEFAULT_OUT),
                    help="Output CSV (default: ./biomes.csv)")
    ap.add_argument("--top", type=int, default=20,
                    help="Print this many top biomes to stdout (default: 20)")
    args = ap.parse_args()

    print(f"Fetching biomes from {MGNIFY_BIOMES_URL} ...")
    with requests.Session() as s:
        s.headers["Accept"] = "application/json"
        biomes = fetch_all_biomes(s)

    df = to_dataframe(biomes)
    df.to_csv(args.out, index=False)
    print(f"\nWrote {len(df)} biomes -> {args.out}")

    print(f"\nTop {args.top} by sample count:")
    for _, r in df.head(args.top).iterrows():
        print(f"  {r['Sample_Count']:>7,}  {r['Environment']:<30}  {r['API_Lineage']}")


if __name__ == "__main__":
    main()
