#!/usr/bin/env python3
"""Aggregate V5-V7 species-level Recall, TP, FP per (mock × method) for Table 1.

Reads each results/<mock>/V5-V7/<method>/<param>/taxonomy_final.tsv and
data/mock_community/<mock>/V5-V7/trueish_taxonomy_<mock>.tsv, computes
species-level overlap and writes a markdown table to stdout.
"""
import argparse, glob, os, sys
from pathlib import Path
import pandas as pd

def load_tax(fp):
    df = pd.read_csv(fp, sep="\t", comment=None)
    # try common id column names
    id_col = None
    for c in ("Feature ID", "feature_id", "id", "sequence_id", "ID", "FeatureID", "#OTU ID"):
        if c in df.columns:
            id_col = c; break
    if id_col is None:
        id_col = df.columns[0]
    tax_col = None
    for c in ("Taxon", "taxonomy", "Taxonomy", "taxon"):
        if c in df.columns:
            tax_col = c; break
    if tax_col is None:
        tax_col = df.columns[1]
    return dict(zip(df[id_col].astype(str), df[tax_col].astype(str)))

def species_of(t):
    if not t or t == "nan": return None
    parts = [p.strip() for p in t.split(";")]
    sp = None
    for p in parts:
        if p.startswith("s__"):
            v = p[3:].strip()
            if v: sp = v
    return sp

def evaluate(true_fp, pred_fp):
    t = load_tax(true_fp)
    p = load_tax(pred_fp)
    true_species = {species_of(v) for v in t.values()} - {None}
    pred_species = set()
    for sid, ptax in p.items():
        sp = species_of(ptax)
        if sp: pred_species.add(sp)
    tp = len(true_species & pred_species)
    fp = len(pred_species - true_species)
    fn = len(true_species - pred_species)
    rec = tp / (tp + fn) if (tp + fn) else 0.0
    return rec, tp, fp, fn, len(true_species), len(pred_species)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-root", default="/path/to/taxobench/results")
    ap.add_argument("--data-root", default="/path/to/taxobench/data")
    ap.add_argument("--region", default="V5-V7")
    ap.add_argument("--mocks", nargs="+",
                    default=["biome_digestive_system","biome_lung","biome_mixed",
                             "mockrobiota-16S","biome_skin","biome_soil","biome_sputum"])
    ap.add_argument("--methods", nargs="+",
                    default=["taxonova","naive-bayes","idtaxa","mapseq","spingo","kraken2","blast","vsearch"])
    args = ap.parse_args()

    rows = []
    for m in args.mocks:
        true_fp = (Path(args.data_root) / "mock_community" / m / args.region /
                   f"trueish_taxonomy_{m}.tsv")
        if not true_fp.exists():
            true_fp = next(iter(Path(args.data_root, "mock_community", m, args.region).glob("trueish*.tsv")), None)
        for meth in args.methods:
            method_dir = Path(args.results_root) / m / args.region / meth
            if not method_dir.is_dir(): continue
            for param_dir in sorted(method_dir.iterdir()):
                pred_fp = param_dir / "taxonomy_final.tsv"
                if not pred_fp.exists(): continue
                if not true_fp or not Path(true_fp).exists(): continue
                rec, tp, fp, fn, n_true, n_pred = evaluate(str(true_fp), str(pred_fp))
                rows.append(dict(mock=m, method=meth, param=param_dir.name,
                                 recall=rec, TP=tp, FP=fp, FN=fn,
                                 n_true=n_true, n_pred=n_pred))

    df = pd.DataFrame(rows)
    if df.empty:
        print("No rows. Run the pipeline first.", file=sys.stderr); return

    df.to_csv(f"v5v7_table1.csv", index=False)
    # markdown table
    pretty = {"taxonova":"Taxonova","naive-bayes":"QIIME 2","idtaxa":"IDTAXA",
              "mapseq":"MAPseq","spingo":"SPINGO","kraken2":"Kraken 2",
              "blast":"BLASTn","vsearch":"VSEARCH"}
    pivoted = {}
    for mock in args.mocks:
        sub = df[df.mock == mock].drop_duplicates(subset=["method"], keep="first")
        if sub.empty: continue
        pivoted[mock] = sub.set_index("method").to_dict("index")

    METHODS = ["taxonova","naive-bayes","idtaxa","mapseq","spingo","kraken2","blast","vsearch"]
    print(f"### {args.region}\n")
    print(f"| Mock Community | Metric | " + " | ".join(pretty[m] for m in METHODS) + " |")
    print("|---|---|" + "|".join(["---"]*len(METHODS)) + "|")
    for mock in args.mocks:
        if mock not in pivoted: continue
        for metric, key in [("Recall","recall"), ("TP","TP"), ("FP","FP")]:
            label = mock.replace("biome_","").replace("mockrobiota-16S","Mockrobiota").replace("digestive_system","Digestive").replace("_"," ").title()
            label = label if metric == "Recall" else ""
            row = [pretty[m] + ":" + str(round(pivoted[mock].get(m,{}).get(key,0), 3)) for m in METHODS]
            vals = [pivoted[mock].get(m,{}).get(key, "—") for m in METHODS]
            cells = []
            for v in vals:
                if isinstance(v, float):
                    cells.append(f"{v:.3f}" if metric == "Recall" else f"{v:.1f}")
                else:
                    cells.append(str(v))
            print(f"| {label} | {metric} | " + " | ".join(cells) + " |")

if __name__ == "__main__":
    main()
