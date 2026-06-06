#!/usr/bin/env python3
"""LOCO step 05 — plot MC and OC vs threshold for each (region, rank).

Produces:
    results/loco/evaluation/curves/<region>_<rank>.png
    results/loco/evaluation/curves/figure_oc_species.png   # all regions, OC at species
    results/loco/evaluation/curves/figure_mc_species.png

Reads results/loco/evaluation/mc_oc_curve.csv.
"""
import csv
from collections import defaultdict
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import _common  # noqa

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REGIONS = ["V4", "V3-V4", "V5-V7", "V1-V9", "RRN"]
RANKS   = ["class", "order", "family", "genus", "species"]
COLORS  = {
    "taxonova":    "#d62728",
    "idtaxa":      "#1f77b4",
    "naive-bayes": "#2ca02c",
    "blast":       "#9467bd",
    "vsearch":     "#ff7f0e",
    "spingo":      "#8c564b",
}


def load_curve(fp):
    rows = []
    with open(fp) as fh:
        for r in csv.DictReader(fh):
            r["threshold"] = float(r["threshold"])
            r["mc_rate"]   = float(r["mc_rate"])
            r["oc_rate"]   = float(r["oc_rate"])
            rows.append(r)
    return rows


def per_region_per_rank_plot(rows, out_dir):
    by = defaultdict(list)
    for r in rows:
        by[(r["region"], r["rank"], r["classifier"])].append(r)
    out_dir.mkdir(parents=True, exist_ok=True)
    for region in REGIONS:
        for rank in RANKS:
            fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharex=True)
            for cls, color in COLORS.items():
                pts = sorted(by.get((region, rank, cls), []), key=lambda x: x["threshold"])
                if not pts:
                    continue
                xs  = [p["threshold"] for p in pts]
                axes[0].plot(xs, [p["mc_rate"] for p in pts], "-", color=color, label=cls, lw=1.7)
                axes[1].plot(xs, [p["oc_rate"] for p in pts], "-", color=color, label=cls, lw=1.7)
            axes[0].set_title(f"{region} — {rank}: MC")
            axes[1].set_title(f"{region} — {rank}: OC")
            for ax in axes:
                ax.set_xlabel("confidence threshold τ")
                ax.set_ylabel("error rate")
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1.0)
                ax.grid(alpha=0.3)
            axes[1].legend(fontsize=8, loc="upper right")
            fig.tight_layout()
            fp = out_dir / f"{region}_{rank}.png"
            fig.savefig(fp, dpi=150)
            plt.close(fig)
            print(f"[loco/05] wrote {fp}")


def overview_figure(rows, out_dir):
    """5x2 grid: rows=regions, cols=(MC species, OC species)."""
    by = defaultdict(list)
    for r in rows:
        if r["rank"] != "species":
            continue
        by[(r["region"], r["classifier"])].append(r)
    fig, axes = plt.subplots(len(REGIONS), 2, figsize=(11, 14), sharex=True)
    for i, region in enumerate(REGIONS):
        for cls, color in COLORS.items():
            pts = sorted(by.get((region, cls), []), key=lambda x: x["threshold"])
            if not pts:
                continue
            xs  = [p["threshold"] for p in pts]
            axes[i, 0].plot(xs, [p["mc_rate"] for p in pts], "-", color=color, label=cls, lw=1.7)
            axes[i, 1].plot(xs, [p["oc_rate"] for p in pts], "-", color=color, label=cls, lw=1.7)
        axes[i, 0].set_ylabel(f"{region}\nrate")
        axes[i, 0].grid(alpha=0.3)
        axes[i, 1].grid(alpha=0.3)
        axes[i, 0].set_xlim(0, 1); axes[i, 1].set_xlim(0, 1)
        axes[i, 0].set_ylim(0, 1.0); axes[i, 1].set_ylim(0, 1.0)
    axes[0, 0].set_title("Mis-classification (species)")
    axes[0, 1].set_title("Over-classification (species)")
    for ax in axes[-1]:
        ax.set_xlabel("confidence threshold τ")
    axes[0, 1].legend(fontsize=9, loc="upper right")
    fig.tight_layout()
    fp = out_dir / "figure_species_mc_oc.png"
    fig.savefig(fp, dpi=200)
    plt.close(fig)
    print(f"[loco/05] wrote {fp}")


def main():
    cfg = _common.load_pipeline_config()
    eval_dir = Path(cfg["results_dir"]) / "loco" / "evaluation"
    curve_fp = eval_dir / "mc_oc_curve.csv"
    if not curve_fp.exists():
        sys.exit(f"missing {curve_fp}; run scripts/loco/04_evaluate.py first")
    rows = load_curve(curve_fp)
    out_dir = eval_dir / "curves"
    per_region_per_rank_plot(rows, out_dir)
    overview_figure(rows, out_dir)


if __name__ == "__main__":
    main()
