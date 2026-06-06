#!/usr/bin/env python3
"""Plot confidence-threshold sweep curves (Reviewer 2 comment 2.3).

Reads results/sweep/evaluation/all_metrics.csv (produced by
04b_evaluate_sweep.py) and writes one supplementary figure per (mock,
metric) showing precision / recall / f-score / detection_rate as a
function of confidence threshold, one line per classifier, faceted by
region and taxonomic level.

Output: results/sweep/evaluation/figures/<mock>_<metric>.pdf

The figure caption should note that confidence semantics differ across
tools (softmax posterior for taxonova / naive-bayes; bootstrap support
for SPINGO; consensus fraction for vsearch / blast; LCA k-mer fraction
for Kraken 2; per-rank score for IDTAXA / MAPseq), so curves are
informative within-tool, not directly comparable across tools at a
fixed numeric threshold.
"""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

CLASSIFIER_ORDER = [
    "taxonova", "idtaxa", "naive-bayes", "vsearch",
    "blast", "kraken2", "mapseq", "spingo",
]

CLASSIFIER_COLOR = {
    "taxonova":    "#d62728",
    "idtaxa":      "#1f77b4",
    "naive-bayes": "#2ca02c",
    "vsearch":     "#9467bd",
    "blast":       "#8c564b",
    "kraken2":     "#e377c2",
    "mapseq":      "#7f7f7f",
    "spingo":      "#bcbd22",
}

LEVELS = ["Family", "Genus", "Species"]
REGIONS = ["V4", "V3-V4", "V5-V7"]
METRIC_LABEL = {
    "precision": "Precision",
    "recall": "Recall",
    "f_score": "F-score",
    "detection_rate": "Taxon detection rate",
}


def _plot_metric(df: pd.DataFrame, mock: str, metric: str, out_fp: Path):
    fig, axes = plt.subplots(
        nrows=len(LEVELS), ncols=len(REGIONS),
        figsize=(4.0 * len(REGIONS), 2.8 * len(LEVELS)),
        sharex=True, sharey=True,
    )
    for r, level in enumerate(LEVELS):
        for c, region in enumerate(REGIONS):
            ax = axes[r][c] if len(LEVELS) > 1 else axes[c]
            sub = df[(df["mock"] == mock) & (df["region"] == region) &
                     (df["level"] == level)]
            if sub.empty:
                ax.set_visible(False)
                continue
            for cls in CLASSIFIER_ORDER:
                cls_sub = sub[sub["classifier"] == cls].sort_values("threshold")
                if cls_sub.empty:
                    continue
                ax.plot(
                    cls_sub["threshold"], cls_sub[metric],
                    marker="o", markersize=3, linewidth=1.2,
                    label=cls, color=CLASSIFIER_COLOR.get(cls, "k"),
                )
            ax.set_ylim(-0.02, 1.02)
            ax.set_xlim(-0.02, 1.02)
            ax.grid(True, linewidth=0.3, alpha=0.4)
            if r == 0:
                ax.set_title(region, fontsize=10)
            if c == 0:
                ax.set_ylabel(f"{level}\n{METRIC_LABEL[metric]}", fontsize=9)
            if r == len(LEVELS) - 1:
                ax.set_xlabel("Confidence threshold", fontsize=9)

    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center",
               ncol=min(len(labels), 8), fontsize=8,
               bbox_to_anchor=(0.5, 1.02), frameon=False)
    fig.suptitle(f"{mock} — {METRIC_LABEL[metric]} vs confidence threshold",
                 y=1.06, fontsize=11)
    fig.tight_layout()
    out_fp.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_fp, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out_fp}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metrics-csv", required=True,
                    help="results/sweep/evaluation/all_metrics.csv")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--metrics", nargs="+",
                    default=["precision", "recall", "f_score", "detection_rate"])
    args = ap.parse_args()

    df = pd.read_csv(args.metrics_csv)
    out_dir = Path(args.out_dir)
    for mock in sorted(df["mock"].unique()):
        for metric in args.metrics:
            _plot_metric(df, mock, metric, out_dir / f"{mock}_{metric}.pdf")


if __name__ == "__main__":
    main()
