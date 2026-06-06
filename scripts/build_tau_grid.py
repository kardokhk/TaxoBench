#!/usr/bin/env python3
"""Build per-classifier per-experiment τ-grids from empirical confidence distributions.

For each classifier, collect the per-prediction species-rank confidence across the
relevant eval set, then emit percentile-based τ-grids so each curve actually shows
its dynamic range. Writes configs/tau_grids.json keyed by classifier and experiment.

Continuous classifiers (taxonova, idtaxa, naive-bayes) get a fine percentile grid;
bimodal classifiers (mapseq, spingo) get a coarser grid sufficient for step curves;
degenerate-at-species classifiers (blast, vsearch, kraken2) are excluded.

Run after sweep and LOCO classification steps have populated their results trees.
"""
import argparse
import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import _common  # noqa


CONTINUOUS = ("taxonova", "idtaxa", "naive-bayes")
BIMODAL    = ("mapseq", "spingo")
EXCLUDED   = ("blast", "vsearch", "kraken2")

PERCENTILES_CONTINUOUS = [10, 25, 40, 50, 60, 70, 80, 85, 90, 93, 95, 97, 99]
PERCENTILES_BIMODAL    = [10, 25, 50, 75, 90, 99]

DEFAULT_TAU = {
    "taxonova":    0.7,
    "idtaxa":      0.6,
    "naive-bayes": 0.7,
    "mapseq":      0.4,
    "spingo":      0.5,
}


def _percentile(sorted_vals, p):
    if not sorted_vals:
        return 0.0
    if len(sorted_vals) == 1:
        return float(sorted_vals[0])
    k = (len(sorted_vals) - 1) * (p / 100.0)
    lo = int(math.floor(k))
    hi = int(math.ceil(k))
    if lo == hi:
        return float(sorted_vals[lo])
    return float(sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * (k - lo))


def _round_grid(values):
    seen = []
    for v in values:
        r = round(float(v), 4)
        r = max(0.0, min(1.0, r))
        if r not in seen:
            seen.append(r)
    return sorted(seen)


def _build_grid(species_confs, classifier, kind):
    """species_confs: list[float]. kind: 'continuous' | 'bimodal'."""
    pcts = PERCENTILES_CONTINUOUS if kind == "continuous" else PERCENTILES_BIMODAL
    sorted_vals = sorted(species_confs)
    grid = [_percentile(sorted_vals, p) for p in pcts]
    if classifier in DEFAULT_TAU:
        grid.append(DEFAULT_TAU[classifier])
    grid.extend([0.0, 1.0])
    return _round_grid(grid)


def _parse_semi(s, n=7):
    """Per-rank confidence parse with single-scalar broadcast (matches
    scripts/loco/04_evaluate.py::_parse_confidence_field)."""
    if not s:
        return [0.0] * n
    s = str(s).strip()
    if ";" in s:
        out = []
        for tok in s.split(";"):
            try:
                out.append(float(tok.strip()))
            except ValueError:
                out.append(0.0)
        if len(out) < n:
            out += [0.0] * (n - len(out))
        return out[:n]
    try:
        v = float(s)
        return [v] * n
    except ValueError:
        return [0.0] * n


def collect_sweep_species(sweep_root: Path, classifier: str):
    """Walk results/sweep/<mock>/<region>/<classifier>/<param>/taxonomy_final.tsv.
    Return list[float] of species-rank confidences."""
    out = []
    if not sweep_root.exists():
        return out
    for mock_dir in sorted(sweep_root.iterdir()):
        if not mock_dir.is_dir() or mock_dir.name == "evaluation":
            continue
        for region_dir in sorted(mock_dir.iterdir()):
            cls_dir = region_dir / classifier
            if not cls_dir.is_dir():
                continue
            for param_dir in sorted(cls_dir.iterdir()):
                fp = param_dir / "taxonomy_final.tsv"
                if not fp.exists():
                    continue
                with open(fp) as fh:
                    fh.readline()
                    for line in fh:
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) >= 3:
                            confs = _parse_semi(parts[2], n=7)
                            out.append(confs[6])
    return out


def collect_loco_species(loco_root: Path, classifier: str):
    """Walk results/loco/<region>/<classifier>/{novel,control}_predictions.tsv.
    Return list[float] of species-rank confidences."""
    out = []
    if not loco_root.exists():
        return out
    for region_dir in sorted(loco_root.iterdir()):
        if region_dir.name == "evaluation" or not region_dir.is_dir():
            continue
        cls_dir = region_dir / classifier
        if not cls_dir.is_dir():
            continue
        for fname in ("control_predictions.tsv", "novel_predictions.tsv"):
            fp = cls_dir / fname
            if not fp.exists():
                continue
            with open(fp) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                cols = reader.fieldnames or []
                conf_col = next((c for c in cols
                                 if c.lower() in ("confidence", "consensus")), None)
                if not conf_col:
                    continue
                # idtaxa lineage has 8 elements (Root + 7); drop Root for species index
                idtaxa = (classifier == "idtaxa")
                n = 8 if idtaxa else 7
                sp_idx = 7 if idtaxa else 6
                for row in reader:
                    confs = _parse_semi(row.get(conf_col, ""), n=n)
                    out.append(confs[sp_idx])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(ROOT / "configs" / "tau_grids.json"))
    ap.add_argument("--sweep-root", default=None)
    ap.add_argument("--loco-root", default=None)
    args = ap.parse_args()

    cfg = _common.load_pipeline_config()
    results = Path(cfg["results_dir"])
    sweep_root = Path(args.sweep_root) if args.sweep_root else results / "sweep"
    loco_root  = Path(args.loco_root)  if args.loco_root  else results / "loco"

    grids = {}
    distributions = {}

    for cls in CONTINUOUS + BIMODAL:
        kind = "continuous" if cls in CONTINUOUS else "bimodal"
        sweep_confs = collect_sweep_species(sweep_root, cls)
        loco_confs  = collect_loco_species(loco_root, cls)
        sweep_grid = _build_grid(sweep_confs, cls, kind) if sweep_confs else []
        loco_grid  = _build_grid(loco_confs,  cls, kind) if loco_confs  else []
        grids[cls] = {"sweep": sweep_grid, "loco": loco_grid, "kind": kind,
                      "default_tau": DEFAULT_TAU.get(cls)}
        distributions[cls] = {
            "sweep": {"n": len(sweep_confs), "p50": _percentile(sorted(sweep_confs), 50)
                                                if sweep_confs else None,
                      "p90": _percentile(sorted(sweep_confs), 90)
                                                if sweep_confs else None,
                      "p99": _percentile(sorted(sweep_confs), 99)
                                                if sweep_confs else None},
            "loco":  {"n": len(loco_confs),  "p50": _percentile(sorted(loco_confs), 50)
                                                if loco_confs else None,
                      "p90": _percentile(sorted(loco_confs), 90)
                                                if loco_confs else None,
                      "p99": _percentile(sorted(loco_confs), 99)
                                                if loco_confs else None},
        }

    for cls in EXCLUDED:
        grids[cls] = {"sweep": [], "loco": [], "kind": "excluded_at_species",
                      "default_tau": None,
                      "note": "Confidence semantics degenerate at species rank; excluded."}

    out_fp = Path(args.out)
    out_fp.parent.mkdir(parents=True, exist_ok=True)
    with open(out_fp, "w") as fh:
        json.dump({"grids": grids, "distributions": distributions}, fh, indent=2)
    print(f"[build_tau_grid] wrote {out_fp}")
    for cls in CONTINUOUS + BIMODAL:
        d = distributions[cls]
        print(f"  {cls:12s} sweep n={d['sweep']['n']:>6d} p50={d['sweep']['p50']}  "
              f"loco n={d['loco']['n']:>6d} p50={d['loco']['p50']}")


if __name__ == "__main__":
    main()
