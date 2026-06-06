#!/usr/bin/env python3
"""Step 04b - Confidence-threshold sweep evaluator (Reviewer 2 comment 2.3).

For every (mock x region x classifier) under a sweep results tree, build a
P / R / F-vs-confidence-threshold curve at every taxonomic level.

Two cases:

* Per-rank-confidence tools (taxonova, idtaxa, kraken2, mapseq, spingo):
  one permissive run is enough — the per-rank Confidence column in
  taxonomy_final.tsv is reinterpreted at every threshold T in the sweep
  grid by truncating the predicted lineage to the deepest rank whose
  reported per-rank confidence is >= T. This is the natural extension of
  what each tool does internally.

* Single-confidence tools (naive-bayes, vsearch, blast):
  classify-sklearn / classify-consensus-{vsearch,blast} already truncate
  ranks below the configured threshold and report a single Confidence /
  Consensus value for the whole prediction. We therefore run them at
  multiple actual confidence settings (one prediction directory per
  threshold) and the param-dir name itself supplies the threshold.

Output: results/sweep/evaluation/all_metrics.csv (long format,
columns: mock, region, classifier, threshold, level, precision, recall,
f-score, detection_rate, n_pred, n_true).

Run after 03_assign with --params-config configs/classifiers_sweep.json
--results-dir results/sweep/<mock>.
"""

import argparse
import json
import re
from collections import Counter
from pathlib import Path

import pandas as pd

import _common

# Import the sophisticated rank-by-rank evaluator used by 04_evaluate.py so
# that classifiers which leave upper ranks empty (notably SPINGO, whose
# output starts populated at genus) are scored consistently with the rest of
# the manuscript.
_common.add_src_to_path()
from evaluate_results import precision_recall_fscore  # noqa: E402


LEVELS = ["Class", "Order", "Family", "Genus", "Species"]
LEVEL_TO_INT = {"Class": 2, "Order": 3, "Family": 4, "Genus": 5, "Species": 6}

# Tools whose taxonomy_final.tsv carries a semicolon-separated per-rank
# Confidence column. We post-hoc threshold-sweep these from a single
# permissive run.
PER_RANK_TOOLS = {"taxonova", "idtaxa", "kraken2", "mapseq", "spingo"}

# Tools that emit one Confidence/Consensus value for the whole prediction.
# We re-run these at several actual thresholds; each param directory
# becomes one curve point.
SINGLE_CONF_TOOLS = {"naive-bayes", "vsearch", "blast"}

# Param-name -> threshold extractor for single-confidence tools.
SINGLE_CONF_PARAM_KEY = {
    "naive-bayes": "p-confidence",
    "vsearch":     "p-min-consensus",
    "blast":       "p-min-consensus",
}

# Fallback sweep grid for per-rank tools when no per-classifier grid is loaded.
# Dense at the low end so tools whose confidence semantics produce small
# values (e.g. Kraken 2's k-mer-fraction LCA score, typically 0.01–0.05
# for short amplicons) still produce non-trivial curves.
SWEEP_GRID = sorted({
    0.000, 0.001, 0.002, 0.005, 0.010, 0.015, 0.020, 0.030, 0.050,
    0.075, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 0.400, 0.450,
    0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.800, 0.850, 0.900,
    0.950, 0.990, 1.000,
})


def _load_tau_grids(path: Path | None):
    """Load configs/tau_grids.json if present. Returns {classifier: [tau, ...]}.
    Falls back to SWEEP_GRID for any classifier not in the file."""
    if path is None or not path.exists():
        return {}
    with open(path) as fh:
        doc = json.load(fh)
    out = {}
    for cls, entry in doc.get("grids", {}).items():
        grid = entry.get("sweep") or []
        if grid:
            out[cls] = sorted({round(float(t), 4) for t in grid})
    return out


def _find_trueish(data_root: Path, mock: str, region: str) -> Path | None:
    candidates = [
        data_root / "mock_community" / mock / region / "trueish-taxonomies.tsv",
        data_root / "mock_community" / mock / region / f"trueish_taxonomy_{mock}.tsv",
        data_root / "mock_community" / mock / region / f"trueish-taxonomy_{mock}.tsv",
    ]
    return next((p for p in candidates if p.exists()), None)


def _read_tax(fp: Path):
    """Return list of (feature_id, [rank_strings], [rank_confidences])."""
    rows = []
    with open(fp) as fh:
        header = fh.readline()
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            fid = parts[0]
            tax = parts[1]
            confs_raw = parts[2] if len(parts) >= 3 else ""
            ranks = [r.strip() for r in tax.split(";")]
            confs = []
            if confs_raw:
                for c in confs_raw.split(";"):
                    try:
                        confs.append(float(c.strip()))
                    except ValueError:
                        confs.append(0.0)
            rows.append((fid, ranks, confs))
    return rows


def _read_true(fp: Path) -> dict[str, list[str]]:
    out = {}
    with open(fp) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] in ("id", "Feature ID", "taxonomy", "Taxon"):
                continue
            if len(parts) >= 2:
                out[parts[0]] = [r.strip() for r in parts[1].split(";")]
    return out


_RANK_PREFIXES = ("d__", "p__", "c__", "o__", "f__", "g__", "s__")


def _strip_prefix(s: str) -> str:
    s = s.strip()
    return s[3:] if s[:3] in _RANK_PREFIXES else s


def _is_assigned(content: str) -> bool:
    if not content:
        return False
    if content in ("Unclassified", "Unassigned", "No blast hit", "AMBIGUOUS"):
        return False
    if content.startswith("unclassified_"):
        return False
    return True


def _truncate_at_threshold(ranks: list[str], confs: list[float], threshold: float) -> list[str]:
    """Apply the threshold to each rank position independently and return a
    length-7 list of rank labels.

    A rank is *retained* iff its label is assigned and its per-rank
    confidence >= threshold. A rank that fails either test is replaced
    with the empty rank prefix (e.g. 'g__'). We do **not** truncate a
    whole tail on the first failure: tools such as SPINGO report empty
    upper ranks (d__..f__) but valid genus/species, and the rank-by-rank
    matcher in evaluate_results.precision_recall_fscore is happy to score
    a non-contiguous lineage as long as the deepest assigned ranks agree
    with the expected lineage.

    This is a strict superset of the simpler "truncate the tail" rule:
    for monotone-confidence tools (taxonova caps child <= parent
    confidence), the two rules coincide because a low-confidence parent
    always implies a low-confidence child.
    """
    out = []
    for i in range(7):
        if i < len(ranks):
            label = ranks[i].strip()
            content = _strip_prefix(label)
        else:
            content = ""
        conf = confs[i] if i < len(confs) else 0.0
        if _is_assigned(content) and conf >= threshold:
            out.append(f"{_RANK_PREFIXES[i]}{content}")
        else:
            out.append(_RANK_PREFIXES[i])
    return out


def _truncate_to_level(ranks: list[str], level_int: int) -> str:
    """Take the lineage truncated by threshold and collapse to level_int.
    Empty content at level_int -> 'unassigned'."""
    if level_int >= len(ranks):
        return "unassigned"
    label = ranks[level_int].strip()
    content = _strip_prefix(label)
    if not _is_assigned(content):
        # use the deepest assigned ancestor
        for j in range(level_int - 1, -1, -1):
            anc = _strip_prefix(ranks[j].strip())
            if _is_assigned(anc):
                return anc
        return "unassigned"
    return content


def _evaluate_threshold(true_map, pred_rows, threshold):
    """Return dict level_name -> (precision, recall, f, det, n_pred, n_true).

    Precision/recall/f-score use evaluate_results.precision_recall_fscore
    (rank-by-rank match / over- / under- / mis-classification, same
    semantics as 04_evaluate.py). detection_rate is unique-taxon-set
    intersection at the deepest non-empty rank up to the evaluated level.
    """
    out = {}
    for level_name, lvl_int in LEVEL_TO_INT.items():
        exp_taxa = []
        obs_taxa = []
        for fid, ranks, confs in pred_rows:
            if fid not in true_map:
                continue
            true_ranks = true_map[fid]
            kept = _truncate_at_threshold(ranks, confs, threshold)
            obs_full = "; ".join(kept[: lvl_int + 1])
            true_padded = list(true_ranks) + [_RANK_PREFIXES[i] for i in range(len(true_ranks), 7)]
            exp_full = "; ".join(t.strip() for t in true_padded[: lvl_int + 1])
            exp_taxa.append(exp_full)
            obs_taxa.append(obs_full)
        _, _, prec, rec, f, _ = precision_recall_fscore(exp_taxa, obs_taxa)
        # detection rate (taxon-set) — collapse each lineage to its deepest
        # non-empty content within ranks 0..lvl_int (this matches what a
        # downstream user gets when querying "which taxa did the tool
        # detect?" and gives SPINGO-style genus/species-only outputs a fair
        # shake at deeper levels).
        def _deepest_assigned(lineage):
            parts = [p.strip() for p in lineage.split(";")]
            for j in range(min(lvl_int, len(parts) - 1), -1, -1):
                content = _strip_prefix(parts[j])
                if _is_assigned(content):
                    return f"{_RANK_PREFIXES[j]}{content}"
            return None

        def _resolve_obs_deepest(exp_d, obs_d):
            """Bracket-aware row-level resolution mirroring
            05_realworld_evaluate.py::_resolve_group_keys: if the truth
            taxon at the deepest-obs rank is a member of an obs bracket,
            collapse to truth; else keep bracket as one bucket."""
            if obs_d is None or "[" not in obs_d:
                return obs_d
            # Expand only the bracketed prefix__[a|b|c] to candidates.
            m = re.match(r"^([dpcofgs]__)\[([^\]]+)\]$", obs_d)
            if not m:
                return obs_d
            prefix = m.group(1)
            members = [x for x in m.group(2).split("|") if x]
            candidates = {prefix + m_ for m_ in members}
            if exp_d in candidates:
                return exp_d
            return obs_d

        exp_deep_list = [_deepest_assigned(e) for e in exp_taxa]
        obs_deep_list = [_deepest_assigned(o) for o in obs_taxa]
        obs_deep_resolved = [
            _resolve_obs_deepest(e, o) for e, o in zip(exp_deep_list, obs_deep_list)
        ]
        exp_set = {x for x in exp_deep_list if x is not None}
        obs_set = {x for x in obs_deep_resolved if x is not None}
        det = round(len(exp_set & obs_set) / max(1, len(exp_set)), 4)
        n_pred = sum(1 for o in obs_deep_list if o is not None)
        n_true = sum(1 for e in exp_deep_list if e is not None)
        # Abstention at this level = predictions with no assigned content at lvl_int.
        n_attempted = len(obs_taxa)
        n_abstained = 0
        for o in obs_taxa:
            parts = [p.strip() for p in o.split(";")]
            if lvl_int >= len(parts) or not _is_assigned(_strip_prefix(parts[lvl_int])):
                n_abstained += 1
        abstain_rate = round(n_abstained / max(1, n_attempted), 4)
        out[level_name] = (prec, rec, f, det, n_pred, n_true, n_abstained, abstain_rate)
    return out


def _param_threshold_singleconf(classifier: str, param_id: str) -> float | None:
    """Parse the threshold value out of a single-confidence param dir name
    like 'p-confidence-0.7' or 'p-maxaccepts-10_p-min-consensus-0.8_p-perc-identity-0.8'."""
    key = SINGLE_CONF_PARAM_KEY.get(classifier)
    if not key:
        return None
    m = re.search(rf"{re.escape(key)}-([0-9.]+)", param_id)
    if not m:
        return None
    return float(m.group(1))


def _walk_sweep_tree(sweep_root: Path, mock: str, regions: list[str] | None):
    """Yield (region, classifier, param_id, pred_fp) tuples."""
    mock_dir = sweep_root / mock
    if not mock_dir.is_dir():
        return
    for region_dir in sorted(mock_dir.iterdir()):
        if not region_dir.is_dir():
            continue
        if regions and region_dir.name not in regions:
            continue
        for cls_dir in sorted(region_dir.iterdir()):
            if not cls_dir.is_dir():
                continue
            for param_dir in sorted(cls_dir.iterdir()):
                pred_fp = param_dir / "taxonomy_final.tsv"
                if pred_fp.exists():
                    yield region_dir.name, cls_dir.name, param_dir.name, pred_fp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep-results-dir", required=True,
                    help="Root holding sweep results, e.g. results/sweep/")
    ap.add_argument("--mocks", nargs="+", required=True)
    ap.add_argument("--regions", nargs="+", default=None)
    ap.add_argument("--out-dir", default=None,
                    help="Output dir; default <sweep-results-dir>/evaluation")
    ap.add_argument("--data-dir", default=None)
    ap.add_argument("--tau-grid", default=str(Path(__file__).resolve().parents[1] /
                                              "configs" / "tau_grids.json"),
                    help="JSON with per-classifier sweep tau grids.")
    args = ap.parse_args()

    cfg = _common.load_pipeline_config()
    data_root = Path(args.data_dir) if args.data_dir else Path(cfg["data_dir"])
    sweep_root = Path(args.sweep_results_dir)
    out_dir = Path(args.out_dir) if args.out_dir else sweep_root / "evaluation"
    out_dir.mkdir(parents=True, exist_ok=True)

    per_cls_grids = _load_tau_grids(Path(args.tau_grid) if args.tau_grid else None)
    if per_cls_grids:
        print(f"[04b] loaded per-classifier tau grids for: {sorted(per_cls_grids)}")

    rows = []
    for mock in args.mocks:
        for region, cls, param_id, pred_fp in _walk_sweep_tree(sweep_root, mock, args.regions):
            true_fp = _find_trueish(data_root, mock, region)
            if true_fp is None:
                print(f"  - skip {mock}/{region}: no trueish file")
                continue
            true_map = _read_true(true_fp)
            pred_rows = _read_tax(pred_fp)

            if cls in PER_RANK_TOOLS:
                # post-hoc sweep — per-classifier percentile grid if available,
                # else the dense fallback SWEEP_GRID.
                thresholds = per_cls_grids.get(cls, SWEEP_GRID)
                for T in thresholds:
                    metrics = _evaluate_threshold(true_map, pred_rows, T)
                    for lvl, (p, r, f, det, np_, nt, n_abs, abs_rate) in metrics.items():
                        rows.append({
                            "mock": mock, "region": region, "classifier": cls,
                            "threshold": T, "level": lvl,
                            "precision": p, "recall": r, "f_score": f,
                            "detection_rate": det,
                            "n_predictions": np_, "n_true_assigned": nt,
                            "n_abstained": n_abs, "abstain_rate": abs_rate,
                            "param_id": param_id, "sweep_kind": "post_hoc",
                        })
            elif cls in SINGLE_CONF_TOOLS:
                T = _param_threshold_singleconf(cls, param_id)
                if T is None:
                    print(f"  - {mock}/{region}/{cls}/{param_id}: cannot parse threshold")
                    continue
                # threshold=T applied at run time; evaluate single point with
                # threshold=0 (all reported labels are already gated by T).
                metrics = _evaluate_threshold(true_map, pred_rows, 0.0)
                for lvl, (p, r, f, det, np_, nt, n_abs, abs_rate) in metrics.items():
                    rows.append({
                        "mock": mock, "region": region, "classifier": cls,
                        "threshold": T, "level": lvl,
                        "precision": p, "recall": r, "f_score": f,
                        "detection_rate": det,
                        "n_predictions": np_, "n_true_assigned": nt,
                        "n_abstained": n_abs, "abstain_rate": abs_rate,
                        "param_id": param_id, "sweep_kind": "rerun",
                    })
            else:
                print(f"  - skip {cls}: not configured for sweep")

    if not rows:
        print("[04b] no rows produced")
        return

    df = pd.DataFrame(rows)
    long_fp = out_dir / "all_metrics.csv"
    df.to_csv(long_fp, index=False)
    print(f"[04b] wrote {long_fp} ({len(df):,} rows)")

    # Tidy per-(metric, level) wide tables for quick inspection
    by_metric = out_dir / "by_metric"
    by_metric.mkdir(exist_ok=True)
    for metric in ["precision", "recall", "f_score", "detection_rate"]:
        for lvl in LEVELS:
            sub = df[df["level"] == lvl]
            if sub.empty:
                continue
            piv = sub.pivot_table(
                index=["mock", "region", "classifier"],
                columns="threshold", values=metric, aggfunc="mean")
            piv.to_csv(by_metric / f"{lvl}_{metric}.csv")
    print(f"[04b] wrote per-(level x metric) wide tables -> {by_metric}")


if __name__ == "__main__":
    main()
