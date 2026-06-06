#!/usr/bin/env python3
"""Step 04 - Evaluate predicted vs trueish taxonomies.

Auto-discovers every (mock x region) under results/ that has predictions and a
trueish-taxonomies file, computes precision / recall / f-score / detection-rate /
sensitivity / specificity / L1 / L2 at each taxonomic level, and writes a single
consolidated tree under `results/evaluation/`:

Metrics:
    precision        — per-feature TP / (TP + FP)
    recall           — per-feature TP / (TP + FN) (matches the manuscript's "Recall")
    f-score          — 2 P R / (P + R)
    detection_rate   — taxon-set TP / (TP + FN), i.e. fraction of true taxa detected
    sensitivity      — taxon-set TP / (TP + FN); mathematically identical to
                       `detection_rate`. Provided as a separate column with the
                       reviewer-requested name (Reviewer 1, comment 1.9).
    specificity      — taxon-set TN / (TN + FP) where the universe is the unique
                       taxa at this level in the per-region reference DB
                       (Greengenes2 backbone for V4/V3-V4/V5-V7/V1-V9; GROND for
                       RRN). Quantifies the rate at which a tool correctly
                       refrains from predicting taxa that are absent from the
                       ground truth (Reviewer 1, comment 1.9).
    l1_error         — sum_t |freq_exp_t - freq_obs_t| (abundance distance)
    l2_error         — sum_t (freq_exp_t - freq_obs_t)^2

    results/evaluation/all_metrics.csv          # tidy long format
    results/evaluation/by_metric/<metric>.csv   # mock,region,method,param x level

Re-run any time params change. The evaluation tree is overwritten in full;
prediction directories are not touched.

Optional --mock-name / --regions flags restrict the scan for debugging.
"""

import argparse
import re
import shutil
from collections import Counter
from pathlib import Path

import pandas as pd

import _common


# Match `<rank>__[a|b|c]` group labels emitted by taxonova (comment-2.9
# species-group models). Matches at any rank prefix (d/p/c/o/f/g/s).
# Identical to the regex used by 05_realworld_evaluate.py — keep in sync.
_GROUP_RE = re.compile(r"([dpcofgs]__)\[([^\]]+)\]")


def _expand_lineage_brackets(tax_str):
    """Expand a truncated taxonomy string containing one or more
    ``<rank>__[a|b|c]`` group labels into the list of single-member
    taxonomy strings. Bracketless strings return ``[tax_str]``. Mirrors
    `05_realworld_evaluate.py::_expand_group_label`."""
    if "[" not in tax_str:
        return [tax_str]
    expansions = [tax_str]
    for m in _GROUP_RE.finditer(tax_str):
        token = m.group(0)
        prefix = m.group(1)
        members = [x for x in m.group(2).split("|") if x]
        if not members:
            continue
        new = []
        for e in expansions:
            for memb in members:
                new.append(e.replace(token, prefix + memb, 1))
        expansions = new
    return expansions


def _resolve_obs_against_truth(exp_full, obs_full):
    """Bracket-aware per-row resolution mirroring
    `05_realworld_evaluate.py::_resolve_group_keys` semantics for the
    mock/biome (1 ASV = 1 prediction unit) case.

    If `obs_full` contains no group brackets, return it unchanged.
    Otherwise expand brackets to the candidate single-member strings:
      * If `exp_full` is among the candidates → collapse to `exp_full`
        (full TP credit goes to the truth taxon).
      * Else keep the original bracket as a single bucket — splitting
        an unmatched bracket into N members would inflate the FP/L1/L2
        count by N for one prediction (comment-2.9 invariant: 1 row
        contributes 1 prediction unit regardless of bracket cardinality).
    """
    if "[" not in obs_full:
        return obs_full
    candidates = _expand_lineage_brackets(obs_full)
    if exp_full in candidates:
        return exp_full
    return obs_full


LEVELS = ["Class", "Order", "Family", "Genus", "Species"]
LEVEL_TO_INT = {"Class": 2, "Order": 3, "Family": 4, "Genus": 5, "Species": 6}
METRIC_TO_COL = {
    "precision":      "Precision",
    "recall":         "Recall",
    "f-score":        "F-measure",
    "detection_rate": "Taxon Detection Rate",
    "sensitivity":    "Sensitivity",
    "specificity":    "Specificity",
    "l1_error":       "L1 Abundance Error",
    "l2_error":       "L2 Abundance Error",
}


def _ref_taxonomy_fp(region):
    """Path to the per-region reference taxonomy TSV used to define the
    'universe' of possible taxa for specificity (TN / (TN + FP))."""
    cfg = _common.load_pipeline_config()
    ref_root = Path(cfg["reference_database_dir"])
    if region == "RRN":
        return ref_root / "RRN" / "grond_clean.tsv"
    return ref_root / region / "taxonomy.tsv"


_UNIVERSE_CACHE = {}


def _universe_size(region, level_int):
    """Number of distinct taxa at `level_int` in the per-region reference DB.
    This is the universe used for specificity (TN = N - |true ∪ pred|)."""
    key = (region, level_int)
    if key in _UNIVERSE_CACHE:
        return _UNIVERSE_CACHE[key]
    fp = _ref_taxonomy_fp(region)
    if not fp.exists():
        _UNIVERSE_CACHE[key] = None
        return None
    seen = set()
    with open(fp) as fh:
        next(fh, None)  # header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            seen.add(_truncate_tax(parts[1], level_int))
    _UNIVERSE_CACHE[key] = len(seen)
    return len(seen)


def _abundance_norms(exp_taxa, obs_taxa):
    """L1 = sum |p_exp - p_obs|, L2 = sum (p_exp - p_obs)^2 over the union of taxa.

    Frequencies are derived by counting truncated-taxonomy occurrences in each
    list (one entry per feature) and dividing by the list length. Matches the
    old taxobench_11 species-level definition (L2 kept as squared, not sqrt).
    """
    n_exp = max(1, len(exp_taxa))
    n_obs = max(1, len(obs_taxa))
    cexp = Counter(exp_taxa)
    cobs = Counter(obs_taxa)
    l1 = 0.0
    l2 = 0.0
    for t in set(cexp) | set(cobs):
        diff = cexp.get(t, 0) / n_exp - cobs.get(t, 0) / n_obs
        l1 += abs(diff)
        l2 += diff * diff
    return round(l1, 6), round(l2, 6)


def _read_tax_tsv(fp):
    out = {}
    with open(fp) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] in ("id", "Feature ID", "taxonomy", "Taxon"):
                continue
            if len(parts) >= 2:
                out[parts[0]] = parts[1]
    return out


def _load_taxa_pairs(true_fp, pred_fp):
    t = _read_tax_tsv(true_fp)
    p = _read_tax_tsv(pred_fp)
    return {k: (v, p.get(k, "Unassigned")) for k, v in t.items()}


def _truncate_tax(tax_str, level_int):
    parts = [x.strip() for x in tax_str.split(";")]
    return ";".join(parts[:level_int + 1])


def _find_trueish(data_root, mock, region):
    candidates = [
        data_root / "mock_community" / mock / region / "trueish-taxonomies.tsv",
        data_root / "mock_community" / mock / region / f"trueish_taxonomy_{mock}.tsv",
        data_root / "mock_community" / mock / region / f"trueish-taxonomy_{mock}.tsv",
    ]
    return next((p for p in candidates if p.exists()), None)


def _discover(res_root, data_root, mock_filter, region_filter):
    """Yield (mock, region, region_dir, true_fp) for every evaluable combo."""
    for mock_dir in sorted(res_root.iterdir()):
        if not mock_dir.is_dir() or mock_dir.name in ("evaluation", "real_world", "runtime"):
            continue
        if mock_filter and mock_dir.name != mock_filter:
            continue
        for region_dir in sorted(mock_dir.iterdir()):
            if not region_dir.is_dir():
                continue
            if region_filter and region_dir.name not in region_filter:
                continue
            true_fp = _find_trueish(data_root, mock_dir.name, region_dir.name)
            if true_fp is None:
                print(f"  - skip {mock_dir.name}/{region_dir.name}: no trueish file")
                continue
            yield mock_dir.name, region_dir.name, region_dir, true_fp


def _evaluate_combo(mock, region, region_dir, true_fp, prf_fn):
    rows = []
    for pred_fp in sorted(region_dir.glob("*/*/taxonomy_final.tsv")):
        method = pred_fp.parent.parent.name
        param  = pred_fp.parent.name
        try:
            pairs = _load_taxa_pairs(str(true_fp), str(pred_fp))
        except Exception as exc:
            print(f"  ! {mock}/{region}/{method}/{param}: load failed ({exc})")
            continue
        for level_name, lvl_int in LEVEL_TO_INT.items():
            exp_taxa, obs_taxa = [], []
            obs_taxa_resolved = []
            for _id, (te, to) in pairs.items():
                exp_full = _truncate_tax(te, lvl_int)
                obs_full = _truncate_tax(to, lvl_int)
                exp_taxa.append(exp_full)
                obs_taxa.append(obs_full)
                # Bracket-aware resolution: a bracket whose members
                # contain the truth collapses to the truth taxon; an
                # unmatched bracket stays as one opaque bucket (1 row =
                # 1 prediction unit). Used for taxon-set ops and L1/L2;
                # NOT used for prf_fn (which is already bracket-aware
                # via evaluate_classification / _expand_group_label).
                obs_taxa_resolved.append(_resolve_obs_against_truth(exp_full, obs_full))
            _, _, p, r, f, _ = prf_fn(exp_taxa, obs_taxa)
            exp_set, obs_set = set(exp_taxa), set(obs_taxa_resolved)
            det = round(len(exp_set & obs_set) / max(1, len(exp_set)), 3)
            l1, l2 = _abundance_norms(exp_taxa, obs_taxa_resolved)

            # Taxon-set sensitivity / specificity (Selten reviewer 1.9):
            # universe = unique taxa at this level in the per-region reference DB
            tp_t = len(exp_set & obs_set)
            fp_t = len(obs_set - exp_set)
            fn_t = len(exp_set - obs_set)
            sens = round(tp_t / (tp_t + fn_t), 4) if (tp_t + fn_t) else 0.0
            N = _universe_size(region, lvl_int)
            if N is None or N <= 0:
                spec = ""
            else:
                tn_t = max(0, N - len(exp_set | obs_set))
                spec = round(tn_t / (tn_t + fp_t), 6) if (tn_t + fp_t) else 0.0
            rows.append({
                "mock": mock, "region": region,
                "method": method, "param": param,
                "level": level_name,
                "precision": p, "recall": r,
                "f-score": f, "detection_rate": det,
                "sensitivity": sens, "specificity": spec,
                "l1_error": l1, "l2_error": l2,
            })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mock-name", default=None,
                    help="Restrict to one mock (default: all under results/)")
    ap.add_argument("--regions", nargs="+", default=None,
                    help="Restrict to specific regions (default: all found)")
    ap.add_argument("--results-dir", default=None,
                    help="Override results root (default: cfg.results_dir)")
    ap.add_argument("--data-dir", default=None,
                    help="Override data root (default: cfg.data_dir)")
    args = ap.parse_args()

    cfg = _common.load_pipeline_config()
    _common.add_src_to_path()
    import evaluate_results

    res_root  = Path(args.results_dir) if args.results_dir else Path(cfg["results_dir"])
    data_root = Path(args.data_dir)    if args.data_dir    else Path(cfg["data_dir"])
    out_root  = res_root / "evaluation"

    print(f"[04_evaluate] scanning {res_root}")
    all_rows = []
    for mock, region, region_dir, true_fp in _discover(
            res_root, data_root, args.mock_name, args.regions):
        rows = _evaluate_combo(mock, region, region_dir, true_fp,
                               evaluate_results.precision_recall_fscore)
        if not rows:
            print(f"  - {mock}/{region}: no predictions")
            continue
        print(f"  + {mock}/{region}: {len({(r['method'], r['param']) for r in rows})} (method,param) combos")
        all_rows.extend(rows)

    if not all_rows:
        print("[04_evaluate] no evaluable combos found")
        return

    # Overwrite the evaluation tree atomically: stage to .tmp, swap, drop old.
    tmp_root = out_root.with_suffix(".tmp")
    if tmp_root.exists():
        shutil.rmtree(tmp_root)
    (tmp_root / "by_metric").mkdir(parents=True)

    df = pd.DataFrame(all_rows)
    long = df.melt(
        id_vars=["mock", "region", "method", "param", "level"],
        value_vars=list(METRIC_TO_COL.keys()),
        var_name="metric", value_name="value",
    )
    long_fp = tmp_root / "all_metrics.csv"
    long.to_csv(long_fp, index=False)
    print(f"[04_evaluate] wrote {long_fp} ({len(long):,} rows)")

    for metric in METRIC_TO_COL:
        pivot = df.pivot_table(
            index=["mock", "region", "method", "param"],
            columns="level", values=metric, aggfunc="mean",
        )
        pivot = pivot.reindex(columns=LEVELS)  # consistent column order
        fp = tmp_root / "by_metric" / f"{metric}.csv"
        pivot.to_csv(fp)
        print(f"[04_evaluate] wrote {fp}")

    if out_root.exists():
        shutil.rmtree(out_root)
    tmp_root.rename(out_root)
    print(f"[04_evaluate] done -> {out_root}")


if __name__ == "__main__":
    main()
