#!/usr/bin/env python3
"""Step 04 - Evaluate predicted vs trueish taxonomies.

Auto-discovers every (mock x region) under results/ that has predictions and a
trueish-taxonomies file, computes precision / recall / f-score / detection-rate /
sensitivity / specificity / L1 / L2 at each taxonomic level, and writes a single
consolidated tree under `results/evaluation/`:

Metrics:
    precision       : per-feature TP / (TP + FP)
    recall          : per-feature TP / (TP + FN) (matches the manuscript's "Recall")
    f-score         : 2 P R / (P + R)
    detection_rate  : taxon-set TP / (TP + FN), i.e. fraction of true taxa detected
    sensitivity     : taxon-set TP / (TP + FN); mathematically identical to
                       `detection_rate`. Provided as a separate column under
                       the conventional name.
    specificity     : taxon-set TN / (TN + FP) where the universe is the unique
                       taxa at this level in the per-region reference DB
                       (Greengenes2 backbone for V4/V3-V4/V5-V7/V1-V9; GROND for
                       RRN). Quantifies the rate at which a tool correctly
                       refrains from predicting taxa that are absent from the
                       ground truth.
    l1_error        : sum_t |freq_exp_t - freq_obs_t| (abundance distance)
    l2_error        : sum_t (freq_exp_t - freq_obs_t)^2

    results/evaluation/all_metrics.csv          # tidy long format
    results/evaluation/by_metric/<metric>.csv   # mock,region,method,param x level

Re-run any time params change. The evaluation tree is overwritten in full;
prediction directories are not touched.

Optional --mock-name / --regions flags restrict the scan for debugging.

Fail closed (2026-10-31). Every taxonomy_final.tsv must hold each truth ID
exactly once (no missing, unexpected or duplicated IDs), and when a
run.status.json from 03_assign.py sits beside it, that run must have succeeded
and the file must be the one it wrote. Any violation, any load error, or
finding nothing to score stops the script with exit status 1 before anything
is written. None of this changes how a valid output is scored.

Scoring rules (--rules; details in src/scoring_v2.py and README.md):
    v2 (default)  an abstention is no call (FN in P/R/F, never a taxon in the
                  taxon-set or L1/L2 metrics); a group label is correct iff
                  the truth is a member; taxon identity is rank prefix + name.
    v2-strict     as v2, but a multi-member group is no call.
    v1            the published rules; with the v1 truth (--truth-dir, see
                  README.md) this reproduces the published values.
    v1-patched    v1 with abstentions dropped from the taxon-set metrics only.
"""

import argparse
import re
import shutil
import sys
from collections import Counter
from pathlib import Path

import pandas as pd

import _common


# Match `<rank>__[a|b|c]` group labels emitted by taxonova (comment-2.9
# species-group models). Matches at any rank prefix (d/p/c/o/f/g/s).
# Identical to the regex used by 05_realworld_evaluate.py: keep in sync.
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
      * Else keep the original bracket as a single bucket: splitting
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
    ref_root = Path(_REF_DIR_OVERRIDE or cfg["reference_database_dir"])
    if not ref_root.exists():
        # The configured path does not exist (for example the placeholder in
        # configs/pipeline_config.json); fall back to this checkout.
        ref_root = _common.project_root() / "data" / "ref_dbs"
    if region == "RRN":
        return ref_root / "RRN" / "grond_clean.tsv"
    return ref_root / region / "taxonomy.tsv"


_UNIVERSE_CACHE = {}
_REF_DIR_OVERRIDE = None   # set by --ref-dir


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


class OutputError(Exception):
    """A prediction file that must not be scored (partial, stale, failed)."""


def _check_output(true_fp, pred_fp):
    """Raise OutputError unless `pred_fp` covers the truth IDs exactly and,
    when 03_assign.py left a run.status.json, comes from a successful run."""
    status = _common.read_run_status(pred_fp.parent)
    if status is not None:
        if not status.get("ok", False):
            raise OutputError("run.status.json reports a failed run: "
                              + "; ".join(status.get("problems") or ["(no detail)"]))
        if pred_fp.stat().st_mtime < float(status.get("started", 0)) - 1:
            raise OutputError("taxonomy_final.tsv predates the run recorded in run.status.json")
    truth_ids = _common.table_ids(true_fp)
    if len(truth_ids) != len(set(truth_ids)):
        raise OutputError(f"truth file {true_fp} has duplicated IDs")
    probs = _common.coverage_problems(_common.table_ids(pred_fp), truth_ids)
    if probs:
        raise OutputError("; ".join(probs))


def _truncate_tax(tax_str, level_int):
    parts = [x.strip() for x in tax_str.split(";")]
    return ";".join(parts[:level_int + 1])


def _is_abstention(tax_str, level_int):
    parts = [x.strip() for x in tax_str.split(";")]
    if len(parts) <= level_int:
        return True
    last = parts[level_int]
    return last in ("", "Unassigned", "Unclassified") or (len(last) == 3 and last.endswith("__"))


def _find_trueish(data_root, mock, region):
    candidates = [
        data_root / "mock_community" / mock / region / "trueish-taxonomies.tsv",
        data_root / "mock_community" / mock / region / f"trueish_taxonomy_{mock}.tsv",
        data_root / "mock_community" / mock / region / f"trueish-taxonomy_{mock}.tsv",
    ]
    return next((p for p in candidates if p.exists()), None)


def _discover(res_root, data_root, mock_filter, region_filter, missing=None):
    """Yield (mock, region, region_dir, true_fp) for every evaluable combo.
    Cells with predictions but no truth file are appended to `missing`."""
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
                if missing is not None and any(region_dir.glob("*/*/taxonomy_final.tsv")):
                    missing.append(f"{mock_dir.name}/{region_dir.name}: predictions "
                                   "but no trueish taxonomy file")
                continue
            yield mock_dir.name, region_dir.name, region_dir, true_fp


def _v1_level_row(pairs, region, lvl_int, prf_fn, drop_abstentions):
    """Published (v1) metric computation for one level. With
    drop_abstentions=False this is the published code exactly; True is the
    2026-10-29 diagnosis patch (abstentions out of taxon-set and L1/L2)."""
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
    if drop_abstentions:
        # An abstention at this level ("...;g__X;s__", or "Unassigned") is not a
        # predicted taxon: keep it out of the taxon-set and abundance metrics so it
        # is not counted as a false-positive taxon (P/R/F already score it as FN).
        obs_taxa_resolved = [o for o in obs_taxa_resolved if not _is_abstention(o, lvl_int)]
    _, _, p, r, f, _ = prf_fn(exp_taxa, obs_taxa)
    exp_set, obs_set = set(exp_taxa), set(obs_taxa_resolved)
    det = round(len(exp_set & obs_set) / max(1, len(exp_set)), 3)
    l1, l2 = _abundance_norms(exp_taxa, obs_taxa_resolved)

    # Taxon-set sensitivity / specificity:
    # universe = unique taxa at this level in the per-region reference DB
    tp_t = len(exp_set & obs_set)
    fp_t = len(obs_set - exp_set)
    fn_t = len(exp_set - obs_set)
    sens = round(tp_t / (tp_t + fn_t), 4) if (tp_t + fn_t) else 0.0
    N = _universe_size(region, lvl_int)
    if N is None or N <= 0:
        spec = float("nan")  # was "": a blank string makes the by_metric pivot crash
    else:
        tn_t = max(0, N - len(exp_set | obs_set))
        spec = round(tn_t / (tn_t + fp_t), 6) if (tn_t + fp_t) else 0.0
    union = exp_set | obs_set
    return {
        "precision": p, "recall": r,
        "f-score": f, "detection_rate": det,
        "sensitivity": sens, "specificity": spec,
        "l1_error": l1, "l2_error": l2,
        "jaccard": round(tp_t / len(union), 4) if union else float("nan"),
        "fp_taxa": fp_t,
    }


def _evaluate_pairs(pairs, region, rules, prf_fn):
    """All levels for one (truth, prediction) pair set. Returns a list of
    per-level metric dicts (without the identifying columns)."""
    import scoring_v2
    out = []
    for level_name, lvl_int in LEVEL_TO_INT.items():
        if rules == "v1":
            row = _v1_level_row(pairs, region, lvl_int, prf_fn, drop_abstentions=False)
        elif rules == "v1-patched":
            row = _v1_level_row(pairs, region, lvl_int, prf_fn, drop_abstentions=True)
        else:
            row = scoring_v2.score_rows(pairs, lvl_int, rules,
                                        universe_size=_universe_size(region, lvl_int))
        row["level"] = level_name
        out.append(row)
    if rules != "v1":
        gs = scoring_v2.group_stats(pairs, 6)
        strict = scoring_v2.score_rows(pairs, 6, "v2-strict",
                                       universe_size=_universe_size(region, 6))
        for row in out:
            if row["level"] == "Species":
                row.update(gs)
                row["strict_species_f"] = strict["f-score"]
    return out


def _evaluate_combo(mock, region, region_dir, true_fp, prf_fn, rules="v2",
                    tie_hook=None, errors=None):
    """``tie_hook(mock, region, method, param, pred_fp, pairs) -> pairs`` may
    rewrite calls before scoring (for example to rebuild search-hit ties as
    group labels in a sensitivity analysis).

    An output that fails `_check_output` or cannot be loaded raises
    OutputError, or, when a list is passed as ``errors``, is appended to it
    and skipped (the caller then decides; main() exits 1)."""
    rows = []
    for pred_fp in sorted(region_dir.glob("*/*/taxonomy_final.tsv")):
        method = pred_fp.parent.parent.name
        param  = pred_fp.parent.name
        try:
            _check_output(Path(true_fp), pred_fp)
            pairs = _load_taxa_pairs(str(true_fp), str(pred_fp))
        except Exception as exc:
            msg = f"{mock}/{region}/{method}/{param}: {exc}"
            if errors is None:
                raise OutputError(msg) from exc
            print(f"  ! {msg}")
            errors.append(msg)
            continue
        if tie_hook is not None:
            pairs = tie_hook(mock, region, method, param, pred_fp, pairs)
        for row in _evaluate_pairs(pairs, region, rules, prf_fn):
            row.update({"mock": mock, "region": region,
                        "method": method, "param": param})
            rows.append(row)
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
    ap.add_argument("--truth-dir", default=None,
                    help="Data root to read trueish files from instead of --data-dir "
                         "(same <mock>/<region>/ layout; e.g. a corrected-truth tree)")
    ap.add_argument("--ref-dir", default=None,
                    help="Reference root holding <region>/taxonomy.tsv and RRN/grond_clean.tsv "
                         "(specificity universe); default cfg.reference_database_dir")
    ap.add_argument("--out-dir", default=None,
                    help="Write the evaluation tree here (default <results-dir>/evaluation)")
    ap.add_argument("--rules", default="v2", choices=["v1", "v1-patched", "v2", "v2-strict"],
                    help="v2 (default): abstention = no call, group label correct iff "
                         "truth is a member, taxon identity by rank name. "
                         "v1: published rules (reproduces the paper with the v1 truth). "
                         "v1-patched: v1 plus the diagnosis patch (abstentions dropped from "
                         "taxon-set metrics). v2-strict: v2 with "
                         "group labels treated as no call. See src/scoring_v2.py.")
    args = ap.parse_args()
    global _REF_DIR_OVERRIDE
    _REF_DIR_OVERRIDE = args.ref_dir

    cfg = _common.load_pipeline_config()
    _common.add_src_to_path()
    import evaluate_results

    res_root  = Path(args.results_dir) if args.results_dir else Path(cfg["results_dir"])
    data_root = Path(args.data_dir)    if args.data_dir    else Path(cfg["data_dir"])
    out_root  = Path(args.out_dir) if args.out_dir else res_root / "evaluation"
    truth_root = Path(args.truth_dir) if args.truth_dir else data_root

    print(f"[04_evaluate] scanning {res_root}")
    if not res_root.is_dir():
        print(f"[04_evaluate] FAILED: results directory {res_root} does not exist", file=sys.stderr)
        sys.exit(1)
    all_rows = []
    errors = []
    for mock, region, region_dir, true_fp in _discover(
            res_root, truth_root, args.mock_name, args.regions, missing=errors):
        rows = _evaluate_combo(mock, region, region_dir, true_fp,
                               evaluate_results.precision_recall_fscore,
                               rules=args.rules, errors=errors)
        if not rows:
            print(f"  - {mock}/{region}: no predictions")
            continue
        print(f"  + {mock}/{region}: {len({(r['method'], r['param']) for r in rows})} (method,param) combos")
        all_rows.extend(rows)

    if errors:
        print(f"[04_evaluate] FAILED: {len(errors)} output(s) cannot be scored; "
              "nothing written:", file=sys.stderr)
        for e in errors:
            print("  " + e, file=sys.stderr)
        sys.exit(1)
    if not all_rows:
        print("[04_evaluate] FAILED: no evaluable combos found", file=sys.stderr)
        sys.exit(1)

    # Overwrite the evaluation tree atomically: stage to .tmp, swap, drop old.
    tmp_root = out_root.with_suffix(".tmp")
    if tmp_root.exists():
        shutil.rmtree(tmp_root)
    (tmp_root / "by_metric").mkdir(parents=True)

    df = pd.DataFrame(all_rows)
    if args.rules == "v1":
        value_vars = list(METRIC_TO_COL.keys())
    else:
        value_vars = [c for c in df.columns
                      if c not in ("mock", "region", "method", "param", "level")]
    long = df.melt(
        id_vars=["mock", "region", "method", "param", "level"],
        value_vars=value_vars,
        var_name="metric", value_name="value",
    )
    long_fp = tmp_root / "all_metrics.csv"
    long.to_csv(long_fp, index=False)
    print(f"[04_evaluate] wrote {long_fp} ({len(long):,} rows)")

    for metric in METRIC_TO_COL:
        df[metric] = pd.to_numeric(df[metric], errors="coerce")
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
