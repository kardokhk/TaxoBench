#!/usr/bin/env python3
"""Step 05b - Evaluate real-world classifier predictions against paired
shotgun-WGS taxa.

For each (study × classifier × param) prediction at
`results/real_world/<study>/<classifier>/<param>/taxonomy_final.tsv`,
this script:

  1. Loads the amplicon ASV count table (samples × ASVs) — DADA2 output.
  2. Loads the classifier's per-ASV taxonomy assignment.
  3. Aggregates ASV counts → per-sample (rank, taxon) frequency vectors.
  4. Loads the paired shotgun OGU (operational genome unit) BIOM table and
     the OGU → taxonomy lookup, aggregates → per-sample (rank, taxon)
     frequency vectors with WGS as ground truth.
  5. Joins amplicon and WGS samples via `paired_samples.csv`
     (`amplicon_16s_run` ↔ `shotgun_wgs_run`).
  6. For each rank L ∈ {Phylum,…,Species}, computes:
       - Jaccard (taxon-set overlap)
       - Bray-Curtis dissimilarity (weighted)
       - precision, recall, F1 (per-sample taxon set)
       - L1, L2 abundance error (per-sample)
     Per-sample metrics are mean-aggregated per (study × classifier × param × level).
  7. Writes:
       - results/evaluation/realworld/all_metrics.csv (tidy long)
       - results/evaluation/realworld/by_metric/<metric>.csv (wide)
       - results/evaluation/realworld/per_sample/<study>_<classifier>_<param>.csv

Re-run any time taxonova (or any other classifier) is updated — produces a
fresh comparison without re-doing the expensive DADA2 / shotgun processing.
"""

import argparse
import json
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

# Match `<rank>__[a|b|c]` group labels emitted by taxonova (comment-2.9
# species-group models). Matches at any rank prefix (d/p/c/o/f/g/s).
_GROUP_RE = re.compile(r"([dpcofgs]__)\[([^\]]+)\]")

import _common
# `_unifrac_metric` lives in ../src (alongside `evaluate_results.py` and
# `process_taxonova_output.py`). Add it to sys.path here so scripts/ stays
# orchestration-only and src/ stays library-style.
import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
import _unifrac_metric


LEVELS = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]
LEVEL_TO_INT = {"Phylum": 1, "Class": 2, "Order": 3,
                "Family": 4, "Genus": 5, "Species": 6}


# --------------------------------------------------------------------------- #
# Loaders
# --------------------------------------------------------------------------- #

def _load_taxonomy_tsv(fp):
    """Return {feature_id: taxonomy_string}. Skips header rows."""
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


def _load_amplicon_table(fp):
    """asv_table_ids.csv is a CSV with first column = sample id, remaining
    columns = ASV_1, ASV_2, … with read counts. Returns DataFrame indexed
    by sample with ASV columns."""
    df = pd.read_csv(fp, index_col=0)
    df.index.name = "sample_id"
    return df


def _load_wgs_ogu_table(biom_fp):
    """Returns DataFrame indexed by OGU (observation), columns = sample IDs."""
    import biom
    t = biom.load_table(str(biom_fp))
    df = pd.DataFrame(
        t.matrix_data.toarray(),
        index=list(t.ids("observation")),
        columns=list(t.ids("sample")),
    )
    df.index.name = "ogu_id"
    return df


def _truncate_tax(tax_str, level_int):
    """Truncate a 'd__…; p__…; …' taxonomy to the given level (1=phylum…6=species)."""
    if not isinstance(tax_str, str) or not tax_str.strip():
        return ""
    parts = [x.strip() for x in tax_str.split(";")]
    return ";".join(parts[:level_int + 1])


def _per_sample_freqs(count_df, feature_to_tax, level_int):
    """count_df rows = samples, cols = features. feature_to_tax maps
    feature_id → taxonomy. Returns dict[sample_id -> dict[truncated_tax -> freq]].

    Group labels (e.g., taxonova's `s__[a|b|c]` from the species-group
    models, comment-2.9) are kept as-is here; truth-aware expansion is
    applied in `_resolve_group_keys` once the WGS truth set is known."""
    out = {}
    for sample_id, row in count_df.iterrows():
        bag = Counter()
        total = float(row.sum())
        if total <= 0:
            out[sample_id] = {}
            continue
        for feat, count in row.items():
            if count <= 0:
                continue
            tax = feature_to_tax.get(feat) or feature_to_tax.get(str(feat))
            if not tax:
                continue
            t = _truncate_tax(tax, level_int)
            if not t:
                continue
            bag[t] += float(count) / total
        out[sample_id] = dict(bag)
    return out


def _expand_group_label(tax_str):
    """Expand a taxonomy string with one or more `<rank>__[a|b|c]` group
    labels into the cartesian product of single-member taxonomy strings.
    A bracketless string returns a 1-element list. Used to resolve
    taxonova's species-group predictions against WGS truth."""
    expansions = [tax_str]
    for m in _GROUP_RE.finditer(tax_str):
        token = m.group(0)              # e.g. 's__[a|b|c]'
        prefix = m.group(1)             # 's__'
        members = [x for x in m.group(2).split("|") if x]
        if not members:
            continue
        new = []
        for e in expansions:
            for memb in members:
                new.append(e.replace(token, prefix + memb, 1))
        expansions = new
    return expansions


def _resolve_group_keys(amp_dict, wgs_dict):
    """Re-key `amp_dict` so taxonova group-label predictions
    (`<rank>__[a|b|c]`) are resolved against the WGS truth set.

    Mirrors the mock/biome evaluator's `evaluate_classification`
    semantics (one ASV = one prediction unit, comment-2.9):

    For an ASV predicting `<rank>__[a|b|c]` at frequency `f`:
      * `matches = {a,b,c} ∩ wgs_truth_set`
      * If `matches` non-empty → `f` is split equally among the
        matching members (full TP credit goes to the truth species).
      * If empty → the bracket is kept as **a single bucket** with
        frequency `f`. This is critical: splitting an unmatched
        bracket into N member buckets would inflate the false-positive
        count by N for that one prediction, unfairly punishing taxonova
        relative to single-name classifiers (mock per-pair PRF counts
        1 ASV = 1 unit regardless of bracket cardinality).

    Bracketless predictions (other classifiers + non-group taxonova
    models) pass through unchanged. WGS truth is bracket-free by
    construction (GG2 lineage strings)."""
    if not amp_dict:
        return amp_dict
    out = {}
    wgs_keys = set(wgs_dict)
    for k, freq in amp_dict.items():
        if "[" not in k:
            out[k] = out.get(k, 0.0) + freq
            continue
        expansions = _expand_group_label(k)
        if not expansions:
            continue
        matches = [e for e in expansions if e in wgs_keys]
        if matches:
            share = freq / len(matches)
            for e in matches:
                out[e] = out.get(e, 0.0) + share
        else:
            # Keep the bracket as one bucket — 1 ASV = 1 prediction unit.
            out[k] = out.get(k, 0.0) + freq
    return out


# --------------------------------------------------------------------------- #
# Metrics
# --------------------------------------------------------------------------- #

def _bray_curtis(a, b):
    """a, b : dict[taxon -> freq]. Returns 0..1 (0 = identical)."""
    keys = set(a) | set(b)
    if not keys:
        return 0.0
    num = sum(abs(a.get(k, 0.0) - b.get(k, 0.0)) for k in keys)
    den = sum(a.values()) + sum(b.values())
    return num / den if den else 0.0


def _jaccard(a, b):
    sa, sb = set(a), set(b)
    if not sa and not sb:
        return 1.0
    return len(sa & sb) / len(sa | sb)


# eps zero-replacement: small constant added to every entry in the union
# of species, after which the composition is re-normalised. eps=1e-6 is the
# de-facto default for compositional data with ~10^2..10^3 categories.
_AITCHISON_EPS = 1e-6


def _aitchison(a, b, eps=_AITCHISON_EPS):
    """Aitchison distance: Euclidean on CLR-transformed compositions.

    a, b : dict[taxon -> freq]. Computed over the union of keys; missing
    or zero entries are replaced with `eps` and the composition is
    re-normalised. Returns >= 0; 0 = identical compositions.

    This is compositional-aware (Aitchison 1986): a 0.001 vs 0.002 gap
    contributes log(2) to the distance just like a 0.4 vs 0.8 gap does,
    so abundance shifts are measured on the log-ratio scale rather than
    absolute. Wrong species predicted at high mass against zero-in-truth
    contribute large CLR differences (eps-replaced truth gives a big
    negative log term).
    """
    import math
    keys = set(a) | set(b)
    if not keys:
        return 0.0
    va = [a.get(k, 0.0) or eps for k in keys]
    vb = [b.get(k, 0.0) or eps for k in keys]
    # Re-normalise to closure on the simplex.
    sa, sb = sum(va), sum(vb)
    if sa <= 0 or sb <= 0:
        return 0.0
    va = [x / sa for x in va]
    vb = [x / sb for x in vb]
    # CLR: log(x_i) - mean(log(x_j))
    la = [math.log(x) for x in va]
    lb = [math.log(x) for x in vb]
    ma = sum(la) / len(la)
    mb = sum(lb) / len(lb)
    s = 0.0
    for li, lj in zip(la, lb):
        d = (li - ma) - (lj - mb)
        s += d * d
    return math.sqrt(s)


def _prf(a, b):
    """Treat WGS taxa (b) as truth, amplicon taxa (a) as prediction."""
    sa, sb = set(a), set(b)
    tp = len(sa & sb)
    fp = len(sa - sb)
    fn = len(sb - sa)
    p = tp / (tp + fp) if (tp + fp) else 0.0
    r = tp / (tp + fn) if (tp + fn) else 0.0
    f = 2 * p * r / (p + r) if (p + r) else 0.0
    return p, r, f


def _l1l2(a, b):
    keys = set(a) | set(b)
    l1 = sum(abs(a.get(k, 0.0) - b.get(k, 0.0)) for k in keys)
    l2 = sum((a.get(k, 0.0) - b.get(k, 0.0)) ** 2 for k in keys)
    return l1, l2


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #

def _per_sample_rows(amp_freqs, wgs_freqs, pairing,
                     study, classifier, param, level,
                     unifrac_pairs=None):
    """Per-paired-sample metric rows.

    `unifrac_pairs`: optional {(amp_id, wgs_id): float} from
    `_unifrac_metric.weighted_unifrac_pairs`. Only meaningful at species
    rank (tree tips are species-resolved). NaN at other ranks or if
    UniFrac was skipped via --no-unifrac.
    """
    rows = []
    for amp_id, wgs_id in pairing.items():
        a = amp_freqs.get(amp_id, {})
        b = wgs_freqs.get(wgs_id, {})
        if not a and not b:
            continue
        bc = _bray_curtis(a, b)
        jc = _jaccard(a, b)
        ait = _aitchison(a, b)
        p, r, f = _prf(a, b)
        l1, l2 = _l1l2(a, b)
        wuf = float("nan")
        if unifrac_pairs is not None:
            wuf = unifrac_pairs.get((amp_id, wgs_id), float("nan"))
        rows.append({
            "study": study, "classifier": classifier, "param": param,
            "level": level,
            "amplicon_sample": amp_id, "wgs_sample": wgs_id,
            "jaccard":    round(jc, 6),
            "aitchison":  round(ait, 6),
            "bray_curtis": round(bc, 6),
            "weighted_unifrac": (round(wuf, 6) if wuf == wuf else wuf),
            "precision":  round(p, 6),
            "recall":     round(r, 6),
            "f-score":    round(f, 6),
            "l1_error":   round(l1, 6),
            "l2_error":   round(l2, 6),
            "n_amplicon_taxa": len(a),
            "n_wgs_taxa":      len(b),
        })
    return rows


def _process_study(study, meta, rw_root, res_root, eval_root,
                   compute_unifrac=False):
    """Returns (long_rows, per_sample_files_written)."""
    if not all(meta.get(k) for k in ("amplicon_table", "wgs_ogu_biom",
                                     "wgs_taxonomy", "paired_samples")):
        print(f"  - {study}: missing wgs/amplicon paths in mocks.json, skipping eval")
        return [], []

    amp_tab_fp = rw_root / meta["amplicon_table"]
    wgs_biom_fp = rw_root / meta["wgs_ogu_biom"]
    wgs_tax_fp  = rw_root / meta["wgs_taxonomy"]
    paired_fp   = rw_root / meta["paired_samples"]
    for fp in (amp_tab_fp, wgs_biom_fp, wgs_tax_fp, paired_fp):
        if not fp.exists():
            print(f"  ! {study}: missing {fp}, skipping eval")
            return [], []

    paired_df = pd.read_csv(paired_fp)
    pairing = dict(zip(paired_df["amplicon_16s_run"], paired_df["shotgun_wgs_run"]))

    amp_count = _load_amplicon_table(amp_tab_fp)
    wgs_count = _load_wgs_ogu_table(wgs_biom_fp).T  # samples x OGUs
    wgs_tax = _load_taxonomy_tsv(wgs_tax_fp)

    # Pre-compute WGS per-sample freqs at all levels (study-wide, classifier-independent).
    wgs_freqs_at = {lvl: _per_sample_freqs(wgs_count, wgs_tax, LEVEL_TO_INT[lvl])
                    for lvl in LEVELS}

    long_rows = []
    per_sample_outs = []

    study_res_dir = res_root / study
    if not study_res_dir.is_dir():
        print(f"  - {study}: no predictions under {study_res_dir}")
        return [], []

    for clf_dir in sorted(study_res_dir.iterdir()):
        if not clf_dir.is_dir():
            continue
        classifier = clf_dir.name
        for param_dir in sorted(clf_dir.iterdir()):
            if not param_dir.is_dir():
                continue
            param = param_dir.name
            pred_fp = param_dir / "taxonomy_final.tsv"
            if not pred_fp.exists():
                continue
            try:
                asv_tax = _load_taxonomy_tsv(pred_fp)
            except Exception as exc:
                print(f"  ! {study}/{classifier}/{param}: load failed ({exc})")
                continue

            # NOTE on group labels (`s__[a|b|c]`): we keep them in
            # `asv_tax` and resolve them per-sample against WGS truth in
            # `_resolve_group_keys` below. The earlier "collapse to
            # first member" behaviour penalised taxonova whenever WGS
            # truth happened to be a non-first member, which is unfair
            # and inconsistent with the mock/biome evaluator's
            # comment-2.9 multi-rank set-membership matching.

            for level in LEVELS:
                amp_freqs = _per_sample_freqs(amp_count, asv_tax, LEVEL_TO_INT[level])
                wgs_freqs = wgs_freqs_at[level]
                # Truth-aware group-label resolution per amplicon sample.
                amp_freqs = {sid: _resolve_group_keys(d, wgs_freqs.get(pairing.get(sid, ""), {}))
                             for sid, d in amp_freqs.items()}
                # Weighted UniFrac: only meaningful at species rank
                # (tips are species-resolved in GG2 phylogeny). Slow
                # (~30 s per (study,classifier,param)) so opt-in via flag.
                unifrac_pairs = None
                if compute_unifrac and level == "Species":
                    try:
                        unifrac_pairs = _unifrac_metric.weighted_unifrac_pairs(
                            amp_freqs, wgs_freqs, pairing,
                        )
                    except Exception as exc:
                        print(f"  ! {study}/{classifier}/{param}: unifrac failed ({exc})")
                        unifrac_pairs = None
                rows = _per_sample_rows(amp_freqs, wgs_freqs, pairing,
                                        study, classifier, param, level,
                                        unifrac_pairs=unifrac_pairs)
                long_rows.extend(rows)

                # one per-sample CSV per (study x classifier x param)
                if level == "Species":
                    df_ps = pd.DataFrame(rows)
                    if not df_ps.empty:
                        ps_fp = (eval_root / "per_sample" /
                                 f"{study}_{classifier}_{param}.csv")
                        ps_fp.parent.mkdir(parents=True, exist_ok=True)
                        df_ps.to_csv(ps_fp, index=False)
                        per_sample_outs.append(str(ps_fp))

    return long_rows, per_sample_outs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--studies", nargs="+", default=None,
                    help="Restrict to specific studies (default: all in mocks.json).")
    ap.add_argument("--results-dir", default=None,
                    help="Default: results/real_world")
    ap.add_argument("--eval-dir", default=None,
                    help="Default: results/evaluation/realworld")
    ap.add_argument("--unifrac", action="store_true", default=False,
                    help="Compute weighted UniFrac at species rank (slow; "
                         "off by default since Aitchison + Jaccard are the "
                         "current promotion-gate metrics)")
    args = ap.parse_args()

    cfg = _common.load_pipeline_config()
    mocks_cfg = _common.load_mocks_config()

    rw_root = Path(cfg["real_world_input_dir"])
    res_root = (Path(args.results_dir) if args.results_dir
                else Path(cfg["results_dir"]) / "real_world")
    eval_root = (Path(args.eval_dir) if args.eval_dir
                 else Path(cfg["results_dir"]) / "evaluation" / "realworld")

    studies = args.studies or list(mocks_cfg["real_world"].keys())

    eval_root.mkdir(parents=True, exist_ok=True)
    (eval_root / "by_metric").mkdir(exist_ok=True)
    (eval_root / "per_sample").mkdir(exist_ok=True)

    all_rows = []
    for study in studies:
        meta = mocks_cfg["real_world"].get(study)
        if not meta:
            print(f"[05_realworld_evaluate] unknown study {study}; skip")
            continue
        print(f"[05_realworld_evaluate] === {study} ===")
        rows, _ = _process_study(study, meta, rw_root, res_root, eval_root,
                                 compute_unifrac=args.unifrac)
        all_rows.extend(rows)
        print(f"  + {len(rows)} per-sample rows")

    if not all_rows:
        print("[05_realworld_evaluate] no rows produced — nothing to write")
        return

    long_df = pd.DataFrame(all_rows)
    long_fp = eval_root / "all_metrics.csv"
    long_df.to_csv(long_fp, index=False)
    print(f"[05_realworld_evaluate] wrote {long_fp} ({len(long_df):,} rows)")

    # Aggregate (study × classifier × param × level) by mean of per-sample.
    metrics = ["jaccard", "aitchison", "bray_curtis",
               "precision", "recall",
               "f-score", "l1_error", "l2_error"]
    if args.unifrac:
        metrics.insert(2, "weighted_unifrac")
    for metric in metrics:
        pivot = long_df.pivot_table(
            index=["study", "classifier", "param"],
            columns="level", values=metric, aggfunc="mean",
        ).reindex(columns=LEVELS)
        fp = eval_root / "by_metric" / f"{metric}.csv"
        pivot.to_csv(fp)
        print(f"[05_realworld_evaluate] wrote {fp}")

    print(f"[05_realworld_evaluate] done -> {eval_root}")


if __name__ == "__main__":
    main()
