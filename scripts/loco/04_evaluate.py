#!/usr/bin/env python3
"""LOCO step 04 — compute MC and OC rates as a function of confidence threshold.

For each (region × classifier × rank × threshold τ):

    MC denominator = number of CONTROL queries (true rank-L taxon present in
                     retained reference)
    MC numerator   = controls with predicted rank-L taxon ≠ true rank-L taxon
                     AND confidence at L ≥ τ

    OC denominator = number of NOVEL queries with first_novel_rank ≤ L
                     (i.e. their true rank-L taxon is absent from retained ref)
    OC numerator   = those queries that received ANY rank-L call with
                     confidence at L ≥ τ

Outputs:
    results/loco/evaluation/mc_oc_curve.csv      # tidy long format
    results/loco/evaluation/mc_oc_default.csv    # at the per-classifier default τ

A "no call" at rank L can come from (a) confidence < τ, or (b) the classifier
returned a shorter lineage (e.g. only domain-phylum-class, no order). Either
counts as "abstained" at rank L.
"""
import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import _common  # noqa


LEVELS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
LEVEL_TO_INT = {n: i for i, n in enumerate(LEVELS)}
RANKS_OF_INTEREST = ["class", "order", "family", "genus", "species"]

# Fallback grid used when no per-classifier grid is loaded from configs/tau_grids.json.
THRESHOLDS = [round(0.05 * i, 2) for i in range(0, 21)]  # 0.0 .. 1.0 step 0.05

# Default operating thresholds quoted from each tool's documentation / typical use.
DEFAULT_TAU = {
    "taxonova":    0.7,
    "idtaxa":      0.6,    # 60 / 100
    "naive-bayes": 0.7,
    "blast":       0.51,   # consensus threshold
    "vsearch":     0.51,
    "spingo":      0.5,
    "kraken2":     0.0,    # kraken2's own --confidence default is 0
    "mapseq":      0.4,    # the cutoff the mock panel used
}


def _load_tau_grids(path):
    """Load configs/tau_grids.json; return {classifier: [tau,...]} for 'loco' kind.
    Falls back to THRESHOLDS for classifiers not in the file."""
    if not path or not Path(path).exists():
        return {}
    with open(path) as fh:
        doc = json.load(fh)
    out = {}
    for cls, entry in doc.get("grids", {}).items():
        grid = entry.get("loco") or []
        if grid:
            out[cls] = sorted({round(float(t), 4) for t in grid})
    return out


def _parse_lineage(tax_str, n_levels=7):
    """Return list[n_levels] of stripped lineage parts, padded with ''."""
    parts = [p.strip() for p in tax_str.split(";")]
    while len(parts) < n_levels:
        parts.append("")
    return parts[:n_levels]


def _strip_prefix(name):
    """'g__Bacillus' -> 'Bacillus'; 'Bacillus' -> 'Bacillus'."""
    return name.split("__", 1)[-1].strip() if name else ""


#: Every way a classifier in this panel says "no call at this rank". MAPseq
#: writes a literal "N/A"; without it here that is scored as a taxon named N/A,
#: i.e. a guaranteed error, which understates MAPseq. Adding it changes no
#: published number: no output from the six previously evaluated classifiers
#: contains the string (checked 2026-08-19).
#: Deliberately NOT here: "root". IDTAXA emits it 214,377 times as the leading
#: element of its lineage, which load_predictions already strips by position.
#: Treating it as an abstention token would reclassify real IDTAXA calls and
#: move published numbers. Nor "na"/"none", which no output contains and which
#: could collide with a genuine taxon name.
ABSTENTION_TOKENS = ("unclassified", "unknown", "unassigned", "n/a",
                     "no blast hit", "ambiguous")


def _is_empty(name):
    """True if a lineage cell is missing (e.g. 'g__' or '')."""
    if not name:
        return True
    s = _strip_prefix(name)
    return s == "" or s.lower() in ABSTENTION_TOKENS


def _parse_confidence_field(field, n_levels=7):
    """Parse a confidence column. Cases:
        '0.99'                               -> [0.99] * n_levels (single conf)
        '0.99; 0.97; 0.80; 0.5'              -> [...] padded with 0
        '0.5; 0.5; ... '                     -> as-is
    Returns list of length n_levels.
    """
    if field is None or field == "":
        return [0.0] * n_levels
    s = str(field).strip()
    if ";" in s:
        vals = []
        for tok in s.split(";"):
            tok = tok.strip()
            try:
                vals.append(float(tok))
            except ValueError:
                vals.append(0.0)
        # pad / trim to n_levels
        if len(vals) < n_levels:
            vals = vals + [0.0] * (n_levels - len(vals))
        return vals[:n_levels]
    try:
        v = float(s)
        return [v] * n_levels
    except ValueError:
        return [0.0] * n_levels


def _normalize_idtaxa_lineage(lineage_parts):
    """IDTAXA's output starts with 'Root;Domain;...'. Strip 'Root' and use raw
    names as-is (they don't carry d__/p__/etc. prefixes by default after our
    R training script). We keep them aligned with the 7-level scheme.
    """
    parts = [p for p in lineage_parts if p and p.lower() != "root"]
    while len(parts) < 7:
        parts.append("")
    return parts[:7]


def load_predictions(tsv_fp, classifier):
    """Yield (id, lineage[7], conf[7])."""
    with open(tsv_fp) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        # Some tools name the columns differently — find the right ones.
        cols = reader.fieldnames or []
        id_col   = next((c for c in cols if c in
                        ("Feature ID", "Feature.ID", "feature_id", "id", "sequence_id")), cols[0])
        tax_col  = next((c for c in cols if c.lower() in ("taxon", "taxonomy")), None)
        conf_col = next((c for c in cols if c.lower() in ("confidence", "consensus")), None)
        for row in reader:
            _id = row[id_col]
            tax = row.get(tax_col, "") if tax_col else ""
            conf = row.get(conf_col, "") if conf_col else ""
            if classifier == "idtaxa":
                # IDTAXA lineage has 8 elements (Root + 7 ranks). Parse as 8,
                # then drop Root so we end up with 7 = [D,P,C,O,F,G,S].
                parts8 = _parse_lineage(tax, n_levels=8)
                parts = parts8[1:]
                while len(parts) < 7:
                    parts.append("")
                confs = _parse_confidence_field(conf, n_levels=8)[1:]
            else:
                parts = _parse_lineage(tax, n_levels=7)
                confs = _parse_confidence_field(conf, n_levels=7)
            yield _id, parts, confs


def load_truth(tsv_fp):
    out = {}
    with open(tsv_fp) as fh:
        for i, line in enumerate(fh):
            parts = line.rstrip("\n").split("\t")
            if i == 0 and parts[0] in ("Feature ID", "id"):
                continue
            if len(parts) >= 2:
                out[parts[0]] = _parse_lineage(parts[1])
    return out


def load_novel_meta(tsv_fp):
    """Returns {id: novel_rank_int} — first rank where this query is novel."""
    out = {}
    with open(tsv_fp) as fh:
        for i, line in enumerate(fh):
            parts = line.rstrip("\n").split("\t")
            if i == 0 and parts[0] in ("Feature ID", "id"):
                continue
            if len(parts) >= 3:
                out[parts[0]] = int(parts[2])
    return out


def lineage_match(true_parts, pred_parts, level_int):
    """Are true and predicted lineages identical at all ranks 0..level_int?

    Compares stripped (no d__/p__/...) names case-insensitively. Empty cells
    on either side are treated as 'no name', not a match (so an empty pred
    cell counts as 'no call', handled separately upstream).
    """
    for lvl in range(level_int + 1):
        tn = _strip_prefix(true_parts[lvl]).lower()
        pn = _strip_prefix(pred_parts[lvl]).lower()
        if tn != pn or _is_empty(true_parts[lvl]) or _is_empty(pred_parts[lvl]):
            return False
    return True


def evaluate_one(region, classifier, predictions_dir, truth_paths, thresholds):
    """Return list of rows: (region, classifier, rank, threshold, MC, OC, n_mc, n_oc, n_abstain_mc, n_abstain_oc)."""
    novel_truth   = load_truth(truth_paths["novel_truth"])
    control_truth = load_truth(truth_paths["control_truth"])
    novel_meta    = load_novel_meta(truth_paths["novel_meta"])

    novel_pred = {}
    for _id, parts, confs in load_predictions(
            predictions_dir / "novel_predictions.tsv", classifier):
        novel_pred[_id] = (parts, confs)

    control_pred = {}
    for _id, parts, confs in load_predictions(
            predictions_dir / "control_predictions.tsv", classifier):
        control_pred[_id] = (parts, confs)

    rows = []
    for rank in RANKS_OF_INTEREST:
        L = LEVEL_TO_INT[rank]
        # --- MC denominator = controls whose true rank-L taxon is non-empty ---
        mc_ids = [_id for _id, true_parts in control_truth.items()
                  if not _is_empty(true_parts[L])]
        # --- OC denominator = novel queries with first_novel_rank <= L ---
        oc_ids = [_id for _id, novel_lvl in novel_meta.items()
                  if novel_lvl <= L]

        for tau in thresholds:
            mc_err = 0; mc_abstain = 0
            for _id in mc_ids:
                pp, cc = control_pred.get(_id, ([""] * 7, [0.0] * 7))
                if _is_empty(pp[L]) or cc[L] < tau:
                    mc_abstain += 1
                    continue
                # has a call at L; is it correct?
                tn = _strip_prefix(control_truth[_id][L]).lower()
                pn = _strip_prefix(pp[L]).lower()
                if tn != pn:
                    mc_err += 1
            oc_err = 0; oc_abstain = 0
            for _id in oc_ids:
                pp, cc = novel_pred.get(_id, ([""] * 7, [0.0] * 7))
                if _is_empty(pp[L]) or cc[L] < tau:
                    oc_abstain += 1
                    continue
                # ANY call at L is OC because true L-taxon is absent
                oc_err += 1
            denom_mc = max(1, len(mc_ids))
            denom_oc = max(1, len(oc_ids))
            rows.append({
                "region": region, "classifier": classifier,
                "rank": rank, "threshold": tau,
                "mc_rate": round(mc_err / denom_mc, 6),
                "oc_rate": round(oc_err / denom_oc, 6),
                "mc_abstain_rate": round(mc_abstain / denom_mc, 6),
                "oc_abstain_rate": round(oc_abstain / denom_oc, 6),
                "n_mc_total": len(mc_ids),
                "n_mc_errors": mc_err,
                "n_mc_abstain": mc_abstain,
                "n_oc_total": len(oc_ids),
                "n_oc_errors": oc_err,
                "n_oc_abstain": oc_abstain,
            })
    return rows


def _region_from_dir(dir_name, tag):
    """Bare region name for a (possibly tagged) region directory.

    Replicates coexist as '<region>__<tag>' (see 03_classify.py::_region_dir).
    The tag is stripped here so the 'region' column keeps reading 'V4' rather
    than 'V4__seed2026', which is what 05_plot.py and the figure scripts match
    on. Provenance is not lost: each replicate writes into its own
    results/<loco-root>/evaluation/ directory.
    """
    if tag and dir_name.endswith("__" + tag):
        return dir_name[: -(len(tag) + 2)]
    return dir_name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regions", nargs="+", default=None,
                    help="Bare region names, e.g. V4. Compared after the --tag "
                         "suffix is stripped from the directory name.")
    ap.add_argument("--classifiers", nargs="+", default=None)
    ap.add_argument("--tau-grid", default=str(ROOT / "configs" / "tau_grids.json"),
                    help="JSON with per-classifier LOCO tau grids.")
    ap.add_argument("--loco-root", default="loco",
                    help="Subdirectory of data/ and results/ holding the run. "
                         "Default 'loco' is the shipped 200-species holdout; "
                         "the 2,000-species replicate is 'loco2k'.")
    ap.add_argument("--tag", default=None,
                    help="Only evaluate region directories carrying this "
                         "suffix, e.g. 'seed2026' for '<region>__seed2026'.")
    args = ap.parse_args()
    per_cls_grids = _load_tau_grids(args.tau_grid)
    if per_cls_grids:
        print(f"[loco/04] loaded per-classifier tau grids for: {sorted(per_cls_grids)}")

    cfg = _common.load_pipeline_config()
    res_root  = Path(cfg["results_dir"]) / args.loco_root
    data_root = Path(cfg["data_dir"]) / args.loco_root
    out_root  = res_root / "evaluation"
    if not res_root.is_dir():
        sys.exit("missing {}; run 03_classify.py --loco-root {} first".format(
            res_root, args.loco_root))
    out_root.mkdir(parents=True, exist_ok=True)

    rows = []
    wrong_tag = []
    for region_dir in sorted(res_root.iterdir()):
        if not region_dir.is_dir() or region_dir.name == "evaluation":
            continue
        # A tagged run must not silently pick up an untagged sibling sitting in
        # the same root, and a mistyped tag must not read as "no results".
        if args.tag and not region_dir.name.endswith("__" + args.tag):
            wrong_tag.append(region_dir.name)
            continue
        region = _region_from_dir(region_dir.name, args.tag)
        if args.regions and region not in args.regions:
            continue
        truth_paths = {
            "novel_truth":   data_root / region_dir.name / "novel_truth.tsv",
            "control_truth": data_root / region_dir.name / "control_truth.tsv",
            "novel_meta":    data_root / region_dir.name / "novel_meta.tsv",
        }
        for cls_dir in sorted(region_dir.iterdir()):
            if not cls_dir.is_dir():
                continue
            classifier = cls_dir.name
            if args.classifiers and classifier not in args.classifiers:
                continue
            if not (cls_dir / "novel_predictions.tsv").exists():
                continue
            if not (cls_dir / "control_predictions.tsv").exists():
                continue
            thresholds = per_cls_grids.get(classifier, THRESHOLDS)
            print(f"[loco/04] evaluating {region}/{classifier} (n_tau={len(thresholds)})",
                  flush=True)
            rows.extend(evaluate_one(region, classifier, cls_dir, truth_paths, thresholds))

    if wrong_tag:
        print("[loco/04] skipped (tag != {!r}): {}".format(
            args.tag, ", ".join(sorted(wrong_tag))), file=sys.stderr)
    if not rows:
        print("[loco/04] no predictions found")
        return

    long_fp = out_root / "mc_oc_curve.csv"
    cols = ["region", "classifier", "rank", "threshold",
            "mc_rate", "oc_rate", "mc_abstain_rate", "oc_abstain_rate",
            "n_mc_total", "n_mc_errors", "n_mc_abstain",
            "n_oc_total", "n_oc_errors", "n_oc_abstain"]
    with open(long_fp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    print(f"[loco/04] wrote {long_fp} ({len(rows)} rows)")

    # Default-tau slice
    short_fp = out_root / "mc_oc_default.csv"
    with open(short_fp, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            tau = DEFAULT_TAU.get(r["classifier"])
            if tau is None:
                continue
            if abs(r["threshold"] - tau) < 0.025:  # closest grid point
                w.writerow(r)
    print(f"[loco/04] wrote {short_fp}")


if __name__ == "__main__":
    main()
