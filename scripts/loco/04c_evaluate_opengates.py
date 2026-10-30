#!/usr/bin/env python3
"""LOCO step 04c: MC/OC sweep for taxonova from unmasked ("open gates") output.

Why this exists
---------------
`04_evaluate.py` reads the *processed* prediction TSVs, in which abstention has
already been applied and propagated downward at classification time. Taxonova's
parent-rank gates were pinned at 0.5 there, so sweeping
`--confidence_threshold` moved only the species gate. That put a floor on the
reachable under-classification rate and truncated taxonova's MC/OC curves.

This script instead reads the raw per-rank output produced with every gate at
0.0, then applies one threshold uniformly to every rank and propagates
abstention downward in software, once per point on the sweep. That is what makes
the curve comparable to a competitor's.

`04_evaluate.py` is deliberately left untouched: it produced the published
numbers and must keep reproducing them.

Usage
-----
    python scripts/loco/04c_evaluate_opengates.py \
        --regions V4 V3-V4 V5-V7 V1-V9 RRN \
        --out results/loco/evaluation/mc_oc_curve.opengates.csv
"""
import argparse
import csv
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# 04_evaluate.py is not importable by name (leading digit), so load it by path
# and reuse its parsing, truth loading and rate arithmetic verbatim.
_spec = importlib.util.spec_from_file_location(
    "loco_eval", Path(__file__).with_name("04_evaluate.py"))
_ev = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_ev)

LEVELS = _ev.LEVELS
LEVEL_TO_INT = _ev.LEVEL_TO_INT
RANKS_OF_INTEREST = _ev.RANKS_OF_INTEREST

# Finer than the default 0.05 grid near the ends, where the curves bend hardest.
THRESHOLDS = sorted({round(v, 4) for v in (
    [0.0] + [0.01 * i for i in range(1, 10)]
    + [0.05 * i for i in range(2, 20)]
    + [0.95, 0.97, 0.99, 0.995, 0.999, 1.0]
)})


def load_opengates_raw(path):
    """-> {id: (labels[7], confs[7])} straight from taxonova's raw TSV.

    No thresholding and no propagation: both are applied later, per tau. A rank
    whose label is an abstention token is recorded as empty with confidence 0,
    because the model declined there regardless of any threshold.
    """
    out = {}
    with open(path) as fh:
        # QUOTE_NONE for symmetry with the truth reader; unexposed on
        # Greengenes2, which carries no quoted names, but the same
        # asymmetry silently corrupted the TAXXI evaluator.
        for row in csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE):
            labels, confs = [], []
            for lvl in LEVELS:
                name = (row.get(lvl) or "").strip()
                try:
                    conf = float(row.get(f"{lvl}_confidence") or 0.0)
                except ValueError:
                    conf = 0.0
                if name.lower() in ("unclassified", "unknown", "unassigned", ""):
                    name, conf = "", 0.0
                labels.append(name)
                confs.append(conf)
            out[row["sequence_id"]] = (labels, confs)
    return out


# Threshold policies. Each maps a sweep parameter tau to a 7-vector of per-rank
# gates, so one curve is still one parameter but the ranks need not share it.
#
#   uniform  every rank gated at tau. The like-for-like comparison:
#            no rank can be gated differently from any other, so the
#            curve is directly comparable to a competitor's.
#   native   what taxonova actually ships. `strict_thresholds_default` in the
#            V4 checkpoint's ood_config is {domain..genus: 0.5, species: 0.7},
#            matching cli/config_builder.py. Only species moves with tau.
#   graded   the manuscript's recommendation (genus 0.7, species 0.8, parents
#            0.5) swept proportionally, so tau scales the whole profile instead
#            of only its deepest rank.
#   ranklocal gate only the rank being scored, at tau, with no propagation.
#            This is how every competitor is evaluated by 04_evaluate.py, so it
#            is the most directly comparable, but it ignores taxonova's
#            hierarchy: it will credit a species call whose genus was never
#            made. Report it alongside `uniform`, not instead of it.
NATIVE_PARENT_GATE = 0.5
GRADED_PROFILE = [0.5, 0.5, 0.5, 0.5, 0.5, 0.7, 0.8]
POLICIES = ("uniform", "native", "graded", "ranklocal")


def gates_for(policy, tau):
    if policy == "ranklocal":
        raise ValueError("ranklocal is applied per scored rank, not as a vector")
    if policy == "uniform":
        return [tau] * 7
    if policy == "native":
        return [NATIVE_PARENT_GATE] * 6 + [tau]
    if policy == "graded":
        scale = tau / GRADED_PROFILE[-1]
        return [min(1.0, g * scale) for g in GRADED_PROFILE]
    raise ValueError(f"unknown policy {policy!r}")


def apply_gates(labels, confs, gates):
    """Gate each rank at its own threshold and propagate abstention downward.

    Mirrors src/process_taxonova_output.py: once a rank is dropped, every rank
    below it is dropped too, regardless of its own confidence. That propagation
    is the whole point of the downward-traversal design, and it is why a parent
    gate set too high silently destroys correct species calls.
    """
    out_labels, out_confs = [], []
    dropped = False
    for name, conf, gate in zip(labels, confs, gates):
        if dropped or not name or conf < gate:
            dropped = True
            out_labels.append("")
            out_confs.append(0.0)
        else:
            out_labels.append(name)
            out_confs.append(conf)
    return out_labels, out_confs


def evaluate_region(region, control_raw, novel_raw, truth_paths, policy):
    novel_truth = _ev.load_truth(truth_paths["novel_truth"])
    control_truth = _ev.load_truth(truth_paths["control_truth"])
    novel_meta = _ev.load_novel_meta(truth_paths["novel_meta"])

    rows = []
    for rank in RANKS_OF_INTEREST:
        L = LEVEL_TO_INT[rank]
        mc_ids = [i for i, t in control_truth.items() if not _ev._is_empty(t[L])]
        oc_ids = [i for i, lvl in novel_meta.items() if lvl <= L]

        for tau in THRESHOLDS:
            if policy == "ranklocal":
                def called(labels, confs):
                    return bool(labels[L]) and confs[L] >= tau
            else:
                gates = gates_for(policy, tau)

                def called(labels, confs, _g=gates):
                    pp, _ = apply_gates(labels, confs, _g)
                    return not _ev._is_empty(pp[L])

            mc_err = mc_abstain = 0
            for _id in mc_ids:
                labels, confs = control_raw.get(_id, ([""] * 7, [0.0] * 7))
                if not called(labels, confs):
                    mc_abstain += 1
                    continue
                if (_ev._strip_prefix(control_truth[_id][L]).lower()
                        != _ev._strip_prefix(labels[L]).lower()):
                    mc_err += 1

            oc_err = oc_abstain = 0
            for _id in oc_ids:
                labels, confs = novel_raw.get(_id, ([""] * 7, [0.0] * 7))
                if not called(labels, confs):
                    oc_abstain += 1
                else:
                    oc_err += 1

            dmc, doc = max(1, len(mc_ids)), max(1, len(oc_ids))
            rows.append({
                "region": region, "classifier": "taxonova",
                "policy": policy,
                "rank": rank, "threshold": tau,
                "mc_rate": round(mc_err / dmc, 6),
                "oc_rate": round(oc_err / doc, 6),
                "mc_abstain_rate": round(mc_abstain / dmc, 6),
                "oc_abstain_rate": round(oc_abstain / doc, 6),
                "n_mc_total": len(mc_ids), "n_mc_errors": mc_err,
                "n_mc_abstain": mc_abstain,
                "n_oc_total": len(oc_ids), "n_oc_errors": oc_err,
                "n_oc_abstain": oc_abstain,
            })
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regions", nargs="+",
                    default=["V4", "V3-V4", "V5-V7", "V1-V9", "RRN"])
    ap.add_argument("--loco-root", default="loco",
                    help="Subdirectory of data/ and results/ holding the run. "
                         "Default 'loco' is the shipped 200-species holdout; "
                         "the 2,000-species replicate is 'loco2k'. "
                         "--results-root/--data-root/--out still override it.")
    ap.add_argument("--tag", default=None,
                    help="Suffix on the region directory, e.g. 'seed2026' for "
                         "'<region>__seed2026'. The region column stays bare.")
    ap.add_argument("--results-root", default=None)
    ap.add_argument("--data-root", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--policies", nargs="+", default=list(POLICIES),
                    choices=list(POLICIES))
    ap.add_argument("--extra-thresholds", nargs="*", type=float, default=[],
                    help="Additional sweep points, merged into the default grid "
                         "and de-duplicated. Purely additive: with none given "
                         "(the default) the grid and the output are unchanged, "
                         "so every published curve reproduces byte-identically. "
                         "Use it where the default 0.05 grid is too coarse to "
                         "reach a target control under-classification, e.g. "
                         "--extra-thresholds 0.9992 0.9995 0.9999 for the RRN "
                         "species curve, whose confidences saturate near 1.")
    ap.add_argument("--arm", default="taxonova",
                    help="Which model directory under results/loco/<region>/ to "
                         "score, e.g. 'taxonova_flat' or 'taxonova_flat584' for "
                         "the M4 ablation. The label written to the classifier "
                         "column follows the directory name.")
    args = ap.parse_args()

    # Additive only. `THRESHOLDS` is a module global read by evaluate_region, so
    # the extra points are merged in place rather than threaded through every
    # call site. An empty --extra-thresholds leaves the grid exactly as shipped.
    if args.extra_thresholds:
        global THRESHOLDS
        THRESHOLDS = sorted({round(v, 8) for v in
                             list(THRESHOLDS) + list(args.extra_thresholds)})
        print(f"[04c] grid extended to {len(THRESHOLDS)} thresholds")

    # --loco-root supplies the defaults; the three explicit roots still win, so
    # existing call sites keep behaving exactly as before.
    results_root = Path(args.results_root or (ROOT / "results" / args.loco_root))
    data_root = Path(args.data_root or (ROOT / "data" / args.loco_root))
    out = Path(args.out or (results_root / "evaluation"
                            / "mc_oc_curve.opengates.csv"))
    all_rows, skipped = [], []

    for region in args.regions:
        # Directory carries the tag; the region column does not, so the output
        # stays joinable with 04_evaluate.py and with the figure scripts.
        region_dir = "{}__{}".format(region, args.tag) if args.tag else region
        pdir = results_root / region_dir / args.arm
        control = pdir / "control_predictions.opengates.raw.tsv"
        novel = pdir / "novel_predictions.opengates.raw.tsv"
        if not control.exists() or not novel.exists():
            skipped.append(region_dir)
            continue
        truth_paths = {
            "novel_truth":   data_root / region_dir / "novel_truth.tsv",
            "control_truth": data_root / region_dir / "control_truth.tsv",
            "novel_meta":    data_root / region_dir / "novel_meta.tsv",
        }
        # Load once, sweep every policy over the same in-memory predictions.
        c_raw, n_raw = load_opengates_raw(control), load_opengates_raw(novel)
        for policy in args.policies:
            rows = evaluate_region(region, c_raw, n_raw, truth_paths, policy)
            for r in rows:
                r["classifier"] = args.arm
            all_rows.extend(rows)
            print(f"[04c] {region}/{policy}: {len(rows)} rows")

    if skipped:
        # Never let a missing region look like a region with no error.
        print(f"[04c] SKIPPED (no opengates output): {', '.join(skipped)}",
              file=sys.stderr)
    if not all_rows:
        sys.exit("[04c] nothing evaluated")

    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(all_rows[0].keys()))
        w.writeheader()
        w.writerows(all_rows)
    print(f"[04c] wrote {out} ({len(all_rows)} rows)")


if __name__ == "__main__":
    main()
