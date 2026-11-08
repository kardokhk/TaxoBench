"""`taxobench report`: scores.tsv -> markdown summary.

Every row shows n, the 95% CI, the group-containing (G1) and strict (G2)
values side by side, and coverage, so a gain bought by grouping or abstaining
is visible.
"""
from __future__ import annotations

import math

import pandas as pd


def _f(x, d=3):
    try:
        x = float(x)
    except (TypeError, ValueError):
        return "NA"
    return "NA" if math.isnan(x) else f"{x:.{d}f}"


def _ci(r):
    return f"{_f(r['value'])} ({_f(r.get('ci_low'))} to {_f(r.get('ci_high'))})"


def to_markdown(scores: pd.DataFrame, ranks=("genus", "species")) -> str:
    s = scores.copy()
    keys = ["tier", "dataset", "region", "ref_id", "rank", "classifier", "version", "run_id"]
    L = ["# TaxoBench scores", "",
         "Rules v2. G1 = group-containing (default), G2 = strict. 95% CI from a cluster "
         "bootstrap. K = queries with a named truth at the rank, L = novel at the rank.", ""]
    for (tier, ds, region, ref_id), cell in s.groupby(["tier", "dataset", "region", "ref_id"], sort=True):
        L += [f"## {tier} / {ds} / {region} / {ref_id}", "",
              "| rank | classifier | version | K | L | F1 G1 (95% CI) | F1 G2 | precision | recall | coverage | OCR | oov share | UC wrong above | homonym-resolved | unexpanded groups |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for (rank, clf, ver, rid), g in cell.groupby(["rank", "classifier", "version", "run_id"], sort=True):
            if rank not in ranks:
                continue
            def get(metric, rule="G1"):
                m = g[(g.metric == metric) & (g.rule == rule)]
                return m.iloc[0] if len(m) else None
            def cnt(name):
                r = get("count_" + name)
                return "NA" if r is None else str(int(r["value"]))
            f1 = get("f1")
            if f1 is None:
                continue
            L.append(f"| {rank} | {clf} | {ver} | {int(f1['n_K'])} | {int(f1['n_L'])} | {_ci(f1)} | "
                     f"{_f(get('f1', 'G2')['value']) if get('f1', 'G2') is not None else 'NA'} | "
                     f"{_f(get('precision')['value'])} | {_f(get('recall')['value'])} | "
                     f"{_f(get('coverage')['value'])} | {_f(get('ocr')['value'])} | "
                     f"{_f(get('oov_share')['value'])} | {cnt('uc_wrong_above')} | "
                     f"{cnt('n_homonym_resolved')} | {cnt('n_group_unexpanded')} |")
        L.append("")
    return "\n".join(L)
