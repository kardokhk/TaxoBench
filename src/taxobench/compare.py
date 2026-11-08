"""`taxobench compare`: release A against release B.

Both runs must share the reference and have predicted the same queries in
every cell they cover (fail closed otherwise).

Per cell (tier x dataset x region x ref_id) x rank x group rule x metric:
A, B, delta = B - A, paired cluster bootstrap 95% CI of delta (the same
resampled clusters for both), bootstrap P, and an exact McNemar test on the
TP outcome. Per tier: the cells pooled (clusters keyed by cell), with the same
bootstrap; plus the median delta over cells, Wilcoxon signed-rank P and the
Hodges-Lehmann CI.

Confirmatory family (declared in a JSON plan before the run): entries
{"tier", "rank", "metric", "rule"}; their pooled-tier bootstrap P values are
Holm-adjusted. Every other P is Benjamini-Hochberg adjusted and labelled
exploratory. Guardrails: {"tier", "rank", "metric", "rule", "margin"}; a
guardrail passes when the lower CI bound of the pooled delta is above -margin
(non-inferiority).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from . import RANKS
from . import metrics as M
from . import stats as S
from .contract import ContractError
from .score import Context, _subset

METRICS = ["f1", "precision", "recall", "edgar_accuracy", "coverage", "mcr", "ucr", "ocr"]
#: counts carried on every row for A and B (sums over the rows compared)
SIDE_COUNTS = ["uc_wrong_above", "ca_wrong_above", "n_homonym_resolved", "n_group_unexpanded"]


def _side(MA, MB):
    i = {c: k for k, c in enumerate(M.COUNT_COLS)}
    out = {}
    for c in SIDE_COUNTS:
        out[f"{c}_A"] = float(MA[:, i[c]].sum())
        out[f"{c}_B"] = float(MB[:, i[c]].sum())
    return out


def load_plan(fp) -> dict:
    if fp is None:
        return {"confirmatory": [], "guardrails": []}
    with open(fp) as fh:
        plan = json.load(fh)
    for k in ("confirmatory", "guardrails"):
        plan.setdefault(k, [])
        for e in plan[k]:
            for f in ("tier", "rank", "metric"):
                if f not in e:
                    raise ContractError(f"plan {fp}: {k} entry {e} lacks {f!r}")
            e.setdefault("rule", "G1")
    return plan


def _find_run(ctx: Context, run_id: str):
    for run, preds in ctx.run_rows():
        if run["run_id"] == run_id:
            return run, preds
    raise ContractError(f"run {run_id!r} not found")


def compare(ctx: Context, run_a: str, run_b: str, plan=None, rules=("G1", "G2"),
            ranks=RANKS, metrics=METRICS, split="eval", n_boot=10000, seed=20261031):
    plan = plan or {"confirmatory": [], "guardrails": []}
    A, pa = _find_run(ctx, run_a)
    B, pb = _find_run(ctx, run_b)
    if A["ref_id"] != B["ref_id"]:
        raise ContractError(f"runs use different references: {A['ref_id']} vs {B['ref_id']}")
    cells_a = {tuple(c.values()): (al, ref) for _, c, al, ref in ctx.cells(A, pa)}
    cells_b = {tuple(c.values()): (al, ref) for _, c, al, ref in ctx.cells(B, pb)}
    if set(cells_a) != set(cells_b):
        raise ContractError(f"runs cover different cells: only A {sorted(set(cells_a) - set(cells_b))}, "
                            f"only B {sorted(set(cells_b) - set(cells_a))}")
    rows, pooled = [], {}
    for key in sorted(cells_a):
        ala, ref = cells_a[key]
        alb, _ = cells_b[key]
        if ala.query_ids != alb.query_ids:
            raise ContractError(f"cell {key}: the runs predicted different query sets")
        ma, mb = ala.split == split, alb.split == split
        if not ma.any():
            continue
        ala, alb = _subset(ala, ma), _subset(alb, mb)
        tier = key[0]
        for rank in ranks:
            for rule in rules:
                da, db = M.outcomes(ala, rank, rule), M.outcomes(alb, rank, rule)
                keep = (da["outcome"] != "excluded").to_numpy()
                ca, MA = M.cluster_counts(da[keep])
                cb, MB = M.cluster_counts(db[keep])
                assert (ca == cb).all()
                mc = S.mcnemar_exact((da.loc[keep, "outcome"] == "TP").to_numpy(),
                                     (db.loc[keep, "outcome"] == "TP").to_numpy())
                pk = (tier, rank, rule)
                pooled.setdefault(pk, ([], []))
                pooled[pk][0].append(MA)
                pooled[pk][1].append(MB)
                side = _side(MA, MB)
                for m in metrics:
                    r = S.paired_cluster_bootstrap(MA, MB, m, n_boot=n_boot, seed=seed)
                    rows.append(dict(level="cell", tier=tier, dataset=key[1], region=key[2],
                                     ref_id=key[3], rank=rank, rule=rule, metric=m,
                                     n=int(keep.sum()), **r, **mc, **side))
    cell_df = pd.DataFrame(rows)
    # pooled per tier, plus cell-level Wilcoxon / Hodges-Lehmann
    trows = []
    for (tier, rank, rule), (las, lbs) in sorted(pooled.items()):
        MA, MB = np.vstack(las), np.vstack(lbs)
        for m in metrics:
            r = S.paired_cluster_bootstrap(MA, MB, m, n_boot=n_boot, seed=seed)
            sub = cell_df[(cell_df.tier == tier) & (cell_df["rank"] == rank) &
                          (cell_df.rule == rule) & (cell_df.metric == m)]
            hl = S.hodges_lehmann(sub["delta"].to_numpy())
            trows.append(dict(level="tier", tier=tier, dataset="*", region="*", ref_id=A["ref_id"],
                              rank=rank, rule=rule, metric=m, n=int(MA[:, 0].sum() + MA[:, 1].sum()),
                              **r, **hl, **_side(MA, MB)))
    tier_df = pd.DataFrame(trows)
    out = pd.concat([tier_df, cell_df], ignore_index=True)
    out["family"] = "exploratory"
    conf_idx = []
    for e in plan["confirmatory"]:
        m = ((out.level == "tier") & (out.tier == e["tier"]) & (out["rank"] == e["rank"]) &
             (out.metric == e["metric"]) & (out.rule == e["rule"]))
        if not m.any():
            raise ContractError(f"confirmatory entry {e} matches no compared tier")
        conf_idx.extend(np.flatnonzero(m.to_numpy()))
    out["p_adj"] = np.nan
    if conf_idx:
        out.loc[conf_idx, "family"] = "confirmatory"
        out.loc[conf_idx, "p_adj"] = S.holm(out.loc[conf_idx, "p_boot"].to_numpy())
    ex = out.index[(out.family == "exploratory") & out["p_boot"].notna()]
    out.loc[ex, "p_adj"] = S.bh(out.loc[ex, "p_boot"].to_numpy())
    guard = []
    for g in plan["guardrails"]:
        m = ((out.level == "tier") & (out.tier == g["tier"]) & (out["rank"] == g["rank"]) &
             (out.metric == g["metric"]) & (out.rule == g["rule"]))
        if not m.any():
            raise ContractError(f"guardrail {g} matches no compared tier")
        r = out[m].iloc[0]
        guard.append(dict(g, delta=r["delta"], ci_low=r["ci_low"],
                          passed=bool(r["ci_low"] > -float(g["margin"]))))
    meta = {"run_a": run_a, "run_b": run_b, "version_a": A["version"], "version_b": B["version"],
            "ref_id": A["ref_id"], "split": split, "n_boot": n_boot, "seed": seed,
            "plan": plan}
    return out, pd.DataFrame(guard), meta


def _f(x, d=3):
    return "NA" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.{d}f}"


def _p(x):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "NA"
    return "< 0.001" if x < 0.001 else f"{x:.2f}" if x >= 0.01 else f"{x:.3f}"


def to_markdown(out: pd.DataFrame, guard: pd.DataFrame, meta: dict) -> str:
    L = [f"# Comparison: {meta['run_b']} (B) against {meta['run_a']} (A)", "",
         f"Reference `{meta['ref_id']}`, split `{meta['split']}`. Delta = B - A. 95% CI from a "
         f"paired cluster bootstrap ({meta['n_boot']} resamples, seed {meta['seed']}). "
         "Confirmatory P values are Holm-adjusted over the declared family; all others are "
         "Benjamini-Hochberg adjusted and exploratory. A query not named at a rank but wrong "
         "higher up counts as UC there (`uc_wrong_above`), so precision judges only calls named "
         "at the rank.", ""]
    conf = out[out.family == "confirmatory"]
    if len(conf):
        L += ["## Confirmatory family", "", "| tier | rank | rule | metric | A | B | delta (95% CI) | P (Holm) |",
              "|---|---|---|---|---|---|---|---|"]
        for _, r in conf.iterrows():
            L.append(f"| {r.tier} | {r['rank']} | {r.rule} | {r.metric} | {_f(r.A)} | {_f(r.B)} | "
                     f"{_f(r.delta)} ({_f(r.ci_low)} to {_f(r.ci_high)}) | {_p(r.p_adj)} |")
        L.append("")
    if len(guard):
        L += ["## Guardrails (non-inferiority)", "", "| tier | rank | rule | metric | margin | delta | CI low | result |",
              "|---|---|---|---|---|---|---|---|"]
        for _, g in guard.iterrows():
            L.append(f"| {g.tier} | {g['rank']} | {g.rule} | {g.metric} | {g.margin} | {_f(g.delta)} | "
                     f"{_f(g.ci_low)} | {'pass' if g.passed else 'FAIL'} |")
        L.append("")
    tier = out[(out.level == "tier") & (out.metric.isin(["f1", "precision", "recall", "coverage"]))]
    L += ["## Per tier (cells pooled)", "",
          "| tier | rank | rule | metric | A | B | delta (95% CI) | P (BH) | cells | HL (95% CI) | UC wrong above A/B | homonym-resolved A/B | unexpanded groups A/B |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in tier.iterrows():
        L.append(f"| {r.tier} | {r['rank']} | {r.rule} | {r.metric} | {_f(r.A)} | {_f(r.B)} | "
                 f"{_f(r.delta)} ({_f(r.ci_low)} to {_f(r.ci_high)}) | {_p(r.p_adj)} | "
                 f"{int(r.n_cells) if r.n_cells == r.n_cells else 'NA'} | "
                 f"{_f(r.hl)} ({_f(r.hl_low)} to {_f(r.hl_high)}) | "
                 f"{int(r.uc_wrong_above_A)}/{int(r.uc_wrong_above_B)} | "
                 f"{int(r.n_homonym_resolved_A)}/{int(r.n_homonym_resolved_B)} | "
                 f"{int(r.n_group_unexpanded_A)}/{int(r.n_group_unexpanded_B)} |")
    L += ["", "Per-cell rows, all metrics and McNemar counts are in the TSV.", ""]
    return "\n".join(L)


def write(out, guard, meta, out_dir):
    d = Path(out_dir)
    d.mkdir(parents=True, exist_ok=True)
    out.to_csv(d / "compare.tsv", sep="\t", index=False, na_rep="NA", float_format="%.10g")
    guard.to_csv(d / "guardrails.tsv", sep="\t", index=False, na_rep="NA")
    (d / "compare.md").write_text(to_markdown(out, guard, meta))
    (d / "compare.json").write_text(json.dumps(meta, indent=1, default=str))
    return d
