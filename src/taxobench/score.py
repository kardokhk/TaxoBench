"""`taxobench score`: tier tables + run tables -> long score table.

A cell is (tier, dataset, region, ref_id). A run scores the cells of its
ref_id whose queries it predicted; a cell it covers only partly is an error
(fail closed), a cell it does not touch is skipped and listed.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from . import RANKS, __version__
from . import metrics as M
from .adapters.base import REGISTRY, get_adapter
from .contract import ContractError, RunTables, TierTables
from .taxonomy import Reference

HEADLINE = ["precision", "recall", "f1", "edgar_accuracy", "mcr", "ucr", "ocr", "coverage"]
SPLIT_DEFAULT = "eval"


def confidence_kind(run: pd.Series) -> str:
    k = run.get("confidence_kind", "") if hasattr(run, "get") else ""
    if k:
        return k
    for col in ("adapter", "classifier"):
        name = run.get(col, "")
        if name:
            try:
                return get_adapter(name).confidence_kind
            except KeyError:
                pass
    return "unknown"


@dataclass
class Context:
    tiers: list
    runs: list
    refs: dict = field(default_factory=dict)

    @classmethod
    def load(cls, tier_dirs, run_dirs, verify=True):
        tiers = [TierTables.load(d, verify) for d in tier_dirs]
        runs = [RunTables.load(d, verify) for d in run_dirs]
        return cls(tiers, runs)

    def reference(self, tier: TierTables, ref_id: str) -> Reference:
        key = (id(tier), ref_id)
        if key not in self.refs:
            self.refs[key] = Reference.from_taxa(tier.taxa, ref_id)
        return self.refs[key]

    def run_rows(self):
        for rt in self.runs:
            for _, run in rt.runs.iterrows():
                preds = rt.predictions[rt.predictions["run_id"] == run["run_id"]]
                if preds.empty:
                    raise ContractError(f"run {run['run_id']}: no predictions")
                yield run, preds

    def cells(self, run, preds, split=None):
        """Yield (tier, cell dict, Aligned) for every cell the run covers."""
        pq = set(preds["query_id"])
        found = False
        for tier in self.tiers:
            tr = tier.truth[tier.truth["ref_id"] == run["ref_id"]]
            if tr.empty:
                continue
            ref = self.reference(tier, run["ref_id"])
            q = tier.queries[tier.queries["query_id"].isin(set(tr["query_id"]))]
            if split:
                q = q[q["split"].isin([split] if isinstance(split, str) else split)]
            for (t, ds, rg), qc in q.groupby(["tier", "dataset", "region"], sort=True):
                ids = sorted(qc["query_id"])
                have = [x for x in ids if x in pq]
                if not have:
                    continue
                if len(have) != len(ids):
                    raise ContractError(
                        f"run {run['run_id']} covers {len(have)} of {len(ids)} queries of cell "
                        f"{t}/{ds}/{rg}/{run['ref_id']} (partial output)")
                found = True
                al = M.align(tr[tr["query_id"].isin(ids)], tier.queries, preds, ref, ids)
                yield tier, {"tier": t, "dataset": ds, "region": rg, "ref_id": run["ref_id"]}, al, ref
        if not found:
            raise ContractError(f"run {run['run_id']} (ref_id {run['ref_id']}): no cell of any "
                                "tier is covered by its predictions")


def _subset(al: M.Aligned, mask) -> M.Aligned:
    idx = np.flatnonzero(mask)
    return M.Aligned(query_ids=[al.query_ids[i] for i in idx], cluster=al.cluster[idx],
                     split=al.split[idx],
                     t_status={r: v[idx] for r, v in al.t_status.items()},
                     t_set={r: [v[i] for i in idx] for r, v in al.t_set.items()},
                     c_call={r: v[idx] for r, v in al.c_call.items()},
                     c_set={r: [v[i] for i in idx] for r, v in al.c_set.items()},
                     conf={r: v[idx] for r, v in al.conf.items()},
                     c_members={r: [v[i] for i in idx] for r, v in al.c_members.items()},
                     homonym={r: v[idx] for r, v in al.homonym.items()},
                     unexpanded={r: v[idx] for r, v in al.unexpanded.items()})


def score(ctx: Context, rules=("G1", "G2", "G3"), ranks=RANKS, split=SPLIT_DEFAULT,
          n_boot=2000, seed=20261031) -> pd.DataFrame:
    rows = []
    for run, preds in ctx.run_rows():
        kind = confidence_kind(run)
        base = {"run_id": run["run_id"], "classifier": run["classifier"],
                "version": run["version"]}
        for tier, cell, al_all, ref in ctx.cells(run, preds):
            eval_mask = al_all.split == split
            if not eval_mask.any():
                continue
            al = _subset(al_all, eval_mask)
            calib = _subset(al_all, al_all.split == "calib") if (al_all.split == "calib").any() else None
            ab = tier.abundance
            for rank in ranks:
                for rule in rules:
                    df = M.outcomes(al, rank, rule)
                    c = M.totals(df)
                    rt = M.rates(c)
                    ci = M.bootstrap_ci(df, n_boot=n_boot, seed=seed) if n_boot else {}
                    keyrow = dict(base, **cell, split=split, rank=rank, rule=rule,
                                  n_K=c["K"], n_L=c["L"],
                                  n_excluded=int((df["outcome"] == "excluded").sum()),
                                  n_clusters=int(df.loc[df["outcome"] != "excluded", "cluster_id"].nunique()))
                    for k in M.COUNT_COLS:
                        rows.append(dict(keyrow, metric="count_" + k, value=c[k]))
                    for k, v in rt.items():
                        lo, hi = ci.get(k, (math.nan, math.nan))
                        rows.append(dict(keyrow, metric=k, value=v, ci_low=lo, ci_high=hi))
                    if rule != "G1":
                        continue
                    # M8 threshold-free and M9 calibration (G1 only)
                    op = str(run["run_id"]).rsplit(":", 1)[-1]
                    if op in ("default", "frozen"):
                        # one run per operating point (run_id suffix :default or :frozen);
                        # threshold-free metrics come from the ungated ':open' run
                        rows.append(dict(keyrow, metric="average_precision", value=math.nan,
                                         note=f"operating-point run ({op}); M8 is read from the :open run"))
                    elif kind in ("probability", "score"):
                        try:
                            cv = M.curve(df)
                            rows.append(dict(keyrow, metric="average_precision", value=M.average_precision(cv)))
                            rows.append(dict(keyrow, metric="aurc", value=M.aurc(cv)))
                            if calib is not None:
                                t = M.fit_threshold(M.outcomes(calib, rank, rule))
                                rows.append(dict(keyrow, metric="f1_threshold_calib", value=t))
                                rows.append(dict(keyrow, metric="f1_at_calib_threshold",
                                                 value=M.f1_at(df, t)["f1"] if t == t else math.nan))
                        except ContractError as e:
                            rows.append(dict(keyrow, metric="average_precision", value=math.nan, note=str(e)))
                    else:
                        rows.append(dict(keyrow, metric="average_precision", value=math.nan,
                                         note=f"confidence kind {kind}: no curve"))
                    if op in ("default", "frozen"):
                        rows.append(dict(keyrow, metric="ece_mass", value=math.nan,
                                         note=f"operating-point run ({op}); M9 is read from the :open run"))
                    elif kind == "probability":
                        try:
                            cal = M.calibration(df)
                            for k in ("ece_mass", "ece_width", "brier", "n_calibration"):
                                rows.append(dict(keyrow, metric=k, value=cal[k]))
                        except ContractError as e:
                            rows.append(dict(keyrow, metric="ece_mass", value=math.nan, note=str(e)))
                    else:
                        rows.append(dict(keyrow, metric="ece_mass", value=math.nan,
                                         note=f"not a probability ({kind})"))
                    if ab is not None and len(ab) and ab["query_id"].isin(set(al.query_ids)).any():
                        pr = M.profiles(al, ref, ab, rank)
                        for k in [c for c in pr.columns if c not in ("sample_id", "reads")]:
                            v = pr[k].dropna()
                            rows.append(dict(keyrow, metric=f"{k}_median", value=float(v.median()) if len(v) else math.nan,
                                             n_samples=int(len(v))))
                            rows.append(dict(keyrow, metric=f"{k}_mean", value=float(v.mean()) if len(v) else math.nan,
                                             n_samples=int(len(v))))
    out = pd.DataFrame(rows)
    out["rules_version"] = "v2"
    out["taxobench_version"] = __version__
    return out


def write_scores(df: pd.DataFrame, out_dir) -> Path:
    from .contract import write_table
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    fp = out / "scores.tsv"
    df.to_csv(fp, sep="\t", index=False, na_rep="NA", float_format="%.10g")
    from .contract import write_sidecar
    write_sidecar(fp, table="scores")
    return fp
