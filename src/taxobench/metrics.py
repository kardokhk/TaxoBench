"""Metric catalogue M1 to M11 (docs/taxobench.md, section Metrics).

Per run and rank r, over the scored queries:
- K = rows whose truth status at r is ``known`` or ``unresolvable`` (a named
  truth; unresolvable truth is a set of several taxa);
- L = rows ``novel`` at r (the reference has no name there);
- ``unscorable`` rows are excluded.

A call is *named* at r when its call there is named, group or oov. At a rank q
where both the truth and the call are named, the call is consistent when its
taxa intersect the truth set (group rule G1/G3), and an ``oov`` call is never
consistent.

Outcomes at r (M1):
- K row, named at r, consistent at every rank <= r where both are named: TP;
- K row, named at r, otherwise: MC (misclassification);
- K row, not named at r: UC (under-classification), even when a higher rank
  is wrong; ``uc_wrong_above`` counts those UC rows;
- L row, named at r: OC (over-classification); not named: CA (correct
  abstention); ``ca_wrong_above`` counts CA rows wrong at a higher rank.
Precision at r judges only calls named at r (Edgar 2018).
So TP + MC + UC = K and OC + CA = L.

Group rules (M7): G1 group-containing (default; the truth needs only to be a
member), G2 strict (a group call is no call at its rank and below), G3
size-weighted (a consistent group call scores |call ∩ truth| / |call| as TP and
the rest as MC). Every function takes ``rule`` in {"G1", "G2", "G3"}.

M2 precision = TP / (TP + MC + OC), NaN when nothing is named.
M3 recall = TP / K.  M4 F1 from unrounded counts.  M5 Edgar accuracy =
TP / (K + OC).  M6 MCR = MC / K, UCR = UC / K, OCR = OC / L, coverage =
(TP + MC + OC) / (K + L).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

from . import RANKS
from .contract import ContractError, split_ids
from .taxonomy import Reference

RULES = ("G1", "G2", "G3")
COUNT_COLS = ["K", "L", "TP", "MC", "UC", "OC", "CA", "uc_wrong_above", "ca_wrong_above",
              "n_oov", "n_group", "group_size_sum", "n_homonym_resolved", "n_group_unexpanded"]
NAN = float("nan")

# ---------------------------------------------------------------------------
# Row-level outcomes


@dataclass
class Aligned:
    """Truth and calls of one run, aligned by query, ready for scoring."""
    query_ids: list
    cluster: np.ndarray            # cluster id per query (object)
    split: np.ndarray
    t_status: dict                 # rank -> np.array of status strings
    t_set: dict                    # rank -> list of frozenset (atoms)
    c_call: dict                   # rank -> np.array of call strings
    c_set: dict                    # rank -> list of frozenset (atoms)
    conf: dict                     # rank -> np.array float
    c_members: dict = None         # rank -> list of tuples of member frozensets
    homonym: dict = None           # rank -> np.array bool (truth or call used homonym resolution)
    unexpanded: dict = None        # rank -> np.array bool (truth or call holds an unexpanded group)


def align(truth: pd.DataFrame, queries: pd.DataFrame, preds: pd.DataFrame,
          ref: Reference, query_ids=None) -> Aligned:
    """Join truth (one ref_id) and the predictions of one run.

    Fails closed when a scored query has no prediction row at some rank, or a
    prediction names a taxon that is not in the reference.
    """
    if query_ids is None:
        query_ids = sorted(set(truth["query_id"]))
    q = queries.set_index("query_id")
    missing = [x for x in query_ids if x not in q.index]
    if missing:
        raise ContractError(f"{len(missing)} truth queries absent from queries table, e.g. {missing[:3]}")
    tw = truth.set_index(["query_id", "rank"])
    pw = preds.set_index(["query_id", "rank"])
    al = Aligned(query_ids=list(query_ids),
                 cluster=q.loc[query_ids, "cluster_id"].to_numpy(object),
                 split=q.loc[query_ids, "split"].to_numpy(object),
                 t_status={}, t_set={}, c_call={}, c_set={}, conf={},
                 c_members={}, homonym={}, unexpanded={})
    for r in RANKS:
        idx = pd.MultiIndex.from_arrays([query_ids, [r] * len(query_ids)])
        try:
            ts = tw.loc[idx]
        except KeyError:
            have = set(tw.index)
            miss = [i for i in idx if i not in have][:3]
            raise ContractError(f"truth rows missing for rank {r}, e.g. {miss}")
        try:
            ps = pw.loc[idx]
        except KeyError:
            have = set(pw.index)
            miss = [i for i in idx if i not in have][:3]
            raise ContractError(f"prediction rows missing for rank {r}, e.g. {miss} "
                                "(predictions must cover every scored query at every rank)")
        al.t_status[r] = ts["status"].to_numpy(object)
        tids = [split_ids(s) for s in ts["truth_set"]]
        cids = [split_ids(s) for s in ps["taxon_ids"]]
        al.t_set[r] = [ref.expand(x) for x in tids]
        al.c_call[r] = ps["call"].to_numpy(object)
        al.c_members[r] = [ref.members(x) for x in cids]
        al.c_set[r] = [frozenset().union(*m) for m in al.c_members[r]]
        al.homonym[r] = np.array([any(t in ref.homonym_groups for t in a + b) for a, b in zip(tids, cids)], bool)
        al.unexpanded[r] = np.array([any(t in ref.unexpanded_groups for t in a + b) for a, b in zip(tids, cids)], bool)
        al.conf[r] = pd.to_numeric(ps["confidence"], errors="coerce").to_numpy(float)
    return al


def outcomes(al: Aligned, rank: str, rule: str = "G1") -> pd.DataFrame:
    """One row per query: the M1 outcome at `rank` and per-row count columns
    (additive, so cell totals and bootstrap resamples are sums)."""
    if rule not in RULES:
        raise ValueError(rule)
    ri = RANKS.index(rank)
    n = len(al.query_ids)
    cols = {c: np.zeros(n) for c in COUNT_COLS}
    outcome = np.empty(n, dtype=object)
    for i in range(n):
        wrong = False
        strict_cut = False
        for q in RANKS[: ri + 1]:
            call = al.c_call[q][i]
            cset = al.c_set[q][i]
            if rule == "G2" and (strict_cut or (call == "group")):
                strict_cut = True
                call = "nocall"
            if call in ("named", "group", "oov") and al.t_status[q][i] in ("known", "unresolvable"):
                if call == "oov" or not (cset & al.t_set[q][i]):
                    wrong = True
        st = al.t_status[rank][i]
        call = al.c_call[rank][i]
        if rule == "G2" and strict_cut:
            call = "nocall"
        named = call in ("named", "group", "oov")
        if st == "unscorable":
            outcome[i] = "excluded"
            continue
        if call == "oov" and named:
            cols["n_oov"][i] = 1
        if call == "group":
            cols["n_group"][i] = 1
            cols["group_size_sum"][i] = len(al.c_members[rank][i])
        if al.homonym is not None and al.homonym[rank][i]:
            cols["n_homonym_resolved"][i] = 1
        if al.unexpanded is not None and al.unexpanded[rank][i]:
            cols["n_group_unexpanded"][i] = 1
        if st in ("known", "unresolvable"):
            cols["K"][i] = 1
            if named and not wrong:
                if rule == "G3":
                    mem = al.c_members[rank][i]
                    w = sum(1 for m in mem if m & al.t_set[rank][i]) / len(mem)
                    cols["TP"][i], cols["MC"][i] = w, 1 - w
                else:
                    cols["TP"][i] = 1
                outcome[i] = "TP"
            elif named:
                cols["MC"][i] = 1
                outcome[i] = "MC"
            else:
                cols["UC"][i] = 1
                cols["uc_wrong_above"][i] = 1 if wrong else 0
                outcome[i] = "UC"
        elif st == "novel":
            cols["L"][i] = 1
            if named:
                cols["OC"][i] = 1
                outcome[i] = "OC"
            else:
                cols["CA"][i] = 1
                cols["ca_wrong_above"][i] = 1 if wrong else 0
                outcome[i] = "CA"
        else:
            raise ContractError(f"unknown truth status {st!r}")
    df = pd.DataFrame(cols)
    df.insert(0, "outcome", outcome)
    df.insert(0, "split", al.split)
    df.insert(0, "cluster_id", al.cluster)
    df.insert(0, "query_id", al.query_ids)
    df["confidence"] = al.conf[rank]
    return df


# ---------------------------------------------------------------------------
# Counts -> rates (M2 to M6)


def _div(a, b):
    return a / b if b else NAN


def rates(c: dict) -> dict:
    named = c["TP"] + c["MC"] + c["OC"]
    p = _div(c["TP"], named)
    r = _div(c["TP"], c["K"])
    if p != p or r != r:
        f = NAN
    else:
        f = 2 * p * r / (p + r) if (p + r) else 0.0
    return {
        "precision": p, "recall": r, "f1": f,
        "edgar_accuracy": _div(c["TP"], c["K"] + c["OC"]),
        "mcr": _div(c["MC"], c["K"]), "ucr": _div(c["UC"], c["K"]),
        "ocr": _div(c["OC"], c["L"]),
        "coverage": _div(named, c["K"] + c["L"]),
        "oov_share": _div(c["n_oov"], named),
        "group_share": _div(c["n_group"], named),
        "mean_group_size": _div(c["group_size_sum"], c["n_group"]),
    }


def rates_vec(C: np.ndarray, cols=COUNT_COLS) -> dict:
    """rates() over rows of a count matrix (bootstrap resamples)."""
    i = {c: k for k, c in enumerate(cols)}
    with np.errstate(divide="ignore", invalid="ignore"):
        TP, MC, OC, K, L, UC = (C[:, i[x]] for x in ("TP", "MC", "OC", "K", "L", "UC"))
        named = TP + MC + OC
        p = np.where(named > 0, TP / named, np.nan)
        r = np.where(K > 0, TP / K, np.nan)
        f = np.where((p + r) > 0, 2 * p * r / (p + r), np.where(np.isnan(p) | np.isnan(r), np.nan, 0.0))
        return {"precision": p, "recall": r, "f1": f,
                "edgar_accuracy": np.where((K + OC) > 0, TP / (K + OC), np.nan),
                "mcr": np.where(K > 0, MC / K, np.nan), "ucr": np.where(K > 0, UC / K, np.nan),
                "ocr": np.where(L > 0, OC / L, np.nan),
                "coverage": np.where((K + L) > 0, named / (K + L), np.nan)}


def totals(df: pd.DataFrame) -> dict:
    return {c: float(df[c].sum()) for c in COUNT_COLS}


def cluster_counts(df: pd.DataFrame):
    g = df.groupby("cluster_id", sort=True)[COUNT_COLS].sum()
    return g.index.to_numpy(object), g.to_numpy(float)


def bootstrap_ci(df: pd.DataFrame, n_boot=2000, seed=20261031, level=0.95) -> dict:
    """Cluster bootstrap percentile CI of every rate, for one run."""
    _, M = cluster_counts(df)
    if len(M) < 2:
        return {}
    rng = np.random.default_rng(seed)
    W = rng.multinomial(len(M), np.full(len(M), 1 / len(M)), size=n_boot).astype(float)
    R = rates_vec(W @ M)
    a = (1 - level) / 2
    out = {}
    for k, v in R.items():
        v = v[~np.isnan(v)]
        if len(v):
            out[k] = (float(np.quantile(v, a)), float(np.quantile(v, 1 - a)))
    return out


# ---------------------------------------------------------------------------
# M8: threshold-free


def curve(df: pd.DataFrame) -> pd.DataFrame:
    """Precision, recall and coverage as the confidence threshold falls,
    one point per distinct confidence (ties enter together). Rows named at
    the rank without a confidence cannot be ranked and raise."""
    named = df[df["outcome"].isin(["TP", "MC", "OC"])]
    if named["confidence"].isna().any():
        raise ContractError("named calls without confidence; curves need a confidence per call")
    K = df["K"].sum()
    N = df["K"].sum() + df["L"].sum()
    g = named.groupby("confidence", sort=True)[["TP", "MC", "OC"]].sum().iloc[::-1]
    tp = g["TP"].cumsum().to_numpy()
    nm = (g["TP"] + g["MC"] + g["OC"]).cumsum().to_numpy()
    return pd.DataFrame({"threshold": g.index.to_numpy(float), "named": nm, "tp": tp,
                         "precision": tp / nm, "recall": tp / K if K else np.nan,
                         "coverage": nm / N if N else np.nan})


def average_precision(cv: pd.DataFrame) -> float:
    """Step-wise AP: sum over thresholds of (R_k - R_{k-1}) * P_k; no interpolation."""
    if cv.empty:
        return NAN
    r = np.concatenate([[0.0], cv["recall"].to_numpy()])
    return float(np.sum(np.diff(r) * cv["precision"].to_numpy()))


def aurc(cv: pd.DataFrame) -> float:
    """Area under the risk-coverage curve (risk = 1 - precision), step-wise."""
    if cv.empty:
        return NAN
    c = np.concatenate([[0.0], cv["coverage"].to_numpy()])
    return float(np.sum(np.diff(c) * (1 - cv["precision"].to_numpy())))


def f1_at(df: pd.DataFrame, t: float) -> dict:
    d = df.copy()
    below = d["outcome"].isin(["TP", "MC", "OC"]) & ~(d["confidence"] >= t)
    # a call below the threshold becomes no call: TP/MC -> UC, OC -> CA
    d.loc[below & (d["K"] == 1), ["TP", "MC", "UC"]] = [0, 0, 1]
    d.loc[below & (d["L"] == 1), ["OC", "CA"]] = [0, 1]
    return rates(totals(d))


def fit_threshold(calib: pd.DataFrame) -> float:
    """F1-optimal threshold on the calibration part (smallest among ties)."""
    cv = curve(calib)
    if cv.empty:
        return NAN
    K = calib["K"].sum()
    p, r = cv["precision"].to_numpy(), cv["tp"].to_numpy() / K if K else np.zeros(len(cv))
    f = np.where(p + r > 0, 2 * p * r / (p + r), 0)
    best = np.flatnonzero(f == f.max())
    return float(cv["threshold"].to_numpy()[best].min())


# ---------------------------------------------------------------------------
# M9: calibration


def calibration(df: pd.DataFrame, n_mass=15, n_width=10) -> dict:
    """Over named calls with a confidence: correct = TP. ECE with equal-mass
    (primary) and equal-width bins, Brier score, and the reliability table."""
    d = df[df["outcome"].isin(["TP", "MC", "OC"]) & df["confidence"].notna()]
    if d.empty:
        return {"ece_mass": NAN, "ece_width": NAN, "brier": NAN, "n_calibration": 0}
    p = d["confidence"].to_numpy(float)
    y = (d["outcome"] == "TP").to_numpy(float)
    if ((p < -1e-9) | (p > 1 + 1e-9)).any():
        raise ContractError("calibration needs confidences in [0, 1]")
    p = np.clip(p, 0.0, 1.0)            # floating-point excess such as 1.0000000000000147
    # equal-mass bins from quantile edges; tied confidences stay in one bin
    inner = np.unique(np.quantile(p, np.linspace(0, 1, n_mass + 1))[1:-1])
    bm = np.searchsorted(inner, p, side="left")
    ece_m, table = 0.0, []
    for k in np.unique(bm):
        idx = bm == k
        gap = abs(p[idx].mean() - y[idx].mean())
        ece_m += idx.sum() / len(p) * gap
        table.append(("mass", int(k), int(idx.sum()), float(p[idx].mean()), float(y[idx].mean())))
    edges = np.linspace(0, 1, n_width + 1)
    b = np.clip(np.digitize(p, edges[1:-1], right=True), 0, n_width - 1)
    ece_w = 0.0
    for k in range(n_width):
        m = b == k
        if m.any():
            ece_w += m.sum() / len(p) * abs(p[m].mean() - y[m].mean())
            table.append(("width", k, int(m.sum()), float(p[m].mean()), float(y[m].mean())))
    return {"ece_mass": ece_m, "ece_width": ece_w, "brier": float(np.mean((p - y) ** 2)),
            "n_calibration": int(len(p)), "reliability": table}


# ---------------------------------------------------------------------------
# M10, M11: read-weighted profiles per sample


def _lca(ref: Reference, ids) -> str:
    paths = [ref.path(t) for t in ids]
    if not paths:
        return ""
    common = ""
    for level in zip(*paths):
        if all(x == level[0] for x in level):
            common = level[0]
        else:
            break
    return common


def node_of(ref: Reference, sets_by_rank: list, calls_by_rank: list | None = None) -> str:
    """Deepest node a call (or truth) assigns: the LCA of the taxa at the
    deepest rank that is resolved without an error above it. "" = root."""
    node = ""
    for k, s in enumerate(sets_by_rank):
        call = calls_by_rank[k] if calls_by_rank is not None else ("named" if s else "nocall")
        if call in ("nocall", "oov") or not s:
            break
        node = _lca(ref, s) if len(s) > 1 else next(iter(s))
    return node


def profile_parts(al: Aligned, ref: Reference, rank: str, side: str, i: int) -> str:
    """Profile key of query i at `rank`: a taxon_id, or one of the explicit
    parts "~unassigned" (no call, group above the rank, novel truth),
    "~oov:<i>" (a name absent from the reference; never merges with others)."""
    if side == "truth":
        st, s = al.t_status[rank][i], al.t_set[rank][i]
        if st == "known" and len(s) == 1:
            return next(iter(s))
        return "~unassigned"
    call, s = al.c_call[rank][i], al.c_set[rank][i]
    if call == "oov":
        return f"~oov:{i}"
    if call == "named" and len(s) == 1:
        return next(iter(s))
    return "~unassigned"


def _normalise(d: dict, keep_unassigned: bool) -> dict:
    if not keep_unassigned:
        d = {k: v for k, v in d.items() if k != "~unassigned"}
    tot = sum(d.values())
    return {k: v / tot for k, v in d.items()} if tot > 0 else {}


def l1_bc(p: dict, q: dict) -> tuple:
    keys = set(p) | set(q)
    l1 = sum(abs(p.get(k, 0) - q.get(k, 0)) for k in keys)
    bc = 1 - sum(min(p.get(k, 0), q.get(k, 0)) for k in keys)
    return l1, bc


def weighted_unifrac(ref: Reference, p: dict, q: dict) -> float:
    """Earth mover's distance on the reference taxonomy, unit branch lengths
    (OPAL style). Keys are node taxon_ids ("" = root); masses sum to 1 each."""
    diff = {}
    for k, v in p.items():
        diff[k] = diff.get(k, 0) + v
    for k, v in q.items():
        diff[k] = diff.get(k, 0) - v
    sub = {}
    for k, v in diff.items():
        node = k
        while node:
            sub[node] = sub.get(node, 0) + v
            node = ref.parent[node]
    return float(sum(abs(v) for v in sub.values()))


def aitchison(p_counts: dict, q_counts: dict, components: list, pseudo=0.5) -> float:
    a = np.array([p_counts.get(c, 0) + pseudo for c in components], float)
    b = np.array([q_counts.get(c, 0) + pseudo for c in components], float)
    ca, cb = np.log(a) - np.log(a).mean(), np.log(b) - np.log(b).mean()
    return float(np.sqrt(np.sum((ca - cb) ** 2)))


def profiles(al: Aligned, ref: Reference, abundance: pd.DataFrame, rank: str) -> pd.DataFrame:
    """Per sample: L1 and Bray-Curtis under convention A (renormalised over
    named mass) and B (unassigned kept as a part), weighted UniFrac (A and B)
    and Aitchison (M11; components = truth taxa, one pooled "other named"
    part and "unassigned"; pseudocount 0.5 read). Group mass sits at the
    members' LCA node and is never resolved with the truth."""
    pos = {q: i for i, q in enumerate(al.query_ids)}
    ri = RANKS.index(rank)
    rows = []
    ab = abundance[abundance["query_id"].isin(pos)]
    for sample, g in ab.groupby("sample_id", sort=True):
        tc, cc, tn, cn = {}, {}, {}, {}
        for qid, cnt in zip(g["query_id"], pd.to_numeric(g["count"])):
            i = pos[qid]
            if al.t_status[rank][i] == "unscorable" or cnt <= 0:
                continue
            kt, kc = profile_parts(al, ref, rank, "truth", i), profile_parts(al, ref, rank, "call", i)
            tc[kt] = tc.get(kt, 0) + cnt
            cc[kc] = cc.get(kc, 0) + cnt
            nt = node_of(ref, [al.t_set[r][i] if al.t_status[r][i] in ("known", "unresolvable") else frozenset()
                               for r in RANKS[: ri + 1]])
            nc = node_of(ref, [al.c_set[r][i] for r in RANKS[: ri + 1]],
                         [al.c_call[r][i] for r in RANKS[: ri + 1]])
            tn[nt] = tn.get(nt, 0) + cnt
            cn[nc] = cn.get(nc, 0) + cnt
        if not tc:
            continue
        row = {"sample_id": sample, "reads": float(sum(tc.values()))}
        for conv, keep in (("A", False), ("B", True)):
            p, q = _normalise(tc, keep), _normalise(cc, keep)
            row[f"l1_{conv}"], row[f"bray_curtis_{conv}"] = l1_bc(p, q) if p and q else (NAN, NAN)
            tn2 = tn if keep else {k: v for k, v in tn.items() if k}
            cn2 = cn if keep else {k: v for k, v in cn.items() if k}
            pt, pc = _normalise(tn2, True), _normalise(cn2, True)
            row[f"unifrac_{conv}"] = weighted_unifrac(ref, pt, pc) if pt and pc else NAN
        truth_taxa = sorted(k for k in tc if not k.startswith("~"))
        comp_c = {k: v for k, v in cc.items() if k in truth_taxa}
        comp_c["~other"] = sum(v for k, v in cc.items() if k not in truth_taxa and k != "~unassigned")
        comp_c["~unassigned"] = cc.get("~unassigned", 0)
        comp_t = {k: v for k, v in tc.items() if k in truth_taxa}
        comp_t["~other"] = 0
        comp_t["~unassigned"] = tc.get("~unassigned", 0)
        row["aitchison"] = aitchison(comp_t, comp_c, truth_taxa + ["~other", "~unassigned"])
        rows.append(row)
    return pd.DataFrame(rows)
