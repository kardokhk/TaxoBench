"""Statistics for release comparison.

- Paired cluster bootstrap of a difference: the same resampled clusters (with
  multiplicity) for A and B; percentile CI; two-sided bootstrap P.
- Exact McNemar on discordant TP rows.
- Wilcoxon signed-rank over cells with a Hodges-Lehmann estimate and CI.
- Holm (confirmatory family) and Benjamini-Hochberg (exploratory) adjustment.
"""
from __future__ import annotations

import math

import numpy as np
from scipy import stats as sps

from .metrics import COUNT_COLS, rates_vec


def paired_cluster_bootstrap(MA: np.ndarray, MB: np.ndarray, metric: str, n_boot=10000,
                             seed=20261031, level=0.95) -> dict:
    """MA, MB: cluster x COUNT_COLS count matrices over the same clusters
    (same order). Returns A, B, delta = B - A, CI of delta and bootstrap P.

    P is two-sided: 2 * min(share of resampled deltas <= 0, share >= 0), with
    the +1 correction, capped at 1.
    """
    if MA.shape != MB.shape:
        raise ValueError("A and B must be aligned on the same clusters")
    a = rates_vec(MA.sum(0, keepdims=True))[metric][0]
    b = rates_vec(MB.sum(0, keepdims=True))[metric][0]
    n = MA.shape[0]
    out = {"A": float(a), "B": float(b), "delta": float(b - a), "n_clusters": int(n),
           "n_boot": int(n_boot), "seed": int(seed)}
    if n < 2:
        out.update(ci_low=math.nan, ci_high=math.nan, p_boot=math.nan)
        return out
    rng = np.random.default_rng(seed)
    d = np.empty(0)
    chunk = max(1, min(n_boot, int(2e7 // max(n, 1))))
    done = 0
    parts = []
    while done < n_boot:
        k = min(chunk, n_boot - done)
        W = rng.multinomial(n, np.full(n, 1 / n), size=k).astype(float)
        da = rates_vec(W @ MA)[metric]
        db = rates_vec(W @ MB)[metric]
        parts.append(db - da)
        done += k
    d = np.concatenate(parts)
    d = d[~np.isnan(d)]
    if len(d) == 0:
        out.update(ci_low=math.nan, ci_high=math.nan, p_boot=math.nan)
        return out
    alpha = (1 - level) / 2
    out["ci_low"] = float(np.quantile(d, alpha))
    out["ci_high"] = float(np.quantile(d, 1 - alpha))
    le = (np.sum(d <= 0) + 1) / (len(d) + 1)
    ge = (np.sum(d >= 0) + 1) / (len(d) + 1)
    out["p_boot"] = float(min(1.0, 2 * min(le, ge)))
    return out


def mcnemar_exact(a_correct: np.ndarray, b_correct: np.ndarray) -> dict:
    """Exact two-sided McNemar test on paired binary outcomes."""
    a = np.asarray(a_correct, bool)
    b = np.asarray(b_correct, bool)
    n01 = int(np.sum(~a & b))       # B right, A wrong
    n10 = int(np.sum(a & ~b))
    n = n01 + n10
    p = 1.0 if n == 0 else float(min(1.0, 2 * sps.binom.cdf(min(n01, n10), n, 0.5)))
    return {"b_only": n01, "a_only": n10, "p_mcnemar": p}


def hodges_lehmann(d, level=0.95) -> dict:
    """Wilcoxon signed-rank P and the Hodges-Lehmann estimate (median of the
    Walsh averages) with its CI from the signed-rank distribution (normal
    approximation of the critical rank for n > 25, exact table otherwise)."""
    d = np.asarray([x for x in d if not math.isnan(x)], float)
    n = len(d)
    out = {"n_cells": n, "median_delta": float(np.median(d)) if n else math.nan}
    if n == 0:
        out.update(hl=math.nan, hl_low=math.nan, hl_high=math.nan, p_wilcoxon=math.nan)
        return out
    w = (d[:, None] + d[None, :]) / 2
    walsh = np.sort(w[np.triu_indices(n)])
    out["hl"] = float(np.median(walsh))
    if np.all(d == 0) or n < 2:
        out.update(hl_low=out["hl"], hl_high=out["hl"], p_wilcoxon=1.0)
        return out
    out["p_wilcoxon"] = float(sps.wilcoxon(d, zero_method="wilcox", alternative="two-sided").pvalue) \
        if np.any(d != 0) else 1.0
    alpha = 1 - level
    N = len(walsh)
    if n <= 25:
        # exact null distribution of the signed-rank statistic
        counts = np.zeros(n * (n + 1) // 2 + 1)
        counts[0] = 1
        for k in range(1, n + 1):
            counts[k:] = counts[k:] + counts[:-k].copy()
        cdf = np.cumsum(counts) / counts.sum()
        # c = T_low + 1, where T_low is the largest t with P(T <= t) <= alpha/2;
        # the CI is [W_(c), W_(N + 1 - c)] (1-indexed order statistics)
        c = int(np.searchsorted(cdf, alpha / 2, side="right"))
    else:
        z = sps.norm.ppf(1 - alpha / 2)
        c = int(round(n * (n + 1) / 4 - z * math.sqrt(n * (n + 1) * (2 * n + 1) / 24)))
    if c < 1:
        out.update(hl_low=float(walsh[0]), hl_high=float(walsh[-1]))
    else:
        out.update(hl_low=float(walsh[c - 1]), hl_high=float(walsh[N - c]))
    return out


def holm(p) -> np.ndarray:
    p = np.asarray(p, float)
    m = len(p)
    order = np.argsort(p)
    adj = np.empty(m)
    run = 0.0
    for k, i in enumerate(order):
        run = max(run, min(1.0, (m - k) * p[i]))
        adj[i] = run
    return adj


def bh(p) -> np.ndarray:
    p = np.asarray(p, float)
    m = len(p)
    if m == 0:
        return p
    order = np.argsort(p)
    adj = np.empty(m)
    prev = 1.0
    for k in range(m - 1, -1, -1):
        i = order[k]
        prev = min(prev, p[i] * m / (k + 1))
        adj[i] = min(1.0, prev)
    return adj


__all__ = ["paired_cluster_bootstrap", "mcnemar_exact", "hodges_lehmann", "holm", "bh",
           "COUNT_COLS"]
