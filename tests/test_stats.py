"""Statistics known answers."""
import numpy as np
import pytest

from taxobench import metrics as M
from taxobench import stats as S

COLS = M.COUNT_COLS


def cluster_matrix(tp, k):
    """Clusters with k known rows each, tp of them TP (rest UC); one row per cluster."""
    m = np.zeros((len(tp), len(COLS)))
    m[:, COLS.index("K")] = k
    m[:, COLS.index("TP")] = tp
    m[:, COLS.index("UC")] = np.asarray(k) - np.asarray(tp)
    return m


def test_identical_runs_ci_contains_zero_and_holm_is_one():
    rng = np.random.default_rng(1)
    k = rng.integers(1, 10, 80)
    tp = rng.binomial(k, 0.7)
    MA = cluster_matrix(tp, k)
    r = S.paired_cluster_bootstrap(MA, MA.copy(), "recall", n_boot=2000, seed=3)
    assert r["delta"] == 0 and r["ci_low"] <= 0 <= r["ci_high"]
    assert r["p_boot"] == 1.0
    assert S.holm([1.0, 1.0])[0] == 1.0


def test_known_shift_coverage():
    """A known shift in recall: the 95% CI covers the true delta in >= 93% of
    200 seeded replicates."""
    rng = np.random.default_rng(20261031)
    true_a, true_b, n, k = 0.70, 0.75, 120, 5
    hits = 0
    for rep in range(200):
        tpa = rng.binomial(k, true_a, n)
        tpb = rng.binomial(k, true_b, n)
        r = S.paired_cluster_bootstrap(cluster_matrix(tpa, np.full(n, k)),
                                       cluster_matrix(tpb, np.full(n, k)), "recall",
                                       n_boot=1000, seed=rep)
        hits += r["ci_low"] <= true_b - true_a <= r["ci_high"]
    assert hits / 200 >= 0.93


def test_cluster_resampling_wider_than_row_resampling():
    rng = np.random.default_rng(7)
    n_cl, size = 40, 10
    pa = rng.uniform(0.3, 0.95, n_cl)             # strong cluster effect
    rows_a = rng.random((n_cl, size)) < pa[:, None]
    effect = rng.normal(0.05, 0.25, n_cl)          # the difference itself varies by cluster
    rows_b = rng.random((n_cl, size)) < np.clip(pa + effect, 0, 1)[:, None]
    cl_a = cluster_matrix(rows_a.sum(1), np.full(n_cl, size))
    cl_b = cluster_matrix(rows_b.sum(1), np.full(n_cl, size))
    row_a = cluster_matrix(rows_a.ravel().astype(int), np.ones(n_cl * size))
    row_b = cluster_matrix(rows_b.ravel().astype(int), np.ones(n_cl * size))
    w_cl = np.subtract(*[S.paired_cluster_bootstrap(cl_a, cl_b, "recall", 2000, 1)[k] for k in ("ci_high", "ci_low")])
    w_row = np.subtract(*[S.paired_cluster_bootstrap(row_a, row_b, "recall", 2000, 1)[k] for k in ("ci_high", "ci_low")])
    # rows of a cluster share their difference, so resampling rows understates
    # the uncertainty of the difference
    assert w_cl > w_row


def test_holm_and_bh_known_values():
    p = [0.01, 0.04, 0.03, 0.005]
    assert np.allclose(S.holm(p), [0.03, 0.06, 0.06, 0.02])
    assert np.allclose(S.bh(p), [0.02, 0.04, 0.04, 0.02])


def test_hodges_lehmann_and_wilcoxon():
    d = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10]
    r = S.hodges_lehmann(d)
    assert r["hl"] == pytest.approx(0.055)
    # exact 95% CI for n = 10 uses the 9th smallest and 9th largest Walsh averages
    w = sorted((a + b) / 2 for i, a in enumerate(d) for b in d[i:])
    assert r["hl_low"] == pytest.approx(w[8]) and r["hl_high"] == pytest.approx(w[-9])
    assert r["p_wilcoxon"] < 0.01


def test_mcnemar_exact():
    a = np.array([1, 1, 0, 0, 1, 0, 0, 0], bool)
    b = np.array([1, 0, 1, 1, 1, 1, 1, 0], bool)
    r = S.mcnemar_exact(a, b)
    assert (r["a_only"], r["b_only"]) == (1, 4)
    assert r["p_mcnemar"] == pytest.approx(0.375)
