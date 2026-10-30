"""Regression: the v1 rules reproduce the published S3 Table (round 3) values for
gutmock V5-V7 (corrected truth of 2026-10-06, all eight classifiers, five ranks,
precision / recall / F / specificity / L1 / L2). Fixture: tests/fixtures/gutmock_V5-V7/."""
import csv
import json
import unittest

import _paths
import evaluate_results as ER

ev = _paths.load_evaluator()
FX = _paths.FIXTURES / "gutmock_V5-V7"


def fmt(metric, v):
    """Published tables print P/R/F to 3 decimals (0 as '0') and specificity to 5."""
    if metric == "specificity":
        return round(float(v), 5)
    return round(float(v), 3)


class TestV1ReproducesPublished(unittest.TestCase):
    def setUp(self):
        ev._UNIVERSE_CACHE.clear()
        with open(FX / "universe_V5-V7.json") as fh:
            universe = json.load(fh)
        for lv, n in universe.items():
            ev._UNIVERSE_CACHE[("V5-V7", int(lv))] = n
        true_fp = FX / "data/mock_community/gutmock/V5-V7/trueish_taxonomy_gutmock.tsv"
        self.rows = ev._evaluate_combo("gutmock", "V5-V7", FX / "results/gutmock/V5-V7", true_fp,
                                       ER.precision_recall_fscore, rules="v1")
        with open(FX / "published_table_S3_rows.csv") as fh:
            self.pub = list(csv.DictReader(fh))

    def test_all_values(self):
        got = {(r["method"], r["level"]): r for r in self.rows}
        n = 0
        for p in self.pub:
            r = got[(p["classifier"], p["rank"])]
            for m_pub, m_ev in (("precision", "precision"), ("recall", "recall"), ("f_score", "f-score"),
                                ("specificity", "specificity")):
                self.assertEqual(fmt(m_pub, r[m_ev]), fmt(m_pub, p[m_pub]), (p["classifier"], p["rank"], m_pub))
                n += 1
            for m in ("l1_error", "l2_error"):
                self.assertAlmostEqual(float(r[m]), float(p[m]), places=3, msg=(p["classifier"], p["rank"], m))
                n += 1
        self.assertEqual(n, 40 * 6)

    def test_v2_changes_only_what_it_should(self):
        # P/R/F are unchanged by v2 on this cell for classifiers that emit no group labels
        # and write full lineages (all except SPINGO's taxon set and Taxonova's groups).
        v2 = ev._evaluate_combo("gutmock", "V5-V7", FX / "results/gutmock/V5-V7",
                                FX / "data/mock_community/gutmock/V5-V7/trueish_taxonomy_gutmock.tsv",
                                ER.precision_recall_fscore, rules="v2")
        a = {(r["method"], r["level"]): r for r in self.rows}
        b = {(r["method"], r["level"]): r for r in v2}
        for k in a:
            for m in ("precision", "recall", "f-score"):
                self.assertEqual(a[k][m], b[k][m], (k, m))


if __name__ == "__main__":
    unittest.main()
