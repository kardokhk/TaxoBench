"""Golden test of the frozen v1 scoring rules over every published cell.

Replaces the one-cell check of test_regression_v1.py as the guard of the paper
numbers. It runs scripts/04_evaluate.py --rules v1 as a user would and compares
the resulting all_metrics.csv, after a stable sort on the key columns, string
for string with a target table. Every value must be identical (difference 0).

Two targets, each skipped when its data are absent:

1. Zenodo deposit (public; doi:10.5281/zenodo.23808921, version 1.0).
   Inputs: paper_classifier_outputs, paper_truth_v1, reference_taxonomies.
   Target: paper_metric_tables/rescore/all_metrics_v1.csv.
   47 published cells x 8 classifiers, 15,040 values.
   Set TAXOBENCH_PAPER_DEPOSIT to the directory the bundles were unpacked into
   (scripts/fetch_paper_deposit.sh does this, with sha256 checks).

2. Archived benchmark run (MeluXina only; not public).
   61 mock cells x 11 methods (the eight classifiers plus three Taxonova
   variants), 26,840 values of the benchmark's own results/evaluation/all_metrics.csv,
   with the truth as archived (before the 2026-10-06 gut mock V5-V7 correction).
   Set TAXOBENCH_ARCHIVE_ROOT to a taxobench_final tree holding results/,
   data/mock_community/ and data/ref_dbs/, and TAXOBENCH_ARCHIVE_TARGET to the
   archived all_metrics.csv.

Determinism: v1 sums the L1 and L2 terms over a Python set of taxon strings,
so their order, and in rare cases the sixth decimal of a value, depends on the
string hash seed (observed once in nine runs: one L2 value 0.070313 instead of
0.070312; the published tables print four decimals and are unaffected). The
test therefore runs the evaluator with PYTHONHASHSEED=0. The v1 code is not
changed; v2 sums in sorted order.

Pinned interpreter and packages for this job: requirements-golden.txt.
"""
import csv
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import _paths

ROOT = _paths.ROOT
KEY = ["mock", "region", "method", "param", "level", "metric"]


def _sorted_rows(fp):
    with open(fp, newline="") as fh:
        rows = list(csv.reader(fh))
    header, body = rows[0], rows[1:]
    idx = [header.index(k) for k in KEY]
    body.sort(key=lambda r: [r[i] for i in idx])      # list.sort is stable
    return header, body


def _rescore(results, truth, ref, out):
    env = dict(os.environ, PYTHONHASHSEED="0")
    r = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "04_evaluate.py"), "--rules", "v1",
         "--results-dir", str(results), "--truth-dir", str(truth),
         "--ref-dir", str(ref), "--out-dir", str(out)],
        capture_output=True, text=True, env=env)
    if r.returncode != 0:
        raise AssertionError(f"04_evaluate.py exit {r.returncode}\n{r.stderr[-3000:]}")
    return out / "all_metrics.csv"


class _Golden:
    n_cells = n_methods = n_values = None

    def _compare(self, got_fp, want_fp):
        hg, got = _sorted_rows(got_fp)
        hw, want = _sorted_rows(want_fp)
        self.assertEqual(hg, hw)
        self.assertEqual(len(got), self.n_values)
        self.assertEqual(len(want), self.n_values)
        i_m, i_r, i_c = hw.index("mock"), hw.index("region"), hw.index("method")
        self.assertEqual(len({(r[i_m], r[i_r]) for r in want}), self.n_cells)
        self.assertEqual(len({r[i_c] for r in want}), self.n_methods)
        diffs = [(g, w) for g, w in zip(got, want) if g != w]
        self.assertEqual(diffs[:10], [], f"{len(diffs)} of {len(want)} rows differ")
        self.assertEqual(got, want)


@unittest.skipUnless(os.environ.get("TAXOBENCH_PAPER_DEPOSIT"),
                     "set TAXOBENCH_PAPER_DEPOSIT to the unpacked Zenodo deposit")
class TestGoldenDeposit(_Golden, unittest.TestCase):
    n_cells, n_methods, n_values = 47, 8, 15040

    def test_all_published_cells(self):
        d = Path(os.environ["TAXOBENCH_PAPER_DEPOSIT"])
        with tempfile.TemporaryDirectory() as tmp:
            got = _rescore(d / "paper_classifier_outputs" / "results", d / "paper_truth_v1",
                           d / "reference_taxonomies" / "ref_dbs", Path(tmp) / "eval")
            self._compare(got, d / "paper_metric_tables" / "rescore" / "all_metrics_v1.csv")


@unittest.skipUnless(os.environ.get("TAXOBENCH_ARCHIVE_ROOT")
                     and os.environ.get("TAXOBENCH_ARCHIVE_TARGET"),
                     "set TAXOBENCH_ARCHIVE_ROOT and TAXOBENCH_ARCHIVE_TARGET (MeluXina only)")
class TestGoldenArchive61(_Golden, unittest.TestCase):
    n_cells, n_methods, n_values = 61, 11, 26840

    def test_all_archived_cells(self):
        a = Path(os.environ["TAXOBENCH_ARCHIVE_ROOT"])
        with tempfile.TemporaryDirectory() as tmp:
            got = _rescore(a / "results", a / "data", a / "data" / "ref_dbs",
                           Path(tmp) / "eval")
            self._compare(got, Path(os.environ["TAXOBENCH_ARCHIVE_TARGET"]))


if __name__ == "__main__":
    unittest.main()
