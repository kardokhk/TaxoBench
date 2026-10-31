"""Fail-closed behaviour of the run and scoring scripts (pipeline audit B1 to B4,
B10, 2026-10-31). No classifier, QIIME 2 or data tree is needed: fake classifier
commands are shell one-liners.
"""
import os
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path

import _paths  # noqa: F401
import _common
import generate_commands as G

ROOT = _paths.ROOT
FX = _paths.FIXTURES / "gutmock_V5-V7"
TRUTH = FX / "data/mock_community/gutmock/V5-V7/trueish_taxonomy_gutmock.tsv"
PRED = FX / "results/gutmock/V5-V7/taxonova/confidence_threshold-0.7/taxonomy_final.tsv"


def _write_fasta(fp, ids):
    with open(fp, "w") as fh:
        for i in ids:
            fh.write(f">{i}\nACGT\n")


def _write_tax(fp, ids, header=True):
    with open(fp, "w") as fh:
        if header:
            fh.write("Feature ID\tTaxon\tConfidence\n")
        for i in ids:
            fh.write(f"{i}\td__Bacteria; p__X\t1.0; 1.0\n")


class TestRunGuard(unittest.TestCase):
    """_common.run_guarded: a failed rerun must not leave the old output."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.out = self.tmp / "out"
        self.out.mkdir()
        self.fa = self.tmp / "q.fasta"
        _write_fasta(self.fa, ["a", "b", "c"])

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _cmd_writing(self, ids, exit_code=0):
        src = self.tmp / f"src_{len(ids)}_{exit_code}.tsv"
        _write_tax(src, ids)
        return f"cp {src} {self.out}/taxonomy_final.tsv && exit {exit_code}"

    def test_success(self):
        _, ok, probs = _common.run_guarded(self._cmd_writing(["a", "b", "c"]),
                                           self.out, self.out / "run.log", self.fa)
        self.assertTrue(ok, probs)
        st = _common.read_run_status(self.out)
        self.assertTrue(st["ok"])
        self.assertEqual(st["returncode"], 0)

    def test_failed_rerun_removes_stale_output(self):
        _common.run_guarded(self._cmd_writing(["a", "b", "c"]), self.out,
                            self.out / "run.log", self.fa)
        self.assertTrue((self.out / "taxonomy_final.tsv").exists())
        # B1: the rerun fails before writing anything (command not found)
        _, ok, probs = _common.run_guarded("no-such-classifier-binary --run", self.out,
                                           self.out / "run.log", self.fa)
        self.assertFalse(ok)
        self.assertFalse((self.out / "taxonomy_final.tsv").exists())
        self.assertIn("exit status 127", probs[0])
        self.assertFalse(_common.read_run_status(self.out)["ok"])

    def test_nonzero_exit_with_output_is_failure(self):
        _, ok, probs = _common.run_guarded(self._cmd_writing(["a", "b", "c"], exit_code=3),
                                           self.out, self.out / "run.log", self.fa)
        self.assertFalse(ok)
        self.assertFalse((self.out / "taxonomy_final.tsv").exists())
        self.assertTrue((self.out / "taxonomy_final.tsv.failed").exists())

    def test_partial_output_is_failure(self):
        # B3: half the queries missing
        _, ok, probs = _common.run_guarded(self._cmd_writing(["a"]), self.out,
                                           self.out / "run.log", self.fa)
        self.assertFalse(ok)
        self.assertTrue(any("2 of 3 query IDs missing" in p for p in probs), probs)

    def test_duplicates_and_unexpected_are_failure(self):
        _, ok, probs = _common.run_guarded(self._cmd_writing(["a", "b", "c", "c", "z"]),
                                           self.out, self.out / "run.log", self.fa)
        self.assertFalse(ok)
        self.assertTrue(any("duplicated" in p for p in probs), probs)
        self.assertTrue(any("unexpected" in p for p in probs), probs)

    def test_header_only_output_is_failure(self):
        # B4: a parser that writes only a header
        _, ok, probs = _common.run_guarded(self._cmd_writing([]), self.out,
                                           self.out / "run.log", self.fa)
        self.assertFalse(ok)


class TestEvaluateFailsClosed(unittest.TestCase):
    """04_evaluate.py refuses partial, duplicated or failed outputs and exits 1."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.res = self.tmp / "results"
        self.truth = self.tmp / "truth"
        tdir = self.truth / "mock_community" / "m" / "V4"
        tdir.mkdir(parents=True)
        shutil.copy(TRUTH, tdir / "trueish_taxonomy_m.tsv")
        self.pdir = self.res / "m" / "V4" / "taxonova" / "p"
        self.pdir.mkdir(parents=True)
        shutil.copy(PRED, self.pdir / "taxonomy_final.tsv")

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _run(self):
        return subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "04_evaluate.py"), "--rules", "v1",
             "--results-dir", str(self.res), "--truth-dir", str(self.truth),
             "--ref-dir", str(self.tmp / "noref"), "--out-dir", str(self.tmp / "eval")],
            capture_output=True, text=True)

    def test_complete_output_scores(self):
        r = self._run()
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertTrue((self.tmp / "eval" / "all_metrics.csv").exists())

    def test_partial_output_exit_1(self):
        lines = PRED.read_text().splitlines()
        with open(self.pdir / "taxonomy_final.tsv", "w") as fh:
            fh.write("\n".join(lines[: len(lines) // 2]) + "\n")
        r = self._run()
        self.assertEqual(r.returncode, 1)
        self.assertIn("query IDs missing", r.stderr)
        self.assertFalse((self.tmp / "eval").exists())

    def test_duplicated_rows_exit_1(self):
        lines = PRED.read_text().splitlines()
        with open(self.pdir / "taxonomy_final.tsv", "w") as fh:
            fh.write("\n".join(lines + lines[1:]) + "\n")
        r = self._run()
        self.assertEqual(r.returncode, 1)
        self.assertIn("duplicated", r.stderr)

    def test_failed_run_status_exit_1(self):
        _common.write_run_status(self.pdir, ok=False, returncode=127, started=time.time(),
                                 problems=["exit status 127"])
        r = self._run()
        self.assertEqual(r.returncode, 1)
        self.assertIn("failed run", r.stderr)

    def test_stale_output_exit_1(self):
        old = time.time() - 3600
        os.utime(self.pdir / "taxonomy_final.tsv", (old, old))
        _common.write_run_status(self.pdir, ok=True, returncode=0, started=time.time(),
                                 problems=[])
        r = self._run()
        self.assertEqual(r.returncode, 1)
        self.assertIn("predates", r.stderr)

    def test_nothing_to_score_exit_1(self):
        shutil.rmtree(self.res / "m")
        (self.res / "x").mkdir()
        r = self._run()
        self.assertEqual(r.returncode, 1)

    def test_predictions_without_truth_exit_1(self):
        d = self.res / "m" / "V9" / "taxonova" / "p"
        d.mkdir(parents=True)
        shutil.copy(PRED, d / "taxonomy_final.tsv")
        r = self._run()
        self.assertEqual(r.returncode, 1)
        self.assertIn("no trueish taxonomy", r.stderr)


class TestCondaBin(unittest.TestCase):
    """B10: the configured conda executable is used by every builder."""

    def setUp(self):
        self._old = os.environ.get("TAXOBENCH_CONDA_BIN")
        os.environ["TAXOBENCH_CONDA_BIN"] = "/opt/mm/bin/micromamba"

    def tearDown(self):
        if self._old is None:
            os.environ.pop("TAXOBENCH_CONDA_BIN", None)
        else:
            os.environ["TAXOBENCH_CONDA_BIN"] = self._old

    def test_builders_use_override(self):
        ref = {"seqs": "/ref/V4/sequences.fasta", "model": "/ref/V4/m.pt",
               "database": "/db", "database_fasta": __file__, "database_tax": __file__}
        cmds = [
            G.build_taxonova_command("/in/q", "/out", {"confidence_threshold": 0.7}, ref,
                                     "/src", "/env", True),
            G.build_idtaxa_command("/in/q", "/out", {}, ref, "/src", "/env", True),
            G.build_kraken2_command("/in/q", "/out", {}, ref, "/src", "/env", 4, True),
            G.build_mapseq_command("/in/q", "/out", {}, ref, "/src", "/env", 4, True),
            G.build_spingo_command("/in/q", "/out", {}, ref, "/src",
                                   {"spingo_env": "/e1", "processing_env": "/e2"}, 4, True),
        ]
        for cmd in cmds:
            self.assertIn("/opt/mm/bin/micromamba run -p", cmd)
            self.assertNotIn("conda run", cmd)

    def test_no_literal_conda_run_left(self):
        src = (ROOT / "src" / "generate_commands.py").read_text()
        self.assertNotIn("conda run -p", src.replace("`<conda> run -p", ""))


class TestScriptsExitNonZero(unittest.TestCase):
    """B2: parsers exit non-zero instead of writing an empty table."""

    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def _py(self, script, *args):
        return subprocess.run([sys.executable, str(ROOT / "src" / script), *map(str, args)],
                              capture_output=True, text=True)

    def test_idtaxa_unreadable(self):
        raw = self.tmp / "raw.tsv"
        raw.write_text("")
        r = self._py("process_idtaxa_output.py", raw, self.tmp / "out.tsv")
        self.assertNotEqual(r.returncode, 0)
        self.assertFalse((self.tmp / "out.tsv").exists())

    def test_mapseq_nothing_parsed(self):
        raw = self.tmp / "raw.mseq"
        raw.write_text("# header only\n")
        r = self._py("process_mapseq_output.py", raw, self.tmp / "out.tsv")
        self.assertNotEqual(r.returncode, 0)
        self.assertFalse((self.tmp / "out.tsv").exists())

    def test_taxonova_empty(self):
        raw = self.tmp / "raw.tsv"
        raw.write_text("sequence_id\tdomain\tdomain_confidence\n")
        r = self._py("process_taxonova_output.py", raw, self.tmp / "out.tsv", "False")
        self.assertNotEqual(r.returncode, 0)
        self.assertFalse((self.tmp / "out.tsv").exists())


if __name__ == "__main__":
    unittest.main()
