"""Known-answer tests for the truth builder's N rule and for the fast exact search."""
import random
import unittest
from pathlib import Path

import _paths
import truth_match as T


def rand_seq(rng, n, alphabet="ACGT"):
    return "".join(rng.choice(alphabet) for _ in range(n))


def brute(query, subjects, n_wildcard, max_mm=3):
    """v1 all-pairs search (numpy.convolve, both strands), as process_mock did it."""
    res = {}
    for k, s in subjects.items():
        m = min(T.brute_force_min(query, s, n_wildcard),
                T.brute_force_min(T.revcomp(query), s, n_wildcard))
        res[k] = m
    m = min(res.values())
    if m > max_mm:
        return None, []
    return m, sorted(k for k, v in res.items() if v == m)


class TestNRule(unittest.TestCase):
    def test_long_n_run_ties_every_asv_only_in_wildcard_mode(self):
        rng = random.Random(1)
        asv = rand_seq(rng, 60)
        true_src = rand_seq(rng, 40) + asv + rand_seq(rng, 40)
        n_src = rand_seq(rng, 50) + "N" * 80 + rand_seq(rng, 50)   # unrelated, N run > ASV
        subj = {"true": true_src, "nrun": n_src}
        self.assertEqual(T.best_hits(asv, subj, n_wildcard=True), (0, ["nrun", "true"]))
        self.assertEqual(T.best_hits(asv, subj, n_wildcard=False), (0, ["true"]))
        tax = {"true": "d__B; g__X; s__X a", "nrun": "d__B; g__X; s__X b"}
        self.assertEqual(T.trueish_label(*T.best_hits(asv, subj, n_wildcard=True), tax), "d__B; g__X")
        self.assertEqual(T.trueish_label(*T.best_hits(asv, subj, n_wildcard=False), tax), "d__B; g__X; s__X a")

    def test_short_n_run_lowers_mismatches_in_wildcard_mode(self):
        rng = random.Random(2)
        asv = rand_seq(rng, 50)
        mutated = list(asv)
        for i in (5, 20, 35, 45):          # 4 mismatches -> 'other' unless N rescues them
            mutated[i] = "A" if mutated[i] != "A" else "C"
        src = "".join(mutated)
        src_n = src[:4] + "NN" + src[6:]   # covers mismatch at 5
        self.assertEqual(T.best_hits(asv, {"s": src}, n_wildcard=False), (None, []))
        self.assertEqual(T.best_hits(asv, {"s": src_n}, n_wildcard=False), (None, []))
        self.assertEqual(T.best_hits(asv, {"s": src_n}, n_wildcard=True), (3, ["s"]))

    def test_v1_functions_match_convolution(self):
        import numpy
        q = numpy.array(list("ACGTA"))
        s = numpy.array(list("TTACNTATT"))
        self.assertEqual(T.gapless_mismatches(q, s, n_wildcard=True), 0)   # ACNTA: N matches G
        self.assertEqual(T.gapless_mismatches(q, s, n_wildcard=False), 1)


class TestFastSearchEqualsBruteForce(unittest.TestCase):
    def test_random(self):
        rng = random.Random(20261029)
        for trial in range(120):
            L = rng.randint(20, 70)
            core = rand_seq(rng, L)
            subjects = {}
            for j in range(rng.randint(2, 5)):
                s = list(rand_seq(rng, rng.randint(0, 30)) + core + rand_seq(rng, rng.randint(0, 30)))
                # 0-5 substitutions, sometimes an N run, sometimes truncated end (partial overlap)
                for _ in range(rng.randint(0, 5)):
                    i = rng.randrange(len(s))
                    s[i] = rng.choice("ACGTN")
                if rng.random() < 0.3:
                    i = rng.randrange(len(s))
                    s[i:i + rng.randint(1, L + 5)] = ["N"] * rng.randint(1, L + 5)
                s = "".join(s)
                if rng.random() < 0.2:
                    s = s[rng.randint(1, 3):]
                if rng.random() < 0.3:
                    s = T.revcomp(s)
                subjects[f"s{j}"] = s
            for nw in (True, False):
                self.assertEqual(T.best_hits(core, subjects, n_wildcard=nw), brute(core, subjects, nw),
                                 f"trial {trial} n_wildcard={nw}")


class TestFFTPath(unittest.TestCase):
    def test_random_fft_equals_brute_force(self):
        old = T.FFT_THRESHOLD
        T.FFT_THRESHOLD = 0
        try:
            TestFastSearchEqualsBruteForce("test_random").test_random()
        finally:
            T.FFT_THRESHOLD = old


class TestGutV5V7KnownAnswer(unittest.TestCase):
    """The rebuild with N as a mismatch must reproduce the 2026-10-06 corrected
    truth exactly, and wildcard mode must reproduce the published truth."""
    # Needs the generated data tree (DATA.md); set TAXOBENCH_DATA_DIR or the
    # data_dir of configs/pipeline_config.json. The corrected truth is the
    # regression fixture.
    import _common
    TB = Path(_common.load_pipeline_config()["data_dir"]) / "mock_community"
    CORR = _paths.FIXTURES / "gutmock_V5-V7/data/mock_community/gutmock/V5-V7/trueish_taxonomy_gutmock.tsv"

    @unittest.skipUnless((TB / "input_data/gutmock_seqs.fasta").exists()
                         and (TB / "gutmock/V5-V7/trueish_taxonomy_gutmock.tsv").exists(),
                         "taxobench data tree (gutmock V5-V7) not found")
    def test_rebuild(self):
        import rebuild_truth as R
        rep = self.TB / "gutmock/V5-V7/rep_seqs_gutmock.fasta"
        if not rep.exists():
            rep = next((self.TB / "gutmock/V5-V7").rglob("rep_seqs_gutmock.fasta"))
        src = self.TB / "input_data/gutmock_seqs.fasta"
        tax = self.TB / "input_data/gutmock_taxa.txt"
        strict = R.rebuild(rep, src, tax, n_wildcard=False)
        wild = R.rebuild(rep, src, tax, n_wildcard=True)
        self.assertEqual(strict, R.read_truth(self.CORR))
        self.assertEqual(wild, R.read_truth(self.TB / "gutmock/V5-V7/trueish_taxonomy_gutmock.tsv"))


if __name__ == "__main__":
    unittest.main()
