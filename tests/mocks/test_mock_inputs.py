"""The published mock inputs in mocks/inputs/ are complete and intact."""
import gzip
import hashlib
import json
import unittest

import _paths

ROOT = _paths.ROOT
INPUTS = ROOT / "mocks" / "inputs"


class TestMockInputs(unittest.TestCase):

    def setUp(self):
        self.sums = dict(reversed(line.split()) for line in
                         (INPUTS / "MD5SUMS").read_text().splitlines())

    def test_every_matrix_mock_has_inputs(self):
        mocks = json.loads((ROOT / "mocks" / "configs" / "mocks.json").read_text())
        matrix = mocks["matrix"]
        for mock in matrix["v_mocks"] + matrix["rrn_mocks"]:
            for suffix in ("seqs.fasta", "taxa.txt"):
                self.assertIn(f"{mock}_{suffix}", self.sums, mock)

    def test_checksums(self):
        for name, digest in self.sums.items():
            data = gzip.decompress((INPUTS / f"{name}.gz").read_bytes())
            self.assertEqual(hashlib.md5(data).hexdigest(), digest, name)


if __name__ == "__main__":
    unittest.main()
