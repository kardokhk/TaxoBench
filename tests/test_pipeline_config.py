"""Configs and command building: no data, no QIIME 2, no classifiers needed.

The paper configuration must keep emitting the command the published runs
used; the default configuration must use the current Taxonova release.
"""
import hashlib
import json
import py_compile
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import _paths  # noqa: F401
import _common
import generate_commands as G

ROOT = _paths.ROOT
CONFIGS = ROOT / "configs"
PANEL = ["taxonova", "vsearch", "blast", "naive-bayes",
         "idtaxa", "spingo", "kraken2", "mapseq"]
REGIONS = ["V4", "V3-V4", "V5-V7", "V1-V9", "RRN"]


def _taxonova_cmd(config_name, region="V4", ref_root="/ref"):
    cfg = _common.load_classifiers_config(config_name)
    settings = G.taxonova_settings(cfg)
    params = {k: v[0] for k, v in
              cfg["method_parameters_combinations"]["taxonova"].items()}
    model = G.resolve_taxonova_model(region, f"{ref_root}/{region}", settings,
                                     ref_root=ref_root)
    ref = {"seqs": f"{ref_root}/{region}/sequences.fasta", "model": model}
    return G.build_taxonova_command(
        rep_seq_base="/in/rep_seqs_m", outdir="/out", param_dict=params,
        reference_dict=ref, scripts_dir="/src", conda_env="/env",
        dry_run=True, settings=settings)


class TestClassifierConfigs(unittest.TestCase):

    def test_configs_are_complete(self):
        named_envs = _common.load_pipeline_config()["conda_envs"]
        for name in ("classifiers_default.json", "classifiers_paper.json",
                     "classifiers_sweep.json"):
            with self.subTest(config=name):
                raw = json.loads((CONFIGS / name).read_text())
                self.assertEqual(raw["classifiers"], PANEL)
                for c in PANEL:
                    self.assertIn(c, raw["method_parameters_combinations"])
                    env = raw["method_conda_envs"][c]
                    keys = env.values() if isinstance(env, dict) else [env]
                    for k in keys:
                        self.assertIn(k, named_envs, f"{name}: {c} env {k!r}")
                cfg = _common.load_classifiers_config(name)
                self.assertEqual(cfg["method_conda_envs"]["blast"],
                                 named_envs["qiime2"])

    def test_default_uses_current_release(self):
        cfg = _common.load_classifiers_config()
        tx = cfg["taxonova"]
        self.assertEqual(tx["cli"], "classify")
        self.assertEqual(tx["requirement"], "taxonova>=26.10.1")
        self.assertEqual(tx["models"], {"V4": "v4", "V3-V4": "v3v4", "V5-V7": "v5v7",
                                        "V1-V9": "16s", "RRN": "operon"})
        self.assertIsNone(tx["parent_gates"])
        self.assertEqual(cfg["method_parameters_combinations"]["taxonova"],
                         {"confidence_threshold": [0.8]})

    def test_paper_pins_published_runs(self):
        cfg = _common.load_classifiers_config("classifiers_paper.json")
        tx = cfg["taxonova"]
        self.assertEqual(tx["cli"], "classification")
        self.assertEqual(tx["scoring_rules"], "v1")
        self.assertIn("@v1.0.0", tx["requirement"])
        self.assertEqual(tx["parent_gates"],
                         {r: 0.5 for r in G.TAXONOVA_PARENT_RANKS})
        self.assertEqual(cfg["method_parameters_combinations"]["taxonova"],
                         {"confidence_threshold": [0.7]})
        self.assertEqual(sorted(tx["models"]), sorted(REGIONS))
        for region, spec in tx["models"].items():
            self.assertTrue(spec["path"].startswith(f"{region}/taxonova/"))
            self.assertIn("/releases/download/v1.0.0/", spec["url"])
            self.assertEqual(len(spec["md5"]), 32)
        # Every other classifier is configured as in the default panel.
        default = _common.load_classifiers_config()
        for c in PANEL[1:]:
            self.assertEqual(cfg["method_parameters_combinations"][c],
                             default["method_parameters_combinations"][c])


class TestTaxonovaCommand(unittest.TestCase):

    def test_paper_command_is_the_published_one(self):
        # The command the published runs used (b609765), parent gates first.
        expected = (
            'conda run -p /env bash -c "mkdir -p /out && taxonova classification '
            '--domain_threshold 0.5 --phylum_threshold 0.5 --class_threshold 0.5 '
            '--order_threshold 0.5 --family_threshold 0.5 --genus_threshold 0.5 '
            '--confidence_threshold 0.7 --input /in/rep_seqs_m.fasta '
            '--model_input /ref/V4/taxonova/gg2_v4_mrgroups_model.pt '
            '--output /out/taxonomy_raw.tsv && python /src/process_taxonova_output.py '
            '/out/taxonomy_raw.tsv /out/taxonomy_final.tsv False "')
        self.assertEqual(_taxonova_cmd("classifiers_paper.json"), expected)

    def test_default_command(self):
        for region, model_id in (("V4", "v4"), ("V1-V9", "16s"), ("RRN", "operon")):
            cmd = _taxonova_cmd("classifiers_default.json", region)
            self.assertIn(f"taxonova classify /in/rep_seqs_m.fasta --model {model_id} "
                          "--confidence 0.8 -o /out/taxonomy_raw.tsv", cmd)
            self.assertNotIn("_threshold", cmd)

    def test_threads_only_with_classify(self):
        for name, has in (("classifiers_default.json", True), ("classifiers_paper.json", False)):
            cfg = _common.load_classifiers_config(name)
            settings = G.taxonova_settings(cfg)
            cmd = G.build_taxonova_command(
                "/in/q", "/out", {"confidence_threshold": 0.8},
                {"seqs": "/r/s.fasta", "model": "/m.pt"}, "/src", "/env", True,
                settings, n_threads=8)
            self.assertEqual("--threads 8" in cmd, has, name)

    def test_param_dict_gate_overrides_config(self):
        settings = {"cli": "classification",
                    "parent_gates": {r: 0.5 for r in G.TAXONOVA_PARENT_RANKS}}
        cmd = G.build_taxonova_command(
            "/in/q", "/out", {"confidence_threshold": 0.0, "genus_threshold": 0.0},
            {"seqs": "/r/s.fasta", "model": "/m.pt"}, "/src", "/env", True, settings)
        self.assertIn("--genus_threshold 0.0 --confidence_threshold 0.0", cmd)
        self.assertEqual(cmd.count("--genus_threshold"), 1)

    def test_bad_settings_rejected(self):
        with self.assertRaises(ValueError):
            G.build_taxonova_command("/q", "/o", {}, {"seqs": "/s", "model": "v4"},
                                     "/src", "/env", True, {"cli": "predict"})
        with self.assertRaises(ValueError):
            G.build_taxonova_command("/q", "/o", {}, {"seqs": "/s", "model": "v4"},
                                     "/src", "/env", True,
                                     {"parent_gates": {"kingdom": 0.5}})

    def test_model_resolution(self):
        with tempfile.TemporaryDirectory() as tmp:
            region_dir = Path(tmp) / "V4"
            (region_dir / "taxonova").mkdir(parents=True)
            for n in ("b_model.pt", "a_model.pt"):
                (region_dir / "taxonova" / n).write_bytes(b"x")
            # No configured entry: first *model.pt in sorted order.
            self.assertEqual(G.resolve_taxonova_model("V4", str(region_dir), {}),
                             str(region_dir / "taxonova" / "a_model.pt"))
            # Model id: passed through.
            self.assertEqual(G.resolve_taxonova_model(
                "V4", str(region_dir), {"models": {"V4": "v4"}}), "v4")
            # Pinned file: md5 checked.
            good = hashlib.md5(b"x").hexdigest()
            spec = {"models": {"V4": {"path": "V4/taxonova/a_model.pt", "md5": good}}}
            self.assertEqual(G.resolve_taxonova_model("V4", str(region_dir), spec,
                                                      ref_root=tmp),
                             str(Path(tmp) / "V4/taxonova/a_model.pt"))
            spec["models"]["V4"]["path"] = "V4/taxonova/b_model.pt"
            spec["models"]["V4"]["md5"] = "0" * 32
            with self.assertRaises(ValueError):
                G.resolve_taxonova_model("V4", str(region_dir), spec, ref_root=tmp)


class TestScriptsLoad(unittest.TestCase):
    """Every script compiles; the orchestration scripts answer --help."""

    def test_compile_all(self):
        for p in sorted(list((ROOT / "scripts").rglob("*.py")) +
                        list((ROOT / "src").rglob("*.py")) +
                        list((ROOT / "mgnify").glob("*.py"))):
            with self.subTest(script=str(p.relative_to(ROOT))):
                py_compile.compile(str(p), doraise=True,
                                   cfile=str(Path(tempfile.gettempdir()) / "tb_compile.pyc"))

    def test_help(self):
        scripts = ["03_assign.py", "04_evaluate.py", "05_realworld_classify.py",
                   "loco/01_make_holdout.py", "loco/02_train_classifier.py",
                   "loco/03_classify.py", "loco/04_evaluate.py",
                   "loco/04c_evaluate_opengates.py",
                   "02d_generate_longread_biomes.py", "02c_generate_v1v9_direct.py",
                   "build_tau_grid.py"]
        for s in scripts:
            with self.subTest(script=s):
                r = subprocess.run([sys.executable, str(ROOT / "scripts" / s), "--help"],
                                   capture_output=True, text=True)
                self.assertEqual(r.returncode, 0, r.stderr[-500:])


if __name__ == "__main__":
    unittest.main()
