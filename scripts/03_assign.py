#!/usr/bin/env python3
"""Step 03 - Assign taxonomy with the 8 supported classifiers.

For each (mock x region x classifier x param-combo), build the command
via `generate_commands` and run them in parallel via joblib.
"""

import argparse
import os
import subprocess
from itertools import product
from os.path import join, exists
from pathlib import Path

from joblib import Parallel, delayed

import _common


SUPPORTED = ("taxonova", "vsearch", "blast", "naive-bayes",
             "idtaxa", "spingo", "kraken2", "mapseq")


def _expand_param_grid(params):
    """{a: [1,2], b: [3]} -> [{a:1,b:3}, {a:2,b:3}]"""
    keys, vals = list(params.keys()), list(params.values())
    vals = [v if isinstance(v, list) else [v] for v in vals]
    return [dict(zip(keys, c)) for c in product(*vals)]


def _run_single(cmd, log_fp):
    Path(os.path.dirname(log_fp)).mkdir(parents=True, exist_ok=True)
    with open(log_fp, "w") as fh:
        ret = subprocess.run(cmd, shell=True, stdout=fh,
                             stderr=subprocess.STDOUT)
    return cmd, ret.returncode


def main():
    ap = argparse.ArgumentParser(description=(
        "Run the classifiers on one mock: one job per region, classifier and "
        "parameter set. --dry-run prints the commands without running them."))
    ap.add_argument("--mock-name", required=True)
    ap.add_argument("--regions", nargs="+", default=["V4", "V3-V4", "V5-V7", "V1-V9", "RRN"])
    ap.add_argument("--classifiers", nargs="+", default=list(SUPPORTED))
    ap.add_argument("--params-config", default=None,
                    help="Classifier config: a JSON path or a name under configs/ "
                         "(default classifiers_default.json; classifiers_paper.json "
                         "pins the published Taxonova runs)")
    ap.add_argument("--ref-root", default=None,
                    help="Reference database root (default: reference_database_dir "
                         "of pipeline_config.json)")
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--n-jobs", type=int, default=1)
    ap.add_argument("--mock-data-dir", default=None,
                    help="Where rep_seqs files live; default data/mock_community/{mock}")
    ap.add_argument("--results-dir", default=None,
                    help="Output dir; default results/{mock}/")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    cfg = _common.load_pipeline_config()
    classifiers_cfg = _common.load_classifiers_config(args.params_config)

    bad = [c for c in args.classifiers if c not in SUPPORTED]
    if bad:
        ap.error(f"Unsupported classifiers: {bad}. Allowed: {SUPPORTED}")

    threads = args.threads or cfg.get("default_n_threads", 4)
    ref_root = args.ref_root or cfg["reference_database_dir"]
    src_dir = cfg["scripts_dir"]
    mock_data_dir = Path(args.mock_data_dir) if args.mock_data_dir \
        else Path(cfg["data_dir"]) / "mock_community" / args.mock_name
    results_root = Path(args.results_dir) if args.results_dir \
        else Path(cfg["results_dir"]) / args.mock_name

    _common.add_src_to_path()
    import generate_commands  # noqa
    tx_settings = generate_commands.taxonova_settings(classifiers_cfg)

    params_all = classifiers_cfg["method_parameters_combinations"]
    envs_all   = classifiers_cfg["method_conda_envs"]

    builders = {
        "taxonova":    lambda **kw: generate_commands.build_taxonova_command(**kw),
        "vsearch":     lambda **kw: generate_commands.build_vsearch_command(**kw),
        "blast":       lambda **kw: generate_commands.build_blast_command(**kw),
        "naive-bayes": lambda **kw: generate_commands.build_naive_bayes_command(**kw),
        "idtaxa":      lambda **kw: generate_commands.build_idtaxa_command(**kw),
        "spingo":      lambda **kw: generate_commands.build_spingo_command(**kw),
        "kraken2":     lambda **kw: generate_commands.build_kraken2_command(**kw),
        "mapseq":      lambda **kw: generate_commands.build_mapseq_command(**kw),
    }

    import glob
    jobs = []
    for region in args.regions:
        # Match the original convention: mockdir/{region}/*/rep_seqs_{mock}.fasta
        # (the '*' is typically 'source/single' or 'source/paired'). For RRN,
        # the file may be at mockdir/RRN/rep_seqs_{mock}.fasta directly.
        candidates = (
            glob.glob(str(mock_data_dir / region / "*" / f"rep_seqs_{args.mock_name}.fasta")) +
            glob.glob(str(mock_data_dir / region / "**" / f"rep_seqs_{args.mock_name}.fasta"), recursive=True) +
            glob.glob(str(mock_data_dir / region / f"rep_seqs_{args.mock_name}.fasta"))
        )
        candidates = sorted(set(candidates))
        if not candidates:
            print(f"[03_assign] skip {region}: no rep_seqs_{args.mock_name}.fasta under {mock_data_dir/region}")
            continue
        rep_base = Path(candidates[0]).with_suffix("")

        for classifier in args.classifiers:
            ref_dict = _common.ref_dict(classifier, ref_root, region, tx_settings)
            if ref_dict is None:
                print(f"[03_assign] skip {classifier}/{region}: no reference dir")
                continue
            param_grid = _expand_param_grid(params_all.get(classifier, {})) or [{}]
            conda_env = envs_all[classifier]
            for param_dict in param_grid:
                param_id = "_".join(f"{k}-{v}" for k, v in sorted(param_dict.items())) or "default"
                outdir = results_root / region / classifier / param_id
                if not args.dry_run:
                    outdir.mkdir(parents=True, exist_ok=True)
                build = builders[classifier]
                kw = dict(rep_seq_base=str(rep_base), outdir=str(outdir),
                          param_dict=param_dict, reference_dict=ref_dict,
                          conda_env=conda_env, dry_run=args.dry_run)
                if classifier in ("taxonova", "idtaxa", "spingo", "kraken2", "mapseq"):
                    kw["scripts_dir"] = src_dir
                if classifier == "taxonova":
                    kw["settings"] = tx_settings
                if classifier in ("vsearch", "blast", "naive-bayes",
                                  "spingo", "kraken2", "mapseq"):
                    kw["n_threads"] = threads
                try:
                    cmd = build(**kw)
                except Exception as exc:
                    print(f"  ! build {classifier}/{region} failed: {exc}")
                    continue
                log_fp = outdir / "run.log"
                jobs.append((cmd, str(log_fp), classifier, region, param_id))

    print(f"[03_assign] queued {len(jobs)} jobs")
    if args.dry_run:
        for cmd, log, c, r, p in jobs:
            print(f"  [dry] {c} | {r} | {p}\n        {cmd}")
        return

    results = Parallel(n_jobs=args.n_jobs, verbose=10)(
        delayed(_run_single)(cmd, log) for cmd, log, *_ in jobs
    )
    fails = [r for r in results if r[1] != 0]
    print(f"[03_assign] done. {len(jobs) - len(fails)} ok / {len(fails)} failed.")


if __name__ == "__main__":
    main()
