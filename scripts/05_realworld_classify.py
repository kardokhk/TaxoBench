#!/usr/bin/env python3
"""Step 05 - Run all 8 classifiers on the 3 real-world FASTAs.

No ground truth -> no evaluation, only `taxonomy_final.tsv` per
(study x classifier x param).
"""

import argparse
import json
import subprocess
from itertools import product
from os.path import join, exists
from pathlib import Path

from joblib import Parallel, delayed

import _common


SUPPORTED = ("taxonova", "vsearch", "blast", "naive-bayes",
             "idtaxa", "spingo", "kraken2", "mapseq")


def _expand(params):
    keys, vals = list(params.keys()), list(params.values())
    vals = [v if isinstance(v, list) else [v] for v in vals]
    return [dict(zip(keys, c)) for c in product(*vals)] or [{}]


def _ref_dict(classifier, ref_root, region):
    # Same logic as 03_assign — duplicated to keep this module standalone.
    region_dir = join(ref_root, region)
    if not Path(region_dir).is_dir():
        return None
    if region == "RRN":
        seqs = join(region_dir, "grond_clean.fasta")
        taxa = join(region_dir, "grond_clean.tsv")
        merged = seqs
    else:
        seqs = join(region_dir, "sequences.fasta")
        taxa = join(region_dir, "taxonomy.tsv")
        merged = join(region_dir, "merged.fasta")
        if not exists(merged):
            merged = None
    base = {"seqs": seqs, "taxa": taxa, "merged": merged}

    import glob
    if classifier == "taxonova":
        m = glob.glob(join(region_dir, "taxonova", "*model.pt"))
        base["model"] = m[0] if m else None
    elif classifier == "naive-bayes":
        m = glob.glob(join(region_dir, "naive-bayes", "*sklearn*qza"))
        base["model"] = m[0] if m else None
    elif classifier == "idtaxa":
        m = (glob.glob(join(region_dir, "idtaxa", "*.rds")) or
             glob.glob(join(region_dir, "idtaxa", "*.RData")))
        base["model"] = m[0] if m else None
    elif classifier == "spingo":
        m = glob.glob(join(region_dir, "spingo", "*.fasta"))
        base["database"] = m[0] if m else None
    elif classifier == "kraken2":
        kdir = join(region_dir, "kraken2")
        ok = all(exists(join(kdir, f)) for f in ("hash.k2d", "opts.k2d", "taxo.k2d"))
        base["database"] = kdir if ok else None
    elif classifier == "mapseq":
        mdir = join(region_dir, "mapseq")
        fa = join(mdir, "reference.fasta")
        tx = join(mdir, "reference.tax")
        base["database_fasta"] = fa if exists(fa) else None
        base["database_tax"]   = tx if exists(tx) else None
    return base


def _run(cmd, log_fp):
    Path(log_fp).parent.mkdir(parents=True, exist_ok=True)
    with open(log_fp, "w") as fh:
        ret = subprocess.run(cmd, shell=True, stdout=fh, stderr=subprocess.STDOUT)
    return cmd, ret.returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--studies", nargs="+",
                    default=["zuo", "mantri", "virtanen"])
    ap.add_argument("--region", default=None,
                    help="override per-study region from mocks.json")
    ap.add_argument("--classifiers", nargs="+", default=list(SUPPORTED))
    ap.add_argument("--params-config", default=None)
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--n-jobs", type=int, default=1)
    ap.add_argument("--results-dir", default=None,
                    help="default: results/real_world")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    cfg = _common.load_pipeline_config()
    classifiers_cfg = _common.load_classifiers_config()
    if args.params_config:
        with open(args.params_config) as fh:
            classifiers_cfg = json.load(fh)
    mocks_cfg = _common.load_mocks_config()

    bad = [c for c in args.classifiers if c not in SUPPORTED]
    if bad:
        ap.error(f"Unsupported classifiers: {bad}")

    threads = args.threads or cfg.get("default_n_threads", 4)
    ref_root = cfg["reference_database_dir"]
    src_dir = cfg["scripts_dir"]
    rw_input = Path(cfg["real_world_input_dir"])
    res_root = Path(args.results_dir) if args.results_dir \
        else Path(cfg["results_dir"]) / "real_world"

    _common.add_src_to_path()
    import generate_commands  # noqa

    params_all = classifiers_cfg["method_parameters_combinations"]
    envs_all   = classifiers_cfg["method_conda_envs"]

    builders = {
        "taxonova":    generate_commands.build_taxonova_command,
        "vsearch":     generate_commands.build_vsearch_command,
        "blast":       generate_commands.build_blast_command,
        "naive-bayes": generate_commands.build_naive_bayes_command,
        "idtaxa":      generate_commands.build_idtaxa_command,
        "spingo":      generate_commands.build_spingo_command,
        "kraken2":     generate_commands.build_kraken2_command,
        "mapseq":      generate_commands.build_mapseq_command,
    }

    jobs = []
    for study in args.studies:
        meta = mocks_cfg["real_world"].get(study)
        if not meta:
            print(f"[05_realworld] unknown study {study}; skip")
            continue
        fasta = rw_input / meta["fasta"]
        if not fasta.exists():
            print(f"[05_realworld] missing fasta {fasta}; skip {study}")
            continue
        rep_base = str(fasta.with_suffix(""))
        region = args.region or meta.get("region", "V3-V4")

        for classifier in args.classifiers:
            ref_dict = _ref_dict(classifier, ref_root, region)
            if ref_dict is None:
                print(f"[05_realworld] skip {classifier}/{region}: no ref")
                continue
            param_grid = _expand(params_all.get(classifier, {}))
            conda_env = envs_all[classifier]
            for param_dict in param_grid:
                param_id = "_".join(f"{k}-{v}" for k, v in sorted(param_dict.items())) or "default"
                outdir = res_root / study / classifier / param_id
                outdir.mkdir(parents=True, exist_ok=True)
                kw = dict(rep_seq_base=rep_base, outdir=str(outdir),
                          param_dict=param_dict, reference_dict=ref_dict,
                          conda_env=conda_env, dry_run=args.dry_run)
                if classifier in ("taxonova", "idtaxa", "spingo", "kraken2", "mapseq"):
                    kw["scripts_dir"] = src_dir
                if classifier in ("vsearch", "blast", "naive-bayes",
                                  "spingo", "kraken2", "mapseq"):
                    kw["n_threads"] = threads
                try:
                    cmd = builders[classifier](**kw)
                except Exception as exc:
                    print(f"  ! build {classifier} failed: {exc}")
                    continue
                jobs.append((cmd, str(outdir / "run.log"), classifier, study))

    print(f"[05_realworld] queued {len(jobs)} jobs")
    if args.dry_run:
        for cmd, log, c, s in jobs:
            print(f"  [dry] {s} | {c}\n        {cmd}")
        return

    results = Parallel(n_jobs=args.n_jobs, verbose=10)(
        delayed(_run)(cmd, log) for cmd, log, *_ in jobs
    )
    fails = [r for r in results if r[1] != 0]
    print(f"[05_realworld] done. {len(jobs) - len(fails)} ok / {len(fails)} failed.")


if __name__ == "__main__":
    main()
