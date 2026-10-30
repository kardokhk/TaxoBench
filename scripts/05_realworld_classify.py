#!/usr/bin/env python3
"""Step 05 - Run all 8 classifiers on the 3 real-world FASTAs.

No ground truth -> no evaluation, only `taxonomy_final.tsv` per
(study x classifier x param).
"""

import argparse
import subprocess
import sys
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


def _run(cmd, log_fp):
    Path(log_fp).parent.mkdir(parents=True, exist_ok=True)
    with open(log_fp, "w") as fh:
        ret = subprocess.run(cmd, shell=True, stdout=fh, stderr=subprocess.STDOUT)
    return cmd, ret.returncode


def main():
    ap = argparse.ArgumentParser(description=(
        "Run the classifiers on the real-world ASV FASTAs (no truth; writes "
        "taxonomy_final.tsv per study, classifier and parameter set). "
        "--dry-run prints the commands without running them."))
    ap.add_argument("--studies", nargs="+",
                    default=["zuo", "mantri", "virtanen"])
    ap.add_argument("--region", default=None,
                    help="override per-study region from mocks.json")
    ap.add_argument("--classifiers", nargs="+", default=list(SUPPORTED))
    ap.add_argument("--params-config", default=None,
                    help="Classifier config: a JSON path or a name under configs/ "
                         "(default classifiers_default.json)")
    ap.add_argument("--ref-root", default=None,
                    help="Reference database root (default: reference_database_dir "
                         "of pipeline_config.json)")
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--n-jobs", type=int, default=1)
    ap.add_argument("--results-dir", default=None,
                    help="default: results/real_world")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    cfg = _common.load_pipeline_config()
    classifiers_cfg = _common.load_classifiers_config(args.params_config)
    mocks_cfg = _common.load_mocks_config()

    bad = [c for c in args.classifiers if c not in SUPPORTED]
    if bad:
        ap.error(f"Unsupported classifiers: {bad}")

    threads = args.threads or cfg.get("default_n_threads", 4)
    ref_root = args.ref_root or cfg["reference_database_dir"]
    src_dir = cfg["scripts_dir"]
    rw_input = Path(cfg["real_world_input_dir"])
    res_root = Path(args.results_dir) if args.results_dir \
        else Path(cfg["results_dir"]) / "real_world"

    _common.add_src_to_path()
    import generate_commands  # noqa
    tx_settings = generate_commands.taxonova_settings(classifiers_cfg)

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
    study_fasta = {}
    for study in args.studies:
        meta = mocks_cfg["real_world"].get(study)
        if not meta:
            print(f"[05_realworld] unknown study {study}; skip")
            continue
        fasta = rw_input / meta["fasta"]
        if not fasta.exists():
            print(f"[05_realworld] missing fasta {fasta}; skip {study}")
            continue
        study_fasta[study] = fasta
        rep_base = str(fasta.with_suffix(""))
        region = args.region or meta.get("region", "V3-V4")

        for classifier in args.classifiers:
            ref_dict = _common.ref_dict(classifier, ref_root, region, tx_settings)
            if ref_dict is None:
                print(f"[05_realworld] skip {classifier}/{region}: no ref")
                continue
            param_grid = _expand(params_all.get(classifier, {}))
            conda_env = envs_all[classifier]
            for param_dict in param_grid:
                param_id = "_".join(f"{k}-{v}" for k, v in sorted(param_dict.items())) or "default"
                outdir = res_root / study / classifier / param_id
                if not args.dry_run:
                    outdir.mkdir(parents=True, exist_ok=True)
                kw = dict(rep_seq_base=rep_base, outdir=str(outdir),
                          param_dict=param_dict, reference_dict=ref_dict,
                          conda_env=conda_env, dry_run=args.dry_run)
                if classifier in ("taxonova", "idtaxa", "spingo", "kraken2", "mapseq"):
                    kw["scripts_dir"] = src_dir
                if classifier == "taxonova":
                    kw["settings"] = tx_settings
                    kw["n_threads"] = threads
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

    # Every classifier must have scored exactly the query set it was given.
    # Feature ids are positional (ASV_1..ASV_n), so a stale imported artifact
    # with a different sequence set would otherwise be joined to the counts by
    # id without any error. A row count alone is not enough; the ids must agree.
    mismatches = []
    for _, _, classifier, study in jobs:
        fasta = study_fasta.get(study)
        if fasta is None or not fasta.exists():
            continue
        want = set()
        with open(fasta) as fh:
            for line in fh:
                if line.startswith(">"):
                    want.add(line[1:].split()[0])
        for out in sorted(res_root.glob(f"{study}/{classifier}/*/taxonomy_final.tsv")):
            got = set()
            with open(out) as fh:
                next(fh, None)
                for line in fh:
                    if line.strip():
                        got.add(line.split("\t")[0])
            if got != want:
                mismatches.append(
                    f"{study}/{classifier}: {len(got)} rows against "
                    f"{len(want)} query sequences, {len(want - got)} missing, "
                    f"{len(got - want)} unexpected ({out})")
    if mismatches:
        print("[05_realworld] FATAL: output does not match the query set:",
              file=sys.stderr)
        for m in mismatches:
            print("  " + m, file=sys.stderr)
        sys.exit(1)
    print("[05_realworld] query-set check passed for all cells")


if __name__ == "__main__":
    main()
