#!/usr/bin/env python3
"""Step 06 - Runtime / memory benchmark for the 8 supported classifiers.

Samples N=[100, 500, 1000, 2000, 5000] sequences from a representative
mock fasta, runs each classifier `--replicates` times, captures wall
time and peak RSS via `psutil`, and writes one CSV per classifier under
results/runtime/{classifier}_runtime.csv with columns:
    n_seqs, replicate, wall_seconds, peak_rss_mb, region
"""

import argparse
from pathlib import Path

import _common


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-fasta", required=False, default=None,
                    help="Representative fasta to sample from")
    ap.add_argument("--mock-name", default="mock-3",
                    help="Mock to use as source if --source-fasta omitted")
    ap.add_argument("--region", default="V4")
    ap.add_argument("--regions", nargs="+", default=None,
                    help="If given, overrides --region (multi-region run)")
    ap.add_argument("--sizes", type=int, nargs="+",
                    default=[100, 500, 1000, 2000, 5000])
    ap.add_argument("--classifiers", nargs="+",
                    default=["taxonova", "vsearch", "blast", "naive-bayes",
                             "idtaxa", "spingo", "kraken2", "mapseq"])
    ap.add_argument("--replicates", type=int, default=3)
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--results-dir", default=None,
                    help="default: results/runtime")
    ap.add_argument("--params-config", default=None)
    args = ap.parse_args()

    cfg = _common.load_pipeline_config()
    classifiers_cfg = _common.load_classifiers_config(args.params_config)

    threads = args.threads or cfg.get("default_n_threads", 4)
    out_dir = Path(args.results_dir) if args.results_dir \
        else Path(cfg["results_dir"]) / "runtime"
    out_dir.mkdir(parents=True, exist_ok=True)

    src_fasta = args.source_fasta
    if src_fasta is None:
        # Try the standard mock layout from step 02
        candidate = Path(cfg["data_dir"]) / "mock_community" / args.mock_name / args.region / "single" / f"rep_seqs_{args.mock_name}.fasta"
        if candidate.exists():
            src_fasta = str(candidate)
        else:
            candidate = Path(cfg["mock_input_dir"]) / f"{args.mock_name}_seqs.fasta"
            if candidate.exists():
                src_fasta = str(candidate)
    if src_fasta is None or not Path(src_fasta).exists():
        ap.error(f"Could not locate source fasta. Pass --source-fasta. Tried: {candidate}")

    regions = args.regions or [args.region]

    flat_envs = {
        k: (v if isinstance(v, str) else v) for k, v in
        classifiers_cfg["method_conda_envs"].items()
    }

    bench_config = {
        "project_dir":            cfg["project_dir"],
        "reference_database_dir": cfg["reference_database_dir"],
        "scripts_dir":            cfg["scripts_dir"],
        "output_dir":             str(out_dir),
        "source_fasta":           src_fasta,
        "classifiers":            args.classifiers,
        "sequence_sizes":         args.sizes,
        "v_regions":              regions,
        "n_replicates":           args.replicates,
        "n_threads":              threads,
        "conda_envs":             flat_envs,
        "taxonova_settings":      classifiers_cfg.get("taxonova", {}),
        "classifier_params": {
            c: {k: (v[0] if isinstance(v, list) else v)
                for k, v in classifiers_cfg["method_parameters_combinations"].get(c, {}).items()}
            for c in args.classifiers
        },
    }

    _common.add_src_to_path()
    from benchmark_classifiers import PerformanceBenchmarker
    bench = PerformanceBenchmarker(bench_config)
    bench.run_benchmarks()
    bench.save_results()
    print(f"[06_runtime] done; outputs in {out_dir}")


if __name__ == "__main__":
    main()
