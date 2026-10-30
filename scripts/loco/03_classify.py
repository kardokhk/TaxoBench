#!/usr/bin/env python3
"""LOCO step 03: classify novel + control queries with a single classifier.

Outputs (per region × classifier):
    results/loco/<region>/<classifier>/novel_predictions.tsv
    results/loco/<region>/<classifier>/control_predictions.tsv

Uses the held-out reference DB built by step 02. We classify with confidence
threshold = 0 (i.e. always emit a call) so we can sweep thresholds in the
evaluation step. Per-rank confidences are preserved in the TSVs.
"""
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import _common  # noqa


#: The conda executable used for `conda run`. `conda run` is needed to
#: activate the qiime2 env properly for R and QIIME 2 (PATH alone is not
#: enough). Set TAXOBENCH_CONDA_BIN or "conda_bin" in pipeline_config.json
#: when conda is not on PATH (for example on some HPC compute nodes).
CONDA_BIN = _common.conda_bin()

VALID = ("taxonova", "idtaxa", "naive-bayes", "blast", "vsearch", "spingo",
         "kraken2", "mapseq")


def _run(cmd):
    print("[loco/03]", cmd, flush=True)
    subprocess.run(cmd, shell=True, check=True)


def _region_dir(args):
    """Region directory name, optionally tagged so replicates coexist.

    An empty --tag counts as no tag: shell callers pass "${TAG:-}" and an
    unset TAG must not silently produce a 'V4__' directory.
    """
    tag = (args.tag or "").strip()
    return f"{args.region}__{tag}" if tag else args.region


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--region", required=True,
                    choices=["V4", "V3-V4", "V5-V7", "V1-V9", "RRN"])
    ap.add_argument("--classifier", required=True, choices=VALID)
    ap.add_argument("--threads", type=int, default=20)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--loco-root", default="loco",
                    help="Subdirectory of data/ holding the holdout. Default "
                         "'loco' is the shipped 200-species run; the 2,000-species "
                         "replicates live under 'loco2k'.")
    ap.add_argument("--tag", default=None,
                    help="Suffix on the region directory, e.g. 'seed2026'.")
    ap.add_argument("--maxrejects", default=None,
                    help="vsearch only: --p-maxrejects. QIIME's default is "
                         "'all' (vsearch --maxrejects 0), an exhaustive scan of "
                         "the reference for every query. That is affordable on "
                         "short amplicons but not on whole 16S or the full rRNA "
                         "operon: measured at 32 threads, the loco2k V1-V9 and "
                         "RRN cells project to 171 h and 1,830 h. Setting 32, "
                         "vsearch's own default, reproduced the exhaustive V5-V7 "
                         "result exactly on the control split and to within "
                         "0.07%% of species calls on the novel split "
                         "(job 5134571). Left unset, behaviour is unchanged.")
    args = ap.parse_args()

    cfg = _common.load_pipeline_config()
    src_dir = cfg["scripts_dir"]
    env = cfg["conda_envs"]["qiime2"]
    spingo_env = cfg["conda_envs"].get("spingo")
    if isinstance(spingo_env, dict):
        spingo_env = spingo_env.get("spingo_env")

    loco_root = Path(cfg["data_dir"]) / args.loco_root / _region_dir(args)
    ref_root  = loco_root / "retained_ref"
    novel_fa  = loco_root / "novel_queries.fasta"
    ctrl_fa   = loco_root / "control_queries.fasta"
    if not novel_fa.exists():
        sys.exit(f"missing {novel_fa}; run 01_make_holdout.py first")

    out_root = Path(cfg["results_dir"]) / args.loco_root / _region_dir(args) / args.classifier
    out_root.mkdir(parents=True, exist_ok=True)

    novel_out = out_root / "novel_predictions.tsv"
    ctrl_out  = out_root / "control_predictions.tsv"

    # Skip if already done.
    if not args.force and novel_out.exists() and ctrl_out.exists():
        print(f"[loco/03] skip {args.classifier}: outputs exist")
        return

    # ---------- per-classifier dispatch ----------
    if args.classifier == "taxonova":
        model = ref_root / "taxonova" / "model.pt"
        for in_fa, out_tsv in ((novel_fa, novel_out), (ctrl_fa, ctrl_out)):
            raw = out_tsv.with_suffix(".raw.tsv")
            cmd = (
                f"{CONDA_BIN} run -p {env} --no-capture-output bash -c \""
                f"taxonova classification "
                f"--input {in_fa} --model_input {model} --output {raw} "
                # Parent gates stated explicitly. These runs took the
                # package default, a flat 0.5 at every parent rank, and
                # 04_evaluate.py reads their processed output to reproduce the
                # published numbers. taxonova now scales the graded profile
                # with --confidence_threshold, under which 0.0 would open every
                # gate; that is what 04c_evaluate_opengates.py wants but not
                # what these outputs were, so it is pinned rather than changed.
                f"--confidence_threshold 0.0 "
                f"--domain_threshold 0.5 --phylum_threshold 0.5 "
                f"--class_threshold 0.5 --order_threshold 0.5 "
                f"--family_threshold 0.5 --genus_threshold 0.5 && "
                f"python {src_dir}/process_taxonova_output.py "
                f"{raw} {out_tsv} False\""
            )
            _run(cmd)

    elif args.classifier == "idtaxa":
        model = ref_root / "idtaxa" / "trainingSet.RData"
        for in_fa, out_tsv in ((novel_fa, novel_out), (ctrl_fa, ctrl_out)):
            raw = out_tsv.with_suffix(".raw.tsv")
            cmd = (
                f"{CONDA_BIN} run -p {env} --no-capture-output bash -c \""
                f"Rscript {src_dir}/run_idtaxa_classification.R "
                f"{in_fa} {model} {raw} 0\""
            )
            _run(cmd)
            shutil.copy(raw, out_tsv)

    elif args.classifier == "naive-bayes":
        model = ref_root / "naive-bayes" / "classifier.qza"
        helper = Path(src_dir) / "qiime_helper.py"
        for in_fa, out_tsv in ((novel_fa, novel_out), (ctrl_fa, ctrl_out)):
            base = in_fa.with_suffix("")
            # Per-classifier path. blast, vsearch and naive-bayes all
            # imported the queries to the SAME "<base>_qiime.qza" and,
            # run concurrently, truncated each other mid-write. The
            # symptom is misleading: QIIME reports "missing the plugin(s)
            # necessary to load it" or "Archive does not contain a
            # correctly formatted VERSION file", both of which read like
            # a broken install rather than a race.
            qza  = Path(f"{base}_qiime.{args.classifier}.qza")
            cls_qza = out_tsv.with_suffix(".qza")
            cmd_imp = (
                f"{CONDA_BIN} run -p {env} --no-capture-output python {helper} "
                f"import_fasta_to_qza {in_fa} {qza}"
            )
            _run(cmd_imp)
            cmd_run = (
                f"{CONDA_BIN} run -p {env} --no-capture-output bash -c \""
                f"qiime feature-classifier classify-sklearn "
                f"--i-reads {qza} --i-classifier {model} "
                f"--p-n-jobs {args.threads} --p-confidence 0.0 "
                f"--o-classification {cls_qza} && "
                f"python {helper} convert_qza_to_tsv {cls_qza} {out_tsv}\""
            )
            _run(cmd_run)

    elif args.classifier == "blast":
        ref_seqs = ref_root / "blast-vsearch" / "sequences.qza"
        ref_tax  = ref_root / "blast-vsearch" / "taxonomy.qza"
        ref_db   = ref_root / "blast-vsearch" / "blastdb.qza"
        helper   = Path(src_dir) / "qiime_helper.py"
        for in_fa, out_tsv in ((novel_fa, novel_out), (ctrl_fa, ctrl_out)):
            base = in_fa.with_suffix("")
            # Per-classifier path. blast, vsearch and naive-bayes all
            # imported the queries to the SAME "<base>_qiime.qza" and,
            # run concurrently, truncated each other mid-write. The
            # symptom is misleading: QIIME reports "missing the plugin(s)
            # necessary to load it" or "Archive does not contain a
            # correctly formatted VERSION file", both of which read like
            # a broken install rather than a race.
            qza  = Path(f"{base}_qiime.{args.classifier}.qza")
            cls_qza = out_tsv.with_suffix(".cls.qza")
            srh_qza = out_tsv.with_suffix(".search.qza")
            cmd_imp = (
                f"{CONDA_BIN} run -p {env} --no-capture-output python {helper} "
                f"import_fasta_to_qza {in_fa} {qza}"
            )
            _run(cmd_imp)
            cmd_run = (
                f"{CONDA_BIN} run -p {env} --no-capture-output bash -c \""
                f"qiime feature-classifier classify-consensus-blast "
                f"--i-query {qza} --i-blastdb {ref_db} --i-reference-taxonomy {ref_tax} "
                f"--p-num-threads {args.threads} --p-strand both --p-query-cov 0.80 "
                f"--p-evalue 0.001 --p-maxaccepts 10 --p-perc-identity 0.80 --p-min-consensus 0.51 "
                f"--o-classification {cls_qza} --o-search-results {srh_qza} && "
                f"python {helper} convert_qza_to_tsv {cls_qza} {out_tsv}\""
            )
            _run(cmd_run)

    elif args.classifier == "vsearch":
        ref_seqs = ref_root / "blast-vsearch" / "sequences.qza"
        ref_tax  = ref_root / "blast-vsearch" / "taxonomy.qza"
        helper   = Path(src_dir) / "qiime_helper.py"
        for in_fa, out_tsv in ((novel_fa, novel_out), (ctrl_fa, ctrl_out)):
            base = in_fa.with_suffix("")
            # Per-classifier path. blast, vsearch and naive-bayes all
            # imported the queries to the SAME "<base>_qiime.qza" and,
            # run concurrently, truncated each other mid-write. The
            # symptom is misleading: QIIME reports "missing the plugin(s)
            # necessary to load it" or "Archive does not contain a
            # correctly formatted VERSION file", both of which read like
            # a broken install rather than a race.
            qza  = Path(f"{base}_qiime.{args.classifier}.qza")
            cls_qza = out_tsv.with_suffix(".cls.qza")
            srh_qza = out_tsv.with_suffix(".search.qza")
            cmd_imp = (
                f"{CONDA_BIN} run -p {env} --no-capture-output python {helper} "
                f"import_fasta_to_qza {in_fa} {qza}"
            )
            _run(cmd_imp)
            maxrej = (f"--p-maxrejects {args.maxrejects} "
                      if args.maxrejects else "")
            cmd_run = (
                f"{CONDA_BIN} run -p {env} --no-capture-output bash -c \""
                f"qiime feature-classifier classify-consensus-vsearch "
                f"--i-query {qza} --i-reference-reads {ref_seqs} --i-reference-taxonomy {ref_tax} "
                f"--p-threads {args.threads} --p-strand both --p-query-cov 0.80 "
                f"--p-maxaccepts 10 --p-perc-identity 0.80 --p-min-consensus 0.51 {maxrej}"
                f"--o-classification {cls_qza} --o-search-results {srh_qza} && "
                f"python {helper} convert_qza_to_tsv {cls_qza} {out_tsv}\""
            )
            _run(cmd_run)

    elif args.classifier == "spingo":
        db = ref_root / "spingo" / "spingo_db.fasta"
        if spingo_env is None:
            sys.exit("missing 'spingo' conda env mapping in pipeline_config")
        for in_fa, out_tsv in ((novel_fa, novel_out), (ctrl_fa, ctrl_out)):
            raw = out_tsv.with_suffix(".raw.tsv")
            cmd = (
                f"bash -c \""
                f"{CONDA_BIN} run -p {spingo_env} --no-capture-output spingo "
                f"--database {db} --input {in_fa} --processors {args.threads} "
                f"--kmersize 8 --bootstrap 10 --subsample 8 > {raw} && "
                f"{CONDA_BIN} run -p {env} --no-capture-output python {src_dir}/process_spingo_output.py "
                f"{raw} {out_tsv}\""
            )
            _run(cmd)

    elif args.classifier == "kraken2":
        db = ref_root / "kraken2"
        for f in ("hash.k2d", "opts.k2d", "taxo.k2d"):
            if not (db / f).exists():
                sys.exit(f"missing {db/f}; run 02_train_classifier.py first")
        for in_fa, out_tsv in ((novel_fa, novel_out), (ctrl_fa, ctrl_out)):
            raw = out_tsv.with_suffix(".raw.tsv")
            report = out_tsv.with_suffix(".report.txt")
            # --report is required, not optional: without it
            # process_kraken2_output.py cannot reconstruct a lineage and fills
            # only the single rank kraken2's LCA landed on, leaving every
            # ancestor blank. The evaluator then scores those blanks as
            # abstentions, which understates kraken2 at every rank above its
            # LCA. Binaries resolved from the env prefix because `conda` is not
            # on PATH on compute nodes.
            cmd = (
                f"bash -c \""
                f"{env}/bin/kraken2 --db {db} --threads {args.threads} "
                f"--output {raw} --report {report} --use-names {in_fa} && "
                f"{env}/bin/python {src_dir}/process_kraken2_output.py "
                f"{raw} {out_tsv} 0.0 {report}\""
            )
            _run(cmd)

    elif args.classifier == "mapseq":
        ref_fa = ref_root / "mapseq" / "reference.fasta"
        ref_tax = ref_root / "mapseq" / "reference.tax"
        for f in (ref_fa, ref_tax):
            if not f.exists():
                sys.exit(f"missing {f}; run 02_train_classifier.py first")
        for in_fa, out_tsv in ((novel_fa, novel_out), (ctrl_fa, ctrl_out)):
            raw = out_tsv.with_suffix(".raw.tsv")
            # Same tophits/minscore/minid as the main panel, so a LOCO point is
            # comparable to the shipped one.
            cmd = (
                f"bash -c \""
                f"{env}/bin/mapseq -nthreads {args.threads} -tophits 20 "
                f"-minscore 30 -minid1 1 -minid2 1 "
                f"{in_fa} {ref_fa} {ref_tax} > {raw} && "
                f"{env}/bin/python {src_dir}/process_mapseq_output.py "
                f"{raw} {out_tsv}\""
            )
            _run(cmd)


if __name__ == "__main__":
    main()
