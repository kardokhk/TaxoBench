#!/usr/bin/env python3
"""LOCO step 03 — classify novel + control queries with a single classifier.

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

VALID = ("taxonova", "idtaxa", "naive-bayes", "blast", "vsearch", "spingo")


def _run(cmd):
    print("[loco/03]", cmd, flush=True)
    subprocess.run(cmd, shell=True, check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--region", required=True,
                    choices=["V4", "V3-V4", "V5-V7", "V1-V9", "RRN"])
    ap.add_argument("--classifier", required=True, choices=VALID)
    ap.add_argument("--threads", type=int, default=20)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    cfg = _common.load_pipeline_config()
    src_dir = cfg["scripts_dir"]
    env = cfg["conda_envs"]["qiime2"]
    spingo_env = cfg["conda_envs"].get("spingo")
    if isinstance(spingo_env, dict):
        spingo_env = spingo_env.get("spingo_env")

    loco_root = Path(cfg["data_dir"]) / "loco" / args.region
    ref_root  = loco_root / "retained_ref"
    novel_fa  = loco_root / "novel_queries.fasta"
    ctrl_fa   = loco_root / "control_queries.fasta"
    if not novel_fa.exists():
        sys.exit(f"missing {novel_fa}; run 01_make_holdout.py first")

    out_root = Path(cfg["results_dir"]) / "loco" / args.region / args.classifier
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
                f"conda run -p {env} --no-capture-output bash -c \""
                f"taxonova classification "
                f"--input {in_fa} --model_input {model} --output {raw} "
                f"--confidence_threshold 0.0 && "
                f"python {src_dir}/process_taxonova_output.py "
                f"{raw} {out_tsv} False\""
            )
            _run(cmd)

    elif args.classifier == "idtaxa":
        model = ref_root / "idtaxa" / "trainingSet.RData"
        for in_fa, out_tsv in ((novel_fa, novel_out), (ctrl_fa, ctrl_out)):
            raw = out_tsv.with_suffix(".raw.tsv")
            cmd = (
                f"conda run -p {env} --no-capture-output bash -c \""
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
            qza  = Path(str(base) + "_qiime.qza")
            cls_qza = out_tsv.with_suffix(".qza")
            cmd_imp = (
                f"conda run -p {env} --no-capture-output python {helper} "
                f"import_fasta_to_qza {in_fa} {qza}"
            )
            _run(cmd_imp)
            cmd_run = (
                f"conda run -p {env} --no-capture-output bash -c \""
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
            qza  = Path(str(base) + "_qiime.qza")
            cls_qza = out_tsv.with_suffix(".cls.qza")
            srh_qza = out_tsv.with_suffix(".search.qza")
            cmd_imp = (
                f"conda run -p {env} --no-capture-output python {helper} "
                f"import_fasta_to_qza {in_fa} {qza}"
            )
            _run(cmd_imp)
            cmd_run = (
                f"conda run -p {env} --no-capture-output bash -c \""
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
            qza  = Path(str(base) + "_qiime.qza")
            cls_qza = out_tsv.with_suffix(".cls.qza")
            srh_qza = out_tsv.with_suffix(".search.qza")
            cmd_imp = (
                f"conda run -p {env} --no-capture-output python {helper} "
                f"import_fasta_to_qza {in_fa} {qza}"
            )
            _run(cmd_imp)
            cmd_run = (
                f"conda run -p {env} --no-capture-output bash -c \""
                f"qiime feature-classifier classify-consensus-vsearch "
                f"--i-query {qza} --i-reference-reads {ref_seqs} --i-reference-taxonomy {ref_tax} "
                f"--p-threads {args.threads} --p-strand both --p-query-cov 0.80 "
                f"--p-maxaccepts 10 --p-perc-identity 0.80 --p-min-consensus 0.51 "
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
                f"conda run -p {spingo_env} --no-capture-output spingo "
                f"--database {db} --input {in_fa} --processors {args.threads} "
                f"--kmersize 8 --bootstrap 10 --subsample 8 > {raw} && "
                f"conda run -p {env} --no-capture-output python {src_dir}/process_spingo_output.py "
                f"{raw} {out_tsv}\""
            )
            _run(cmd)


if __name__ == "__main__":
    main()
