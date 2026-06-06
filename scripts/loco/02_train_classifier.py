#!/usr/bin/env python3
"""LOCO step 02 — train ONE classifier on the held-out reference DB.

Designed for SLURM: one job = one (region, classifier). This keeps the GPU
training (taxonova) separable from CPU-bound work.

Supported classifiers:
  taxonova       — GPU; smart_full_data training
  idtaxa         — CPU; R DECIPHER::LearnTaxa
  naive-bayes    — CPU; QIIME2 fit-classifier-naive-bayes
  blast-vsearch  — CPU; build QZA artifacts + BLAST DB (used by both BLAST and VSEARCH)
  spingo         — CPU; reformat fasta to SPINGO header layout (no separate training)
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import _common  # noqa

VALID = ("taxonova", "idtaxa", "naive-bayes", "blast-vsearch", "spingo")


def _spingo_format(merged_fasta, taxonomy_tsv, out_fasta):
    """Convert merged.fasta (>id\\ttax_string) to SPINGO format.

    SPINGO expects ``>id\\tspecies_label`` where species_label has the form
    "Genus species" (no rank prefixes). When species is missing we fall back
    to the deepest-available rank, prepended by 'g__' / 'f__' so the parser
    doesn't crash. SPINGO ignores levels above species but we keep the full
    lineage as a 3rd tab-separated field for traceability.
    """
    tax_map = {}
    with open(taxonomy_tsv) as fh:
        for i, line in enumerate(fh):
            parts = line.rstrip("\n").split("\t")
            if i == 0 and parts[0] in ("Feature ID", "id"):
                continue
            if len(parts) >= 2:
                tax_map[parts[0]] = parts[1]
    n = 0
    with open(out_fasta, "w") as out, open(merged_fasta) as fh:
        cur_id, seq_lines = None, []

        def emit():
            nonlocal n
            if cur_id is None:
                return
            tax = tax_map.get(cur_id, "")
            tparts = [p.strip() for p in tax.split(";")]
            sp_label = ""
            genus = ""
            for p in tparts:
                if p.startswith("g__") and p[3:]:
                    genus = p[3:]
                if p.startswith("s__") and p[3:]:
                    sp_label = p[3:]
                    break
            if not sp_label and genus:
                sp_label = genus + " sp"
            if not sp_label:
                sp_label = (tparts[-1].split("__", 1)[-1] or "Unclassified")
            if not genus:
                genus = sp_label.split()[0] if sp_label else "Unclassified"
            # SPINGO format: ">id<TAB>genus species<TAB>genus"
            out.write(f">{cur_id}\t{sp_label}\t{genus}\n")
            out.write("".join(seq_lines) + "\n")
            n += 1

        for line in fh:
            if line.startswith(">"):
                emit()
                cur_id = line[1:].rstrip("\n").split("\t", 1)[0].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line.rstrip("\n"))
        emit()
    print(f"[spingo] wrote {n} sequences -> {out_fasta}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--region", required=True,
                    choices=["V4", "V3-V4", "V5-V7", "V1-V9", "RRN"])
    ap.add_argument("--classifier", required=True, choices=VALID)
    ap.add_argument("--threads", type=int, default=20)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    cfg = _common.load_pipeline_config()
    base = Path(cfg["data_dir"]) / "loco" / args.region / "retained_ref"
    fa  = base / "sequences.fasta"
    tax = base / "taxonomy.tsv"
    merged = base / "merged.fasta"
    if not fa.exists():
        sys.exit(f"missing {fa}; run 01_make_holdout.py first")

    cls_dir = base / args.classifier
    cls_dir.mkdir(parents=True, exist_ok=True)

    if args.classifier == "taxonova":
        out_model = cls_dir / "model.pt"
        if out_model.exists() and not args.force:
            print(f"[loco/02] skip taxonova: {out_model} exists")
            return
        env = cfg["conda_envs"]["qiime2"]
        # Canonical tx_wip_09 per-region recipe. The previous --smart_full_data
        # invocation underfit at ~15 epochs / loss 3.4 / overall F1 0.14, which
        # collapsed species/genus/family softmax mass to ~0 across all controls.
        # See corrections/comments_addressed/diagnosis_sweep_loco_2026-05-20.md.
        recipe = {
            "V4":    {"epochs": 45, "tau": None},
            "V3-V4": {"epochs": 20, "tau": 0.6},
            "V5-V7": {"epochs": 90, "tau": 0.5},
            "V1-V9": {"epochs": 60, "tau": 0.5},
            "RRN":   {"epochs": 60, "tau": None},
        }[args.region]
        raw_model = cls_dir / "model.raw.pt"
        cmd = [
            "conda", "run", "-p", env, "--no-capture-output",
            "taxonova", "training", "--verbose",
            "--train_db", str(merged),
            "--model_output", str(raw_model),
            "--epochs", str(recipe["epochs"]),
            "--batch_size", "256",
            "--learning_rate", "1e-3",
            "--kmer-length", "8",
            "--val_split", "0",
            "--early_stopping",
        ]
        print("[loco/02 taxonova:train]", " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True)
        if recipe["tau"] is not None:
            cmd_tau = [
                "conda", "run", "-p", env, "--no-capture-output",
                sys.executable, "-m", "taxonova.calibration.tau_norm",
                "--in_ckpt", str(raw_model),
                "--out_ckpt", str(out_model),
                "--tau", str(recipe["tau"]),
            ]
            print("[loco/02 taxonova:tau-norm]", " ".join(cmd_tau), flush=True)
            subprocess.run(cmd_tau, check=True)
        else:
            os.replace(raw_model, out_model)

    elif args.classifier == "idtaxa":
        out_model = cls_dir / "trainingSet.RData"
        if out_model.exists() and not args.force:
            print(f"[loco/02] skip idtaxa: {out_model} exists")
            return
        env = cfg["conda_envs"]["qiime2"]
        cmd = [
            "conda", "run", "-p", env, "--no-capture-output",
            "Rscript", str(ROOT / "src" / "loco" / "train_idtaxa.R"),
            str(fa), str(tax), str(out_model),
        ]
        print("[loco/02 idtaxa]", " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True)

    elif args.classifier == "naive-bayes":
        out_model = cls_dir / "classifier.qza"
        if out_model.exists() and not args.force:
            print(f"[loco/02] skip naive-bayes: {out_model} exists")
            return
        env = cfg["conda_envs"]["qiime2"]
        cmd = [
            "conda", "run", "-p", env, "--no-capture-output",
            sys.executable,
            str(ROOT / "src" / "loco" / "train_naive_bayes.py"),
            str(fa), str(tax), str(out_model),
        ]
        print("[loco/02 naive-bayes]", " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True)

    elif args.classifier == "blast-vsearch":
        seqs_qza = cls_dir / "sequences.qza"
        tax_qza  = cls_dir / "taxonomy.qza"
        bdb_qza  = cls_dir / "blastdb.qza"
        if seqs_qza.exists() and tax_qza.exists() and bdb_qza.exists() and not args.force:
            print(f"[loco/02] skip blast-vsearch: artifacts exist in {cls_dir}")
            return
        env = cfg["conda_envs"]["qiime2"]
        cmd = [
            "conda", "run", "-p", env, "--no-capture-output",
            sys.executable,
            str(ROOT / "src" / "loco" / "build_blast_vsearch.py"),
            str(fa), str(tax), str(cls_dir),
        ]
        print("[loco/02 blast-vsearch]", " ".join(cmd), flush=True)
        subprocess.run(cmd, check=True)

    elif args.classifier == "spingo":
        out_fa = cls_dir / "spingo_db.fasta"
        if out_fa.exists() and not args.force:
            print(f"[loco/02] skip spingo: {out_fa} exists")
            return
        _spingo_format(str(merged), str(tax), str(out_fa))


if __name__ == "__main__":
    main()
