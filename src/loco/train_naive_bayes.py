#!/usr/bin/env python3
"""Fit a QIIME2 naive-Bayes classifier on a held-out reference DB.

Usage: python train_naive_bayes.py <ref_fasta> <ref_taxonomy_tsv> <out_classifier.qza>
"""
import sys
from qiime2 import Artifact
from qiime2.plugins.feature_classifier.methods import fit_classifier_naive_bayes

if len(sys.argv) != 4:
    sys.stderr.write(
        "Usage: train_naive_bayes.py <ref_fasta> <ref_tax.tsv> <out.qza>\n"
    )
    sys.exit(1)

fa, tx, out = sys.argv[1:4]

print(f"[naive-bayes] importing sequences {fa}", flush=True)
seqs = Artifact.import_data("FeatureData[Sequence]", fa)

print(f"[naive-bayes] importing taxonomy {tx}", flush=True)
tax = Artifact.import_data(
    "FeatureData[Taxonomy]", tx, "HeaderlessTSVTaxonomyFormat"
)

print("[naive-bayes] fit_classifier_naive_bayes (this can take ~1 hour)...", flush=True)
res = fit_classifier_naive_bayes(reference_reads=seqs, reference_taxonomy=tax)
res.classifier.save(out)
print(f"[naive-bayes] saved {out}")
