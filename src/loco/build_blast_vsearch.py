#!/usr/bin/env python3
"""Build BLAST/VSEARCH QZA artifacts for the held-out LOCO reference DB.

Produces under <out_dir>:
    sequences.qza, taxonomy.qza   (used by both BLAST + VSEARCH)
    blastdb.qza                   (BLAST only)

Usage: python build_blast_vsearch.py <ref_fasta> <ref_tax_tsv> <out_dir>
"""
import sys
from pathlib import Path
from qiime2 import Artifact
from qiime2.plugins.feature_classifier.methods import makeblastdb

if len(sys.argv) != 4:
    sys.stderr.write("Usage: build_blast_vsearch.py <ref.fasta> <ref.tsv> <out_dir>\n")
    sys.exit(1)

fa, tx, out_dir = sys.argv[1:4]
out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)

print(f"[bv] importing sequences from {fa}")
seqs = Artifact.import_data("FeatureData[Sequence]", fa)
seqs.save(str(out / "sequences.qza"))

print(f"[bv] importing taxonomy from {tx}")
tax = Artifact.import_data(
    "FeatureData[Taxonomy]", tx, "HeaderlessTSVTaxonomyFormat"
)
tax.save(str(out / "taxonomy.qza"))

print("[bv] building BLAST DB")
db = makeblastdb(seqs)
db.database.save(str(out / "blastdb.qza"))
print("[bv] done")
