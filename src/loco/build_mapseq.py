#!/usr/bin/env python3
"""Build a MAPseq reference from a LOCO retained reference.

MAPseq wants two files: a plain FASTA, and a taxonomy file whose format is
copied here from the existing `data/ref_dbs/*/mapseq/reference.tax`:

    #cutoff: 0.00:0.08 0.70:0.35 0.70:0.35 0.70:0.35 0.80:0.25 0.92:0.08 0.95:0.05
    #name: <name>
    #levels: Domain Phylum Class Order Family Genus Species
    <id>\tDomain;Phylum;Class;Order;Family;Genus;Species

Rank prefixes (`d__`, `p__`, ...) are stripped; absent ranks are left empty, as
in the shipped reference. MAPseq builds its own `.mscluster` index on first use,
so there is no separate index step.

    python src/loco/build_mapseq.py <sequences.fasta> <taxonomy.tsv> <outdir>
"""
import shutil
import sys
from pathlib import Path

# Copied verbatim from the shipped reference so LOCO and the main panel score
# identically. Per-rank "identity:width" pairs MAPseq uses for its confidence.
CUTOFF_LINE = ("#cutoff: 0.00:0.08 0.70:0.35 0.70:0.35 0.70:0.35 "
               "0.80:0.25 0.92:0.08 0.95:0.05")
LEVELS_LINE = "#levels: Domain Phylum Class Order Family Genus Species"


def strip_prefix(cell):
    cell = cell.strip()
    return cell.split("__", 1)[-1].strip() if "__" in cell else cell


def write_tax(taxonomy_tsv, out_tax, name):
    n = 0
    with open(taxonomy_tsv) as fh, open(out_tax, "w") as out:
        out.write(f"{CUTOFF_LINE}\n#name: {name}\n{LEVELS_LINE}\n")
        for i, line in enumerate(fh):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            if i == 0 and parts[0] in ("Feature ID", "id", "#OTU ID"):
                continue
            ranks = [strip_prefix(c) for c in parts[1].split(";")]
            # Pad or trim to exactly seven ranks; MAPseq indexes them by position.
            ranks = (ranks + [""] * 7)[:7]
            out.write(f"{parts[0]}\t{';'.join(ranks)}\n")
            n += 1
    return n


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    seqs, tax, outdir = (Path(a) for a in sys.argv[1:])
    outdir.mkdir(parents=True, exist_ok=True)

    ref_fa = outdir / "reference.fasta"
    ref_tax = outdir / "reference.tax"

    # MAPseq reads the FASTA directly, so a copy is enough. Copied rather than
    # symlinked because MAPseq writes its .mscluster index alongside it.
    shutil.copyfile(seqs, ref_fa)
    n_seq = sum(1 for line in open(ref_fa) if line.startswith(">"))
    n_tax = write_tax(tax, ref_tax, f"LOCO_{outdir.parents[1].name}")

    print(f"[mapseq] {ref_fa}: {n_seq} sequences")
    print(f"[mapseq] {ref_tax}: {n_tax} taxonomy rows")
    if n_seq != n_tax:
        # Not fatal for MAPseq, but it means the reference and the taxonomy
        # disagree, which would silently drop calls. Say so loudly.
        print(f"[mapseq] WARNING: {n_seq} sequences but {n_tax} taxonomy rows",
              file=sys.stderr)


if __name__ == "__main__":
    main()
