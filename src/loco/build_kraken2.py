#!/usr/bin/env python3
"""Build a Kraken2 database from a LOCO retained reference.

Wraps the same two steps the main panel used
(`training/kraken2/build_kraken2_database.sh`), with the same k-mer parameters,
so a LOCO database is comparable to the shipped one:

    1. training/kraken2/kraken2_data_converter.py  -> library.fna, nodes.dmp,
       names.dmp, seqid2taxid.map
    2. kraken2-build --add-to-library / --build

    python src/loco/build_kraken2.py <sequences.fasta> <taxonomy.tsv> <outdir> [threads]

Only hash.k2d, opts.k2d and taxo.k2d are kept; the intermediate library is
removed, because it is large and `data/ref_dbs/*/kraken2/` holds only those
three files.
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CONVERTER = ROOT / "training" / "kraken2" / "kraken2_data_converter.py"

# Same as build_kraken2_database.sh, so LOCO matches the main panel.
KMER_LEN = 31
MINIMIZER_LEN = 15
MINIMIZER_SPACES = 3


#: kraken2-build lives beside the interpreter running this script, and is not
#: necessarily on PATH on a compute node.
ENV_BIN = Path(sys.executable).parent
KRAKEN2_BUILD = ENV_BIN / "kraken2-build"


def run(cmd, **kw):
    print("[kraken2]", " ".join(str(c) for c in cmd), flush=True)
    env = kw.pop("env", None) or dict(os.environ)
    env["PATH"] = f"{ENV_BIN}:{env.get('PATH', '')}"
    subprocess.run([str(c) for c in cmd], check=True, env=env, **kw)


def main():
    if len(sys.argv) not in (4, 5):
        sys.exit(__doc__)
    seqs, tax, outdir = (Path(a) for a in sys.argv[1:4])
    threads = int(sys.argv[4]) if len(sys.argv) == 5 else 16

    outdir.mkdir(parents=True, exist_ok=True)
    staging = outdir / "_build"
    converted = staging / "converted"
    dbdir = staging / "db"
    for d in (converted, dbdir / "taxonomy"):
        d.mkdir(parents=True, exist_ok=True)

    run([sys.executable, CONVERTER, "--taxonomy", tax,
         "--sequences", seqs, "--output", converted])

    for required in ("library.fna", "nodes.dmp", "names.dmp"):
        if not (converted / required).exists():
            sys.exit(f"[kraken2] converter did not produce {required}")

    for dmp in ("nodes.dmp", "names.dmp"):
        shutil.copyfile(converted / dmp, dbdir / "taxonomy" / dmp)

    run([KRAKEN2_BUILD, "--add-to-library", converted / "library.fna",
         "--db", dbdir, "--threads", threads])
    run([KRAKEN2_BUILD, "--build", "--db", dbdir,
         "--kmer-len", KMER_LEN, "--minimizer-len", MINIMIZER_LEN,
         "--minimizer-spaces", MINIMIZER_SPACES, "--threads", threads])

    produced = []
    for k2d in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        src = dbdir / k2d
        if not src.exists():
            sys.exit(f"[kraken2] build did not produce {k2d}")
        shutil.copyfile(src, outdir / k2d)
        produced.append(k2d)

    # The staged library is several times the size of the finished index and
    # nothing downstream reads it. Inodes are the binding quota here.
    shutil.rmtree(staging, ignore_errors=True)
    print(f"[kraken2] wrote {', '.join(produced)} to {outdir}")


if __name__ == "__main__":
    main()
