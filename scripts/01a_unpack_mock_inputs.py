#!/usr/bin/env python3
"""Step 01a: place the published mock inputs in the input directory.

mocks/inputs/ holds, gzipped, the <mock>_seqs.fasta and <mock>_taxa.txt files
behind the published cells: Mockrobiota-16S, MBARC-26, the gut and vaginal
mocks, the nine biome mocks and their RRN counterparts. This step writes them,
uncompressed, to mock_input_dir (data/mock_community/input_data by default) and
checks each one against mocks/inputs/MD5SUMS.

Existing files are left alone unless --force is given.
"""
import argparse
import gzip
import hashlib
from pathlib import Path

import _common


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--mock-name", nargs="*", default=None,
                    help="only these mocks (for example mbarc26 biome_soil); default all")
    ap.add_argument("--out-dir", default=None,
                    help="default: mock_input_dir from pipeline_config.json")
    ap.add_argument("--force", action="store_true", help="overwrite existing files")
    ap.add_argument("--dry-run", action="store_true", help="list what would be written")
    args = ap.parse_args()

    src = _common.project_root() / "mocks" / "inputs"
    out = Path(args.out_dir or _common.load_pipeline_config()["mock_input_dir"])
    sums = {}
    for line in (src / "MD5SUMS").read_text().splitlines():
        digest, name = line.split()
        sums[name] = digest

    names = sorted(sums)
    if args.mock_name:
        wanted = set(args.mock_name)
        names = [n for n in names if n.rsplit("_", 1)[0] in wanted]
        missing = wanted - {n.rsplit("_", 1)[0] for n in names}
        if missing:
            raise SystemExit(f"[01a] not in mocks/inputs: {', '.join(sorted(missing))}")

    if not args.dry_run:
        out.mkdir(parents=True, exist_ok=True)
    n_written = 0
    for name in names:
        dest = out / name
        if dest.exists() and not args.force:
            print(f"[01a] keep {dest}")
            continue
        if args.dry_run:
            print(f"[01a] would write {dest}")
            continue
        data = gzip.decompress((src / f"{name}.gz").read_bytes())
        if hashlib.md5(data).hexdigest() != sums[name]:
            raise SystemExit(f"[01a] md5 mismatch for {name}")
        dest.write_bytes(data)
        n_written += 1
    if not args.dry_run:
        print(f"[01a] {n_written} files written to {out}")


if __name__ == "__main__":
    main()
