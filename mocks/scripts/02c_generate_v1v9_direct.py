#!/usr/bin/env python3
"""Step 02c: V1-V9 generator for *curated* mocks (no MGnify, no backbone).

Mirrors the long-read recipe of 02d_generate_longread_biomes.py, except that
the source pool is the curated mock species directly
(input_data/<mock>_seqs.fasta), not the Greengenes2 backbone with MGnify
abundances.

Why a separate script: the biome generator's cell list and MGnify abundance
lookup are biome-specific. For curated mocks such as MBARC-26 the source
pool is simply the species FASTA that 02_generate.py would give to ART.

Recipe:
  1. Read data/mock_community/input_data/<mock>_seqs.fasta: already
     full-length 16S, one record per species (or per rrn copy for Zymo).
  2. For each of 10 replicates, draw a random subsample of `--seqs-per-rep`
     and apply per-base substitution at rate uniform [0.002, 0.02].
  3. Dereplicate union by MD5(seq) → write
     <mock>/V1-V9/source/single/rep_seqs_<mock>.fasta with DADA2-style
     hashed IDs.
  4. Write trueish-taxonomies.tsv mapping each MD5 ID to the source
     taxonomy from input_data/<mock>_taxa.txt.
  5. Build a matching .qza by importing the rep_seqs fasta into qiime2
     (FeatureData[Sequence]).

Outputs (mirror env-biome layout):
  data/mock_community/<mock>/V1-V9/env-seqs_<mock>{00..09}.fasta
  data/mock_community/<mock>/V1-V9/source/single/rep_seqs_<mock>.fasta
  data/mock_community/<mock>/V1-V9/source/single/rep_seqs_<mock>.qza
  data/mock_community/<mock>/V1-V9/trueish-taxonomies.tsv
"""
from __future__ import annotations

import argparse
import hashlib
import random
import shutil
import subprocess
import zlib
from pathlib import Path

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _common  # noqa: E402

_CFG = _common.load_pipeline_config()
INPUT_DIR = Path(_CFG["mock_input_dir"])
MOCK_DIR = Path(_CFG["data_dir"]) / "mock_community"

N_REPLICATES = 10
DEFAULT_SEQS_PER_REP = 75
# Match mockrobiota-16S: env-seqs are unmutated copies of source (verified
# 2026-05-10: mockrobiota env-seqs are 100% identical to input fasta after
# primer extraction). Only ART simulation downstream introduces noise.
# Keeping rates as configurable knobs for future use.
MUT_RATE_MIN = 0.0
MUT_RATE_MAX = 0.0
RNG_SEED = 20260510


def load_fasta(path: Path) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    cur, buf = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if cur is not None:
                    out.append((cur, "".join(buf).upper()))
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if cur is not None:
            out.append((cur, "".join(buf).upper()))
    return out


def load_taxa(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with open(path) as f:
        next(f)  # header "Feature ID\tTaxon"
        for line in f:
            line = line.rstrip("\n")
            if "\t" in line:
                fid, taxon = line.split("\t", 1)
                out[fid] = taxon
    return out


def mutate(seq: str, rate: float, rng: random.Random) -> str:
    if rate <= 0:
        return seq
    bases = "ACGT"
    n = len(seq)
    n_mut = round(rate * n)
    if n_mut <= 0:
        return seq
    positions = rng.sample(range(n), min(n_mut, n))
    s = list(seq)
    for pos in positions:
        b = s[pos]
        if b not in bases:
            continue
        s[pos] = rng.choice([x for x in bases if x != b])
    return "".join(s)


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        for h, s in records:
            f.write(f">{h}\n{s}\n")


def md5(seq: str) -> str:
    return hashlib.md5(seq.upper().encode()).hexdigest()


def import_to_qza(fasta: Path, qza: Path) -> None:
    """Import a fasta as FeatureData[Sequence] qza via qiime2."""
    cmd = ["qiime", "tools", "import",
           "--input-path", str(fasta),
           "--output-path", str(qza),
           "--type", "FeatureData[Sequence]"]
    print(f"  + qiime tools import → {qza.name}", flush=True)
    subprocess.run(cmd, check=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mock", required=True,
                    help="mock name; e.g. mbarc26 "
                         "(must have input_data/<mock>_seqs.fasta)")
    ap.add_argument("--seqs-per-rep", type=int, default=DEFAULT_SEQS_PER_REP)
    ap.add_argument("--n-replicates", type=int, default=N_REPLICATES)
    ap.add_argument("--no-qza", action="store_true",
                    help="skip qza import (useful for dry-run)")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    seqs_fa = INPUT_DIR / f"{args.mock}_seqs.fasta"
    taxa_tsv = INPUT_DIR / f"{args.mock}_taxa.txt"
    if not seqs_fa.exists() or not taxa_tsv.exists():
        ap.error(f"missing input: {seqs_fa} or {taxa_tsv}")

    out_root = MOCK_DIR / args.mock / "V1-V9"
    out_root.mkdir(parents=True, exist_ok=True)
    rep_dir = out_root / "source" / "single"
    rep_dir.mkdir(parents=True, exist_ok=True)

    print(f"=== {args.mock} V1-V9 (direct mutate) ===")
    pool = load_fasta(seqs_fa)
    taxa = load_taxa(taxa_tsv)
    print(f"  pool: {len(pool)} sequences from {seqs_fa.name}")
    print(f"  taxa: {len(taxa)} entries from {taxa_tsv.name}")
    if not pool:
        raise RuntimeError(f"empty input fasta: {seqs_fa}")

    # zlib.crc32 gives a stable cross-process seed; Python's built-in hash() is
    # process-local randomized and would produce different rep_seqs every run,
    # silently invalidating downstream classification results.
    rng = random.Random(RNG_SEED ^ zlib.crc32(args.mock.encode()))
    seqs_per_rep = min(args.seqs_per_rep, len(pool))

    all_records: list[tuple[str, str]] = []
    for i in range(args.n_replicates):
        # sample WITH replacement: pool may be smaller than seqs_per_rep
        # for small mocks (e.g. Zymo has 46, MBARC has 26).
        if seqs_per_rep <= len(pool):
            subsample = rng.sample(pool, seqs_per_rep)
        else:
            subsample = [rng.choice(pool) for _ in range(seqs_per_rep)]
        rate = rng.uniform(MUT_RATE_MIN, MUT_RATE_MAX)
        records = []
        for fid, seq in subsample:
            new_seq = mutate(seq, rate, rng)
            records.append((fid, new_seq))
            all_records.append((fid, new_seq))
        out_path = out_root / f"env-seqs_{args.mock}{i:02d}.fasta"
        write_fasta(out_path, records)
        avg_len = sum(len(s) for _, s in records) / len(records)
        print(f"  + replicate {i:02d}: rate={rate:.4f}, n={len(records)}, "
              f"avg_len={avg_len:.0f} bp -> {out_path.name}")

    # Dereplicate union by MD5
    seen: dict[str, str] = {}
    md5_to_src: dict[str, str] = {}
    for src_fid, seq in all_records:
        h = md5(seq)
        if h not in seen:
            seen[h] = seq
            md5_to_src[h] = src_fid
    rep_records = [(h, s) for h, s in seen.items()]
    rep_fa = rep_dir / f"rep_seqs_{args.mock}.fasta"
    write_fasta(rep_fa, rep_records)
    avg_rep = sum(len(s) for _, s in rep_records) / len(rep_records)
    print(f"  + rep_seqs: {len(rep_records)} unique ASVs "
          f"(avg_len={avg_rep:.0f} bp) -> {rep_fa}")

    # trueish-taxonomies.tsv: MD5 -> source taxonomy (truncated to species)
    tax_out = out_root / "trueish-taxonomies.tsv"
    n_with_tax = 0
    with open(tax_out, "w") as f:
        f.write("id\ttaxonomy\n")
        for h in seen:
            src = md5_to_src[h]
            t = taxa.get(src, "")
            if t:
                n_with_tax += 1
            f.write(f"{h}\t{t}\n")
    print(f"  + trueish-taxonomies: {n_with_tax}/{len(seen)} mapped -> "
          f"{tax_out}")

    if not args.no_qza:
        rep_qza = rep_dir / f"rep_seqs_{args.mock}.qza"
        if rep_qza.exists() and not args.force:
            print(f"  - qza already exists: {rep_qza}")
        else:
            if rep_qza.exists():
                rep_qza.unlink()
            import_to_qza(rep_fa, rep_qza)

    print("\ndone.")


if __name__ == "__main__":
    main()
