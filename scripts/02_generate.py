#!/usr/bin/env python3
"""Step 02 - Generate per-region mock community.

For each (mock x region):
    - extract V-region sequences via primers (V4 / V3-V4 / V5-V7 / V1-V9), or
      use full-length GROND sequences directly (RRN)
    - simulate Illumina reads with ART
    - run DADA2-style processing (qiime_demultiplex + denoise_dada2)
    - generate trueish taxonomy

Outputs go under {data_dir}/mock_community/{mock_name}/{region}/.
"""

import argparse
import os
import shutil
from os.path import join, exists
from pathlib import Path

import _common


REGIONS_ALL = ["V4", "V3-V4", "V5-V7", "V1-V9", "RRN"]


def main():
    ap = argparse.ArgumentParser(description="Generate mock-community per-region datasets")
    ap.add_argument("--mock-name", required=True,
                    help="Mock name (mock-3, mock-12, ...) or biome name")
    ap.add_argument("--regions", nargs="+", default=REGIONS_ALL,
                    choices=REGIONS_ALL)
    ap.add_argument("--out-dir", default=None,
                    help="Output dir (default: data/mock_community/{mock}/)")
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--n-samples", type=int, default=1)
    ap.add_argument("--count", type=int, default=1,
                    help="ART reads per input sequence (default 1; set ~30 for "
                         "biome mocks to match 02b_generate_v5v7 density)")
    ap.add_argument("--length", type=int, default=200,
                    help="ART read length (default 200). For V3-V4 paired use ~250 so "
                         "the 460 bp amplicon can be merged.")
    ap.add_argument("--paired", action="store_true",
                    help="Simulate paired-end reads (default single-end). "
                         "Don't use for V1-V9: 200+200bp can't merge ~1500bp.")
    ap.add_argument("--n-wildcard", action="store_true",
                    help="Truth builder: treat N in a source sequence as a wildcard, as in the published v1 truth (default: N is a mismatch)")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    # Guard: V1-V9 and RRN are long-read regions and must NOT go through
    # ART/DADA2. ART caps at ~300 bp Illumina reads, so applying it to
    # V1-V9 (~1500 bp) or RRN (~5000 bp) silently produces 200 bp rep_seqs
    # the V1-V9/RRN models cannot classify. Use the long-read generators
    # instead: 02d_generate_longread_biomes.py (biome mocks) or
    # 02c_generate_v1v9_direct.py (curated mocks).
    long_read_regions = {"V1-V9", "RRN"}
    bad = sorted(set(args.regions) & long_read_regions)
    if bad:
        ap.error(
            f"02_generate.py must not be used for long-read regions {bad}. "
            "ART produces ≤200 bp reads and the V1-V9/RRN models require full-length "
            "inputs. Use scripts/02d_generate_longread_biomes.py (biome mocks) "
            "or scripts/02c_generate_v1v9_direct.py (curated mocks), which mutate "
            "full-length sequences directly."
        )

    cfg = _common.load_pipeline_config()
    _common.add_src_to_path()
    threads = args.threads or cfg.get("default_n_threads", 4)

    project_data = Path(cfg["data_dir"])
    input_dir = Path(cfg["mock_input_dir"])
    out_root = Path(args.out_dir) if args.out_dir \
        else project_data / "mock_community" / args.mock_name
    out_root.mkdir(parents=True, exist_ok=True)

    # Locate input expected sequences + taxonomy
    seq_fp = input_dir / f"{args.mock_name}_seqs.fasta"
    tax_fp = input_dir / f"{args.mock_name}_taxa.txt"
    if not seq_fp.exists():
        raise FileNotFoundError(f"Input fasta not found: {seq_fp}")

    print(f"[02_generate] mock = {args.mock_name}")
    print(f"[02_generate] regions = {args.regions}")
    print(f"[02_generate] out = {out_root}")

    if args.dry_run:
        for region in args.regions:
            print(f"  [dry] would generate region={region}")
        return

    # Heavy imports only when actually running
    import utils                 # noqa
    import process_mock          # noqa
    import extract_reads         # noqa
    from Bio import SeqIO        # noqa

    regions_cfg = _common.load_regions_config()

    for region in args.regions:
        region_dir = out_root / region
        region_dir.mkdir(parents=True, exist_ok=True)
        print(f"\n[02_generate] === region {region} ===")

        # 1. Extract V-region (or use full-length for RRN)
        region_fasta = region_dir / f"{args.mock_name}_{region}.fasta"
        if region_fasta.exists() and not args.force:
            print(f"  - region fasta already exists: {region_fasta}")
        else:
            if regions_cfg[region].get("full_length"):
                shutil.copy(seq_fp, region_fasta)
                print(f"  + copied full-length seqs (RRN) -> {region_fasta.name}")
            else:
                primers = regions_cfg[region]["primers"]
                lengths = regions_cfg[region]["length"]
                # Use extract_reads.extract_amplicon if available; otherwise
                # fall back to a simple primer-trim heuristic.
                if hasattr(extract_reads, "extract_amplicon"):
                    extract_reads.extract_amplicon(
                        input_fasta=str(seq_fp),
                        output_fasta=str(region_fasta),
                        forward_primer=primers["f"],
                        reverse_primer=primers["r"],
                        min_len=lengths["min"],
                        max_len=lengths["max"],
                    )
                else:
                    # Naive: pass through; downstream tools may handle it
                    shutil.copy(seq_fp, region_fasta)
                print(f"  + extracted {region} -> {region_fasta.name}")

        # 2. ART simulation
        art_dir = region_dir / "art"
        art_dir.mkdir(exist_ok=True)
        print(f"  + ART simulation -> {art_dir}")
        process_mock.art_simulation(
            input_fasta=str(region_fasta),
            output_prefix=args.mock_name,
            dir_out=str(art_dir),
            n_samples=args.n_samples,
            count=args.count,
            length=args.length,
            paired=args.paired,
            parallel_jobs=threads,
            force=args.force,
        )

        # 3. QIIME demux + DADA2 denoise
        try:
            fq_dir = str(art_dir / ("paired" if args.paired else "single"))
            demux, _, qiime_out = process_mock.qiime_demultiplex(
                out_dir=str(region_dir), fq_out_dir=fq_dir,
                paired=args.paired, force=args.force)
            process_mock.denoise_dada2(
                qiime_out=qiime_out, mockname=args.mock_name,
                demux=demux, paired=args.paired,
                threads=threads, force=args.force)
        except Exception as exc:
            print(f"  ! QIIME/DADA2 step failed: {exc}")

        # 4. Trueish taxonomy
        feature_table = join(region_dir, ("paired" if args.paired else "single"),
                             f"table_{args.mock_name}.qza")
        rep_seqs_fa   = join(region_dir, ("paired" if args.paired else "single"),
                             f"rep_seqs_{args.mock_name}.fasta")
        true_tax_out = region_dir / f"trueish_taxonomy_{args.mock_name}.tsv"
        if exists(feature_table) and exists(rep_seqs_fa) and tax_fp.exists():
            try:
                process_mock.generate_trueish_taxonomies(
                    data_dir=str(region_dir),
                    feature_table=feature_table,
                    expected_sequences=str(seq_fp),
                    expected_taxonomy=str(tax_fp),
                    observed_sequences=rep_seqs_fa,
                    trueish_taxonomy=str(true_tax_out),
                    force=args.force,
                    n_wildcard=args.n_wildcard,
                )
            except Exception as exc:
                print(f"  ! trueish taxonomy failed: {exc}")
        else:
            print("  - skipped trueish taxonomy (missing inputs)")

    print("\n[02_generate] done")


if __name__ == "__main__":
    main()
