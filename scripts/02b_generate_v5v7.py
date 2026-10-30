#!/usr/bin/env python3
"""V5-V7 mock generator that mirrors the existing per-replicate convention.

For each biome, we have 10 full-length env-seqs_<mock>XX.fasta files under
data/mock_community/<mock>/V1-V9/. We extract V5-V7 from each, run ART per
replicate (n_samples=5 like the existing layout), then run QIIME demux + DADA2
exactly like 02_generate.py.
"""
import argparse, glob, os, shutil, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
import _common

PRIMERS = {"f": "AACMGGATTAGATACCCKG", "r": "ACGTCATCCCCACCTTCC"}
LEN = (300, 500)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mock-name", required=True)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--n-samples", type=int, default=5,
                    help="ART replicates per env-seqs file (default 5, matches existing V4/V3-V4)")
    ap.add_argument("--count", type=int, default=30,
                    help="ART reads per input sequence (default 30, matches existing mocks)")
    ap.add_argument("--n-wildcard", action="store_true",
                    help="Truth builder: treat N in a source sequence as a wildcard, as in the published v1 truth (default: N is a mismatch)")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    import extract_reads   # needs the QIIME 2 environment (scikit-bio)
    import process_mock

    cfg = _common.load_pipeline_config()
    data_root = Path(cfg["data_dir"])
    mock_dir = data_root / "mock_community" / args.mock_name
    src_full_dir = mock_dir / "V1-V9"   # source of full-length env-seqs files
    region_dir = mock_dir / "V5-V7"
    region_dir.mkdir(parents=True, exist_ok=True)

    seq_fp = data_root / "mock_community" / "input_data" / f"{args.mock_name}_seqs.fasta"
    tax_fp = data_root / "mock_community" / "input_data" / f"{args.mock_name}_taxa.txt"

    # 1. Extract V5-V7 from each env-seqs_*XX.fasta (full-length)
    env_files = sorted(glob.glob(str(src_full_dir / f"env-seqs_{args.mock_name}*.fasta")))
    if not env_files:
        # fallback: just one round on the input fasta
        print(f"[warn] no per-rep env-seqs files in {src_full_dir}; falling back to {seq_fp}")
        env_files = [str(seq_fp)]

    extracted = []
    for ef in env_files:
        rep_id = Path(ef).stem.replace(f"env-seqs_{args.mock_name}", "")  # e.g. "00"
        out_fa = region_dir / f"env-seqs_{args.mock_name}{rep_id}.fasta"
        if out_fa.exists() and not args.force:
            print(f"  exists: {out_fa.name}")
        else:
            try:
                extract_reads.extract_amplicon(
                    input_fasta=ef, output_fasta=str(out_fa),
                    forward_primer=PRIMERS["f"], reverse_primer=PRIMERS["r"],
                    min_len=LEN[0], max_len=LEN[1], identity=0.8,
                    n_jobs=args.threads,
                )
            except Exception as exc:
                print(f"  ! extract failed for {ef}: {exc}")
                continue
        extracted.append((rep_id, str(out_fa)))

    # also assemble env_seqs.fasta (concat) for downstream tools that want it
    env_concat = region_dir / "env_seqs.fasta"
    with open(env_concat, "w") as out:
        for _, ef in extracted:
            with open(ef) as fh:
                shutil.copyfileobj(fh, out)
    print(f"[+] {len(extracted)} env-seqs V5-V7 files; concat -> {env_concat}")

    # also expose env_taxa.txt symlink (use the input one)
    env_tax_link = region_dir / "env_taxa.txt"
    if not env_tax_link.exists() and tax_fp.exists():
        os.symlink(str(tax_fp), env_tax_link)

    # 2. ART per replicate, with output_prefix = "<mock><rep_id>"
    art_dir = region_dir / "art"
    art_dir.mkdir(exist_ok=True)
    for rep_id, env_fa in extracted:
        prefix = f"{args.mock_name}{rep_id}"
        print(f"[ART] {prefix}")
        process_mock.art_simulation(
            input_fasta=env_fa,
            output_prefix=prefix,
            dir_out=str(art_dir),
            n_samples=args.n_samples,
            count=args.count,
            paired=True,
            parallel_jobs=args.threads,
            force=args.force,
        )

    # 3. QIIME demux + DADA2 (single output combining all art replicates)
    fq_dir = str(art_dir / "paired")
    try:
        demux, _, qiime_out = process_mock.qiime_demultiplex(
            out_dir=str(region_dir), fq_out_dir=fq_dir,
            paired=True, force=args.force)
        process_mock.denoise_dada2(
            qiime_out=qiime_out, mockname=args.mock_name,
            demux=demux, paired=True, threads=args.threads, force=args.force)
    except Exception as exc:
        print(f"  ! QIIME/DADA2 step failed: {exc}")
        raise

    # 4. Trueish taxonomy
    feature_table = region_dir / "paired" / f"table_{args.mock_name}.qza"
    rep_seqs_fa = region_dir / "paired" / f"rep_seqs_{args.mock_name}.fasta"
    true_tax_out = region_dir / f"trueish_taxonomy_{args.mock_name}.tsv"
    if feature_table.exists() and rep_seqs_fa.exists() and tax_fp.exists():
        try:
            process_mock.generate_trueish_taxonomies(
                data_dir=str(region_dir),
                feature_table=str(feature_table),
                expected_sequences=str(seq_fp),
                expected_taxonomy=str(tax_fp),
                observed_sequences=str(rep_seqs_fa),
                trueish_taxonomy=str(true_tax_out),
                force=args.force,
                n_wildcard=args.n_wildcard,
            )
        except Exception as exc:
            print(f"  ! trueish taxonomy failed: {exc}")
    else:
        print("  - skipped trueish taxonomy (missing inputs)")

    print(f"[done] {args.mock_name} V5-V7: rep_seqs has {sum(1 for l in open(rep_seqs_fa) if l.startswith('>'))} ASVs")


if __name__ == "__main__":
    main()
