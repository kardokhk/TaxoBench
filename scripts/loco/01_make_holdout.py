#!/usr/bin/env python3
"""LOCO step 01: generate stratified holdout for a region.

Wraps src/loco/make_holdout.py with the standard taxobench paths.
"""
import argparse
import subprocess
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import _common  # noqa


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--region", required=True,
                    choices=["V4", "V3-V4", "V5-V7", "V1-V9", "RRN"])
    ap.add_argument("--seed", type=int, default=2026)
    ap.add_argument("--n-species",  type=int, default=200)
    ap.add_argument("--n-genera",   type=int, default=100)
    ap.add_argument("--n-families", type=int, default=50)
    ap.add_argument("--n-orders",   type=int, default=25)
    ap.add_argument("--n-classes",  type=int, default=10)
    ap.add_argument("--n-controls", type=int, default=2000)
    ap.add_argument("--fixed-higher-ranks", default=None, metavar="MANIFEST",
                    help="Reuse the genus, family, order and class hold-outs from a "
                         "previous replicate's holdout_manifest.json and redraw only "
                         "the species stratum, so replicates stay comparable in "
                         "reference size.")
    ap.add_argument("--loco-root", default="loco",
                    help="Subdirectory of data/ to write into. The default "
                         "'loco' is the shipped 200-species holdout that every "
                         "published MC/OC number depends on; use a different "
                         "root (e.g. 'loco2k') for the enlarged 2,000-species replicates "
                         "so they cannot overwrite it.")
    ap.add_argument("--tag", default=None,
                    help="Optional suffix on the region directory, e.g. "
                         "'seed7', so replicates of the same region coexist.")
    args = ap.parse_args()

    cfg = _common.load_pipeline_config()
    ref_dir  = Path(cfg["reference_database_dir"]) / args.region
    region_dir = args.region if args.tag is None else f"{args.region}__{args.tag}"
    out_dir  = Path(cfg["data_dir"]) / args.loco_root / region_dir
    if out_dir.exists() and any(out_dir.iterdir()):
        # Overwriting a holdout silently invalidates every result derived from
        # it, and the derived results do not change their timestamps.
        raise SystemExit(
            f"[loco/01] refusing to overwrite non-empty {out_dir}. "
            "Pass --loco-root/--tag to write somewhere new, or remove it "
            "deliberately.")

    cmd = [
        sys.executable, str(ROOT / "src" / "loco" / "make_holdout.py"),
        "--region", args.region,
        "--ref-dir", str(ref_dir),
        "--out-dir", str(out_dir),
        "--seed", str(args.seed),
        "--n-species", str(args.n_species),
        "--n-genera", str(args.n_genera),
        "--n-families", str(args.n_families),
        "--n-orders", str(args.n_orders),
        "--n-classes", str(args.n_classes),
        "--n-controls", str(args.n_controls),
    ]
    if args.fixed_higher_ranks:
        cmd += ["--fixed-higher-ranks", str(args.fixed_higher_ranks)]
    print("[loco/01]", " ".join(cmd))
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
