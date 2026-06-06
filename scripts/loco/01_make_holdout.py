#!/usr/bin/env python3
"""LOCO step 01 — generate stratified holdout for a region.

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
    args = ap.parse_args()

    cfg = _common.load_pipeline_config()
    ref_dir  = Path(cfg["reference_database_dir"]) / args.region
    out_dir  = Path(cfg["data_dir"]) / "loco" / args.region

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
    print("[loco/01]", " ".join(cmd))
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
