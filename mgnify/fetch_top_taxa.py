#!/usr/bin/env python3
"""Fetch the top-N taxa per rank for one or more MGnify biomes.

Walks: biome -> samples -> runs -> analyses -> taxonomy entries, then
aggregates abundances per (biome, rank) and writes one CSV per rank.

Usage
-----
  # by short-name (resolved against biomes.csv)
  python fetch_top_taxa.py --biomes marine soil --max-samples 50

  # by full MGnify lineage
  python fetch_top_taxa.py --lineage 'root:Environmental:Aquatic:Marine' \\
                           --label marine --max-samples 50

Output
------
  results/<run-tag>/<biome>__<rank>__top<N>.csv
      columns: rank, taxon, taxa_level, total_abundance
"""

import argparse
import sys
import time
import urllib.parse
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import requests


API = "https://www.ebi.ac.uk/metagenomics/api/v1"
RANKS = ["species", "genus", "family", "order", "class", "phylum", "kingdom", "domain"]
SKIP_NAMES = {"unclassified", "unknown", "unidentified", "uncultured",
              "environmental samples", "metagenome", ""}
TAXONOMY_ENDPOINTS = ("taxonomy/ssu", "taxonomy", "taxonomy/lsu", "taxonomy/unite")
HERE = Path(__file__).parent


# --------------------------------------------------------------------------- #
# HTTP helpers
# --------------------------------------------------------------------------- #

def make_session() -> requests.Session:
    s = requests.Session()
    s.headers.update({
        "Accept": "application/json",
        "User-Agent": "taxobench-mgnify/1.0",
    })
    return s


def api_get(session: requests.Session, url: str, params: dict = None,
            retries: int = 3, verbose: bool = False) -> Optional[dict]:
    """GET with exponential back-off on 429. Returns parsed JSON or None."""
    for attempt in range(retries):
        try:
            r = session.get(url, params=params, timeout=30)
            if r.status_code == 200:
                return r.json()
            if r.status_code == 404:
                return None
            if r.status_code == 429:
                wait = 2 ** attempt
                if verbose:
                    print(f"    429 rate-limited; sleeping {wait}s", file=sys.stderr)
                time.sleep(wait)
                continue
            if verbose:
                print(f"    HTTP {r.status_code} for {url}", file=sys.stderr)
            return None
        except requests.RequestException as e:
            if verbose:
                print(f"    request error ({e}); attempt {attempt + 1}", file=sys.stderr)
            time.sleep(1)
    return None


def paginate(session: requests.Session, url: str, params: dict = None,
             limit: int = None, verbose: bool = False):
    """Iterate over MGnify paginated `data` arrays."""
    params = dict(params or {})
    params.setdefault("page_size", 100)
    page, count = 1, 0
    while True:
        params["page"] = page
        body = api_get(session, url, params, verbose=verbose)
        if not body:
            return
        items = body.get("data", [])
        if not items:
            return
        for item in items:
            yield item
            count += 1
            if limit and count >= limit:
                return
        if not body.get("links", {}).get("next"):
            return
        page += 1
        time.sleep(0.05)


# --------------------------------------------------------------------------- #
# Domain logic
# --------------------------------------------------------------------------- #

def load_biome_lineages(csv_path: Path) -> Dict[str, str]:
    """Return shortname -> lineage mapping from biomes.csv."""
    if not csv_path.exists():
        return {}
    df = pd.read_csv(csv_path)
    out = {}
    for _, row in df.iterrows():
        key = row["Environment"].lower().replace(" ", "_").replace("-", "_")
        out[key] = row["API_Lineage"]
    return out


def resolve_biomes(biome_args: List[str], lineage_arg: Optional[str],
                   label_arg: Optional[str], biomes_csv: Path
                   ) -> List[tuple]:
    """Return list of (label, lineage) to process."""
    if lineage_arg:
        return [(label_arg or lineage_arg.split(":")[-1].lower(), lineage_arg)]
    table = load_biome_lineages(biomes_csv)
    pairs = []
    for name in biome_args:
        if name not in table:
            raise SystemExit(f"biome '{name}' not in {biomes_csv}. "
                             f"Run list_biomes.py first or use --lineage.")
        pairs.append((name, table[name]))
    return pairs


def fetch_samples(session, lineage: str, max_samples: int, verbose: bool):
    encoded = urllib.parse.quote(lineage, safe="")
    url = f"{API}/biomes/{encoded}/samples"
    params = {"experiment_type": "metagenomic", "ordering": "accession"}
    return list(paginate(session, url, params, limit=max_samples, verbose=verbose))


def fetch_taxonomy_for_analysis(session, analysis_id: str, verbose: bool):
    """Try the SSU endpoint first, then fall back to other taxonomy endpoints."""
    for ep in TAXONOMY_ENDPOINTS:
        body = api_get(session, f"{API}/analyses/{analysis_id}/{ep}", verbose=verbose)
        if body and body.get("data"):
            return body["data"]
    return None


def aggregate_taxa(entries: List[dict]) -> Dict[str, Counter]:
    """Roll up an analysis's taxonomy entries into Counter per rank."""
    out: Dict[str, Counter] = defaultdict(Counter)
    for entry in entries:
        attrs = entry.get("attributes", {})
        count = (attrs.get("count") or attrs.get("abundance")
                 or attrs.get("reads") or attrs.get("frequency") or 0)
        try:
            count = int(count)
        except (TypeError, ValueError):
            continue
        if count <= 0:
            continue
        rank, name = attrs.get("rank", ""), (attrs.get("name") or "").strip()
        hierarchy = attrs.get("hierarchy", {}) or {}
        for target in RANKS:
            taxon = hierarchy.get(target) or (name if rank == target else None)
            if not taxon:
                continue
            taxon = taxon.strip()
            if len(taxon) > 1 and taxon.lower() not in SKIP_NAMES:
                out[target][taxon] += count
    return out


def process_biome(session, label: str, lineage: str, max_samples: int,
                  verbose: bool) -> Dict[str, Counter]:
    print(f"\n=== {label}  ({lineage}) ===")
    samples = fetch_samples(session, lineage, max_samples, verbose)
    print(f"  {len(samples)} samples")

    totals: Dict[str, Counter] = defaultdict(Counter)
    n_analyses = 0
    for i, sample in enumerate(samples, 1):
        sid = sample["id"]
        for run in paginate(session, f"{API}/samples/{sid}/runs", verbose=verbose):
            rid = run["id"]
            for analysis in paginate(session, f"{API}/runs/{rid}/analyses", verbose=verbose):
                entries = fetch_taxonomy_for_analysis(session, analysis["id"], verbose)
                if not entries:
                    continue
                for rank, c in aggregate_taxa(entries).items():
                    totals[rank].update(c)
                n_analyses += 1
        if i % 10 == 0 or i == len(samples):
            print(f"  sample {i}/{len(samples)}  analyses={n_analyses}  "
                  f"unique-genera={len(totals.get('genus', ()))}")
    return totals


def write_results(label: str, totals: Dict[str, Counter], out_dir: Path, top_n: int):
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for rank, counter in totals.items():
        if not counter:
            continue
        rows = counter.most_common(top_n)
        df = pd.DataFrame(rows, columns=["taxon", "total_abundance"])
        df.insert(0, "rank", range(1, len(df) + 1))
        df["taxa_level"] = rank
        df = df[["rank", "taxon", "taxa_level", "total_abundance"]]
        fp = out_dir / f"{label}__{rank}__top{top_n}.csv"
        df.to_csv(fp, index=False)
        written.append(fp)
    return written


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--biomes", nargs="+", default=[],
                    help="Biome short-names (resolved via biomes.csv).")
    ap.add_argument("--lineage", default=None,
                    help="Alternative: raw MGnify lineage "
                         "(e.g. 'root:Environmental:Aquatic:Marine').")
    ap.add_argument("--label", default=None,
                    help="Output prefix when using --lineage.")
    ap.add_argument("--max-samples", type=int, default=100,
                    help="Cap samples fetched per biome (default: 100).")
    ap.add_argument("--top-n", type=int, default=100,
                    help="Top-N taxa per rank to write (default: 100).")
    ap.add_argument("--biomes-csv", default=str(HERE / "biomes.csv"),
                    help="Lineage table from list_biomes.py.")
    ap.add_argument("--out-root", default=str(HERE / "results"),
                    help="Results root (a timestamped subdir is created here).")
    ap.add_argument("--run-tag", default=None,
                    help="Override the timestamped output dir name.")
    ap.add_argument("--verbose", action="store_true",
                    help="Log HTTP and pagination details to stderr.")
    args = ap.parse_args()

    if not args.biomes and not args.lineage:
        ap.error("provide --biomes <name> [<name> ...] or --lineage <path>")

    pairs = resolve_biomes(args.biomes, args.lineage, args.label,
                           Path(args.biomes_csv))
    run_tag = args.run_tag or datetime.now().strftime("%Y-%m-%d_%H%M%S")
    out_dir = Path(args.out_root) / run_tag
    print(f"Output -> {out_dir}")

    t0 = time.time()
    with make_session() as session:
        for label, lineage in pairs:
            totals = process_biome(session, label, lineage,
                                   args.max_samples, args.verbose)
            files = write_results(label, totals, out_dir, args.top_n)
            print(f"  wrote {len(files)} CSV(s) for {label}")

    print(f"\nDone in {time.time() - t0:.1f}s -> {out_dir}")


if __name__ == "__main__":
    main()
