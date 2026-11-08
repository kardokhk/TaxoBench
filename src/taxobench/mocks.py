"""Mock community tier: convert the output of the mock builders (mocks/scripts/)
into tier tables.

The builders write, per mock and region, a truth table
``<data>/mock_community/<mock>/<region>/trueish*.tsv`` (query id, lineage) and
the query sequences ``rep_seqs_<mock>.fasta`` below the same directory. This
module writes the seven-table tier (references, taxa, queries, truth and an
empty abundance table) and one query FASTA per cell, whose IDs are the tier's
query IDs, so a classifier run on that FASTA can be scored directly.

Truth status at rank r:
- unscorable when the truth label is `other`, `Unassigned` or names no rank;
- novel when the truth names no taxon at r, or names one the region's
  reference does not contain;
- unresolvable when the name stands for a group of reference taxa;
- known otherwise.

Every query is in split `eval` and is its own cluster (mock queries have no
per-sample read counts, so the abundance table is empty).
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from . import RANKS
from .contract import write_table
from .taxonomy import Reference, is_nocall, load_reference_lineages, parse_lineage

#: reference taxonomy file under the reference root, per region
REF_FILES = {"RRN": ("RRN", "grond_clean.tsv"), "V1-V9": ("V1-V9", "taxonomy.tsv"),
             "V3-V4": ("V3-V4", "taxonomy.tsv"), "V4": ("V4", "taxonomy.tsv"),
             "V5-V7": ("V5-V7", "taxonomy.tsv")}
#: reference id per region (Greengenes2 2024.09 for 16S regions, GROND r207 for the operon)
REF_IDS = {"RRN": "grond-r207/RRN", "V1-V9": "gg2-2024.09/V1-V9", "V3-V4": "gg2-2024.09/V3-V4",
           "V4": "gg2-2024.09/V4", "V5-V7": "gg2-2024.09/V5-V7"}
SEP = "__"


def query_id(mock: str, region: str, qid: str) -> str:
    return f"{mock}{SEP}{region}{SEP}{qid}"


def find_truth(mock_root: Path, mock: str, region: str):
    d = mock_root / mock / region
    for n in ("trueish-taxonomies.tsv", f"trueish_taxonomy_{mock}.tsv", f"trueish-taxonomy_{mock}.tsv"):
        if (d / n).exists():
            return d / n
    return None


def find_fasta(mock_root: Path, mock: str, region: str):
    hits = sorted((mock_root / mock / region).rglob(f"rep_seqs_{mock}.fasta"))
    return hits[0] if hits else None


def read_truth_tsv(fp) -> dict:
    out = {}
    with open(fp) as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if parts[0] in ("id", "Feature ID", "taxonomy", "Taxon"):
                continue
            if len(parts) >= 2:
                if parts[0] in out:
                    raise ValueError(f"{fp}: duplicated ID {parts[0]}")
                out[parts[0]] = parts[1]
    return out


def read_fasta(fp) -> dict:
    seqs, cur = {}, None
    with open(fp) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                cur = line[1:].split()[0]
                if cur in seqs:
                    raise ValueError(f"{fp}: duplicated ID {cur}")
                seqs[cur] = []
            elif cur is not None and line:
                seqs[cur].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def build(data_dir, ref_root, out_dir, cells=None, tier="mocks", verbose=print):
    """Write the tier tables and per-cell query FASTA files for every
    (mock, region) under <data_dir>/mock_community that has a truth table and
    query sequences, or for the given list of (mock, region) pairs."""
    mock_root, ref_root, out = Path(data_dir) / "mock_community", Path(ref_root), Path(out_dir)
    if cells is None:
        cells = sorted((m.name, r.name) for m in mock_root.iterdir() if m.is_dir()
                       for r in m.iterdir() if r.is_dir() and r.name in REF_IDS)
    refs, taxa_rows, q_rows, t_rows, r_rows = {}, [], [], [], []
    (out / "queries").mkdir(parents=True, exist_ok=True)
    n_cells = 0
    for mock, region in cells:
        if region not in REF_IDS:
            raise ValueError(f"unknown region {region!r} (known: {', '.join(REF_IDS)})")
        truth_fp, fasta_fp = find_truth(mock_root, mock, region), find_fasta(mock_root, mock, region)
        if truth_fp is None or fasta_fp is None:
            verbose(f"skip {mock} {region}: no {'truth table' if truth_fp is None else 'query FASTA'}")
            continue
        truth, seqs = read_truth_tsv(truth_fp), read_fasta(fasta_fp)
        if set(truth) != set(seqs):
            raise ValueError(f"{mock} {region}: truth and FASTA IDs differ "
                             f"({len(set(truth) - set(seqs))} only in truth, "
                             f"{len(set(seqs) - set(truth))} only in FASTA)")
        ref_id = REF_IDS[region]
        if ref_id not in refs:
            sub, name = REF_FILES[region]
            refs[ref_id] = load_reference_lineages(ref_root / sub / name, ref_id)
            r_rows.append({"ref_id": ref_id, "release": ref_id.split("/")[0],
                           "region": region, "source": f"{sub}/{name}"})
        ref: Reference = refs[ref_id]
        with open(out / "queries" / f"{mock}{SEP}{region}.fasta", "w") as fh:
            for qid in sorted(seqs):
                fh.write(f">{query_id(mock, region, qid)}\n{seqs[qid]}\n")
        for qid in sorted(truth):
            lab, qq = truth[qid], query_id(mock, region, qid)
            q_rows.append({"query_id": qq, "tier": tier, "dataset": mock, "region": region,
                           "split": "eval", "cluster_id": qq, "source_id": qid})
            try:
                lin = parse_lineage(lab)
            except Exception:           # noqa: BLE001  (an unparsable label is unscorable)
                lin = {}
            unscorable = not lin or lab.strip().lower() in ("other", "unassigned")
            res = ref.resolve_lineage(lin) if not unscorable else None
            for r in RANKS:
                if unscorable:
                    st, ts = "unscorable", ""
                elif is_nocall(lin.get(r, "")):
                    st, ts = "novel", ""
                elif res[r].call == "named":
                    st, ts = "known", "|".join(res[r].taxon_ids)
                elif res[r].call == "group":
                    st, ts = "unresolvable", "|".join(res[r].taxon_ids)
                else:
                    st, ts = "novel", ""
                t_rows.append({"query_id": qq, "ref_id": ref_id, "rank": r, "truth_set": ts, "status": st,
                               "provenance": "mock truth" + ("" if st != "novel" or is_nocall(lin.get(r, ""))
                                                             else "; name absent from reference")})
        n_cells += 1
    if not q_rows:
        raise ValueError(f"no mock cell with a truth table and query FASTA under {mock_root}")
    for ref_id, ref in refs.items():
        for tid in ref.rank:
            taxa_rows.append({"ref_id": ref_id, "rank": ref.rank[tid], "taxon_id": tid,
                              "name": ref.name[tid], "parent_taxon_id": ref.parent[tid],
                              "is_group": int(ref.is_group[tid])})
    write_table(pd.DataFrame(r_rows), out / "references.tsv", "references")
    write_table(pd.DataFrame(taxa_rows), out / "taxa.tsv", "taxa")
    write_table(pd.DataFrame(q_rows), out / "queries.tsv", "queries")
    write_table(pd.DataFrame(t_rows), out / "truth.tsv", "truth")
    write_table(pd.DataFrame(columns=["query_id", "sample_id", "count"]), out / "abundance.tsv", "abundance")
    verbose(f"{tier}: {len(q_rows)} queries in {n_cells} cells; query FASTA in {out / 'queries'}")
    return out
