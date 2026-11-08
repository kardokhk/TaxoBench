"""The TaxoBench data model (docs/taxobench.md, section Data model).

Seven long tables, TSV with a header, UTF-8, tab separated, empty string for a
missing value. Each file has a JSON sidecar carrying its sha256. A table may be
split over several files named ``<table>.tsv`` or ``<table>.<part>.tsv`` (for
example ``truth.gg2.tsv``); the reader concatenates them. Extra columns are
kept and ignored by the scorer.

Required columns (aliases accepted on reading in brackets):

references   ref_id, release
taxa         ref_id, rank, taxon_id, name, parent_taxon_id
queries      query_id, tier, dataset, region, split, cluster_id [cluster]
abundance    query_id, sample_id, count
truth        query_id, ref_id, rank, truth_set, status
runs         run_id, classifier, version, ref_id
predictions  run_id, query_id, rank, raw_label, taxon_ids, confidence, call

Recommended columns: references: retrieval_date, taxonomy_sha256,
sequences_sha256; queries: seq_sha1; truth: provenance (or builder columns);
runs: model_md5, parameters (JSON), thresholds (JSON), command, date,
code_commit, adapter, adapter_version.

Vocabularies:
- rank: domain, phylum, class, order, family, genus, species;
- split: calib, eval;
- truth status: known, novel, unresolvable, unscorable;
- prediction call: named, group, nocall, oov.

Multi-valued fields (truth_set, taxon_ids) join taxon_id values with ``|``, so
a taxon_id must not contain ``|``. ``truth_set`` is empty for novel and
unscorable rows; ``taxon_ids`` is empty for nocall and oov.

Sidecar: ``<file>.json`` (also accepted: ``<file>.sha256.json``, or a
``MANIFEST.json`` / ``manifest.json`` in the directory mapping file names to
objects) holding at least ``{"sha256": "<hex>"}``. The reader verifies it and
refuses a file whose sidecar is missing or does not match, unless asked not to.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from . import RANKS

SCHEMA_VERSION = "1"
STATUSES = ("known", "novel", "unresolvable", "unscorable")
CALLS = ("named", "group", "nocall", "oov")
SPLITS = ("calib", "eval")

REQUIRED = {
    "references": ["ref_id", "release"],
    "taxa": ["ref_id", "rank", "taxon_id", "name", "parent_taxon_id"],
    "queries": ["query_id", "tier", "dataset", "region", "split", "cluster_id"],
    "abundance": ["query_id", "sample_id", "count"],
    "truth": ["query_id", "ref_id", "rank", "truth_set", "status"],
    "runs": ["run_id", "classifier", "version", "ref_id"],
    "predictions": ["run_id", "query_id", "rank", "raw_label", "taxon_ids",
                    "confidence", "call"],
}
KEYS = {
    "references": ["ref_id"],
    "taxa": ["ref_id", "taxon_id"],
    "queries": ["query_id"],
    "abundance": ["query_id", "sample_id"],
    "truth": ["query_id", "ref_id", "rank"],
    "runs": ["run_id"],
    "predictions": ["run_id", "query_id", "rank"],
}
ALIASES = {"queries": {"cluster": "cluster_id"}}


class ContractError(ValueError):
    """A table violates the data model. Always fatal for scoring."""


def sha256_file(fp) -> str:
    h = hashlib.sha256()
    with open(fp, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_sidecar(fp, generator="taxobench", **extra):
    fp = Path(fp)
    n = sum(1 for _ in open(fp, encoding="utf-8")) - 1
    meta = {"file": fp.name, "sha256": sha256_file(fp), "rows": max(n, 0),
            "schema_version": SCHEMA_VERSION, "generator": generator}
    meta.update(extra)
    with open(str(fp) + ".json", "w") as fh:
        json.dump(meta, fh, indent=1, sort_keys=True)
    return meta


def _sidecar_sha(fp: Path):
    for cand in (Path(str(fp) + ".json"), Path(str(fp) + ".sha256.json")):
        if cand.exists():
            with open(cand) as fh:
                d = json.load(fh)
            if "sha256" in d:
                return d["sha256"], cand
    for name in ("MANIFEST.json", "manifest.json"):
        m = fp.parent / name
        if m.exists():
            with open(m) as fh:
                d = json.load(fh)
            files = d.get("files", d)
            entry = files.get(fp.name) if isinstance(files, dict) else None
            if isinstance(entry, dict) and "sha256" in entry:
                return entry["sha256"], m
            if isinstance(entry, str):
                return entry, m
    return None, None


def table_files(root, table) -> list[Path]:
    root = Path(root)
    files = sorted(set(root.glob(f"{table}.tsv")) | set(root.glob(f"{table}.*.tsv")))
    return [f for f in files if f.is_file()]


def write_table(df: pd.DataFrame, fp, table=None, generator="taxobench"):
    """Write a contract table (columns validated when `table` is given) plus sidecar."""
    if table is not None:
        missing = [c for c in REQUIRED[table] if c not in df.columns]
        if missing:
            raise ContractError(f"{table}: missing columns {missing}")
    fp = Path(fp)
    fp.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(fp, sep="\t", index=False, na_rep="")
    return write_sidecar(fp, generator=generator, table=table or fp.stem)


def read_table(root_or_file, table, verify=True, required=True) -> pd.DataFrame:
    """Read every file of `table` under a directory (or one file) as strings.

    Fails loudly on: no file (when required), a missing or mismatching sha256
    sidecar (when verify), missing columns, duplicated keys, bad vocabulary.
    """
    p = Path(root_or_file)
    files = [p] if p.is_file() else table_files(p, table)
    if not files:
        if required:
            raise ContractError(f"no {table}.tsv or {table}.*.tsv under {p}")
        return pd.DataFrame(columns=REQUIRED[table])
    frames = []
    for fp in files:
        if verify:
            want, src = _sidecar_sha(fp)
            if want is None:
                raise ContractError(f"{fp}: no sha256 sidecar ({fp.name}.json or MANIFEST.json)")
            got = sha256_file(fp)
            if got != want:
                raise ContractError(f"{fp}: sha256 {got} does not match {src} ({want})")
        df = pd.read_csv(fp, sep="\t", dtype=str, keep_default_na=False,
                         na_values=[], quoting=3)
        frames.append(df.rename(columns=ALIASES.get(table, {})))
    df = pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]
    validate(df, table, where=", ".join(str(f) for f in files))
    return df


def validate(df: pd.DataFrame, table: str, where=""):
    missing = [c for c in REQUIRED[table] if c not in df.columns]
    if missing:
        raise ContractError(f"{table} ({where}): missing columns {missing}")
    key = KEYS[table]
    dup = df.duplicated(key, keep=False)
    if dup.any():
        ex = df.loc[dup, key].head(3).to_dict("records")
        raise ContractError(f"{table} ({where}): {int(dup.sum())} rows with duplicated key {key}, e.g. {ex}")
    if "rank" in df.columns and table in ("taxa", "truth", "predictions"):
        bad = set(df["rank"]) - set(RANKS)
        if bad:
            raise ContractError(f"{table} ({where}): unknown rank(s) {sorted(bad)[:5]}")
    if table == "truth":
        bad = set(df["status"]) - set(STATUSES)
        if bad:
            raise ContractError(f"truth ({where}): unknown status {sorted(bad)}")
        named = df["status"].isin(["known", "unresolvable"])
        empty = named & (df["truth_set"] == "")
        if empty.any():
            raise ContractError(f"truth ({where}): {int(empty.sum())} known/unresolvable rows "
                                "with an empty truth_set")
        multi = df["truth_set"].str.contains("|", regex=False)
        if (multi & (df["status"] == "known")).any():
            raise ContractError(f"truth ({where}): status known with several taxa; use unresolvable")
    if table == "predictions":
        bad = set(df["call"]) - set(CALLS)
        if bad:
            raise ContractError(f"predictions ({where}): unknown call {sorted(bad)}")
        has = df["taxon_ids"] != ""
        if (has & df["call"].isin(["nocall", "oov"])).any() or \
           (~has & df["call"].isin(["named", "group"])).any():
            raise ContractError(f"predictions ({where}): taxon_ids inconsistent with call")
    if table == "queries":
        bad = set(df["split"]) - set(SPLITS)
        if bad:
            raise ContractError(f"queries ({where}): unknown split {sorted(bad)}")
    if table == "abundance" and len(df):
        c = pd.to_numeric(df["count"], errors="coerce")
        if c.isna().any() or (c < 0).any():
            raise ContractError(f"abundance ({where}): count must be a non-negative number")


def split_ids(s: str) -> tuple:
    return tuple(x for x in s.split("|") if x) if s else ()


@dataclass
class TierTables:
    """The tables of one tier directory (truth side)."""
    root: Path
    references: pd.DataFrame
    taxa: pd.DataFrame
    queries: pd.DataFrame
    truth: pd.DataFrame
    abundance: pd.DataFrame

    @classmethod
    def load(cls, root, verify=True):
        root = Path(root)
        return cls(root=root,
                   references=read_table(root, "references", verify),
                   taxa=read_table(root, "taxa", verify),
                   queries=read_table(root, "queries", verify),
                   truth=read_table(root, "truth", verify),
                   abundance=read_table(root, "abundance", verify, required=False))


@dataclass
class RunTables:
    """runs + predictions (classifier side)."""
    root: Path
    runs: pd.DataFrame
    predictions: pd.DataFrame

    @classmethod
    def load(cls, root, verify=True):
        root = Path(root)
        return cls(root=root, runs=read_table(root, "runs", verify),
                   predictions=read_table(root, "predictions", verify))
