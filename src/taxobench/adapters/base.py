"""Adapter interface: raw tool output -> contract `predictions` rows.

An adapter is versioned. Changing what it emits for the same raw input needs a
new ``version`` and a golden-test update (tests/adapters/). Parsing produces a
`Parsed` record per query; resolution against the run's reference is shared
(`to_predictions`) so that every tool's names are validated the same way.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable

import pandas as pd

from .. import PREFIX, RANKS
from ..contract import ContractError
from ..taxonomy import Reference, Resolved, is_nocall

#: What the confidence column means. Only "probability" enters calibration
#: (M9); "probability" and "score" enter threshold-free curves (M8).
CONFIDENCE_KINDS = ("probability", "score", "run_share", "none")


@dataclass
class Parsed:
    query_id: str
    lineage: dict                 # rank -> token (no prefix); missing = no call
    confidence: dict = field(default_factory=dict)   # rank -> float
    raw_label: str = ""


class Adapter:
    name = "base"
    version = "0"
    tool = ""
    confidence_kind = "none"
    confidence_note = ""
    #: fill ranks above the deepest resolved rank from the reference path
    complete_from_reference = False

    def parse(self, raw: Path | dict, **params) -> list[Parsed]:
        raise NotImplementedError

    def repair(self, rank: str, name: str, ref: Reference) -> Iterable[str]:
        """Alternative spellings to try when `name` is not in the reference.
        Accepted only if they lead to exactly one taxon (Reference._lookup)."""
        return ()

    def describe(self) -> dict:
        return {"adapter": self.name, "adapter_version": self.version,
                "confidence_kind": self.confidence_kind,
                "confidence_note": self.confidence_note}

    # ------------------------------------------------------------------
    def to_predictions(self, parsed: list[Parsed], ref: Reference, run_id: str,
                       query_ids: Iterable[str] | None = None) -> pd.DataFrame:
        """Resolve parsed calls and build the predictions table.

        Fails closed: duplicated query IDs, and (when query_ids is given) any
        missing or unexpected ID raise ContractError.
        """
        seen = {}
        for p in parsed:
            if p.query_id in seen:
                raise ContractError(f"{self.name}: query {p.query_id!r} appears twice in the output")
            seen[p.query_id] = p
        if query_ids is not None:
            want = list(query_ids)
            missing = set(want) - set(seen)
            extra = set(seen) - set(want)
            if missing or extra:
                raise ContractError(
                    f"{self.name}: output covers {len(seen)} IDs against {len(want)} queries "
                    f"({len(missing)} missing, {len(extra)} unexpected)")
        repair = self.repair if type(self).repair is not Adapter.repair else None
        rows = []
        for qid, p in seen.items():
            res = ref.resolve_lineage(p.lineage, repair=repair)
            if self.complete_from_reference:
                res = _complete(res, ref)
            for r in RANKS:
                rv: Resolved = res[r]
                tok = p.lineage.get(r, "")
                conf = p.confidence.get(r)
                if rv.call in ("nocall",) or conf is None or (isinstance(conf, float) and math.isnan(conf)):
                    conf_s = ""
                else:
                    conf_s = repr(float(conf))
                rows.append((run_id, qid, r, ("" if is_nocall(tok) else PREFIX[r] + tok),
                             "|".join(rv.taxon_ids), conf_s, rv.call, rv.note))
        return pd.DataFrame(rows, columns=["run_id", "query_id", "rank", "raw_label",
                                           "taxon_ids", "confidence", "call", "note"])


def _complete(res: dict, ref: Reference) -> dict:
    """Fill no-call ranks above the deepest named rank from the reference
    path of that taxon (only when it is a single taxon)."""
    deepest = None
    for r in RANKS:
        if res[r].call == "named" and len(res[r].taxon_ids) == 1:
            deepest = res[r].taxon_ids[0]
    if deepest is None:
        return res
    out = dict(res)
    for tid in ref.path(deepest):
        r = ref.rank[tid]
        if out[r].call == "nocall":
            out[r] = Resolved("named", (tid,), "completed from reference")
    return out


def to_float(x) -> float:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return float("nan")
    return v


def split_prefixed(label: str) -> tuple[dict, list]:
    """Prefixed lineage -> ({rank: token}, [ranks in order]); tolerant of a
    trailing empty rank such as ``s__``."""
    out, order = {}, []
    for part in (label or "").split(";"):
        part = part.strip()
        if not part:
            continue
        pre = part[:3]
        if pre in PREFIX.values():
            r = RANKS[list(PREFIX.values()).index(pre)]
            out[r] = part[3:].strip()
            order.append(r)
        elif len(part) >= 3 and part[1:3] == "__":
            raise ContractError(f"unknown rank prefix {pre!r} in {label!r}")
        elif is_nocall(part):
            continue
        else:
            raise ContractError(f"token without rank prefix {part!r} in {label!r}")
    return out, order


REGISTRY: dict[str, Callable[[], Adapter]] = {}


def register(cls):
    REGISTRY[cls.name] = cls
    return cls


def get_adapter(name: str) -> Adapter:
    from . import tools  # noqa: F401  (registers the built-in adapters)
    try:
        from importlib.metadata import entry_points
        for ep in entry_points(group="taxobench.adapters"):
            if ep.name not in REGISTRY:
                REGISTRY[ep.name] = ep.load()
    except Exception:
        pass
    if name not in REGISTRY:
        raise KeyError(f"unknown adapter {name!r}; known: {sorted(REGISTRY)}")
    return REGISTRY[name]()
