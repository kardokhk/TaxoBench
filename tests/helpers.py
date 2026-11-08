"""Build small contract tables by hand for known-answer tests."""
import pandas as pd

from taxobench import RANKS
from taxobench.adapters.base import Adapter, Parsed
from taxobench.metrics import align
from taxobench.taxonomy import Reference, is_nocall, parse_lineage

LINEAGES = [
    "d__B; p__P1; c__C1; o__O1; f__F1; g__G1; s__G1 a",
    "d__B; p__P1; c__C1; o__O1; f__F1; g__G1; s__G1 b",
    "d__B; p__P1; c__C1; o__O1; f__F1; g__G2; s__G2 c",
    "d__B; p__P1; c__C1; o__O1; f__F2; g__G3; s__G3 d",
]
REF_ID = "toy"


def toy_ref(lineages=LINEAGES, ref_id=REF_ID):
    return Reference.from_lineages(lineages, ref_id)


def truth_table(ref, truth: dict, ref_id=REF_ID):
    """truth: {query_id: lineage string or "other"}. Status per rank as in
    the contract: unscorable for other/empty, novel where unnamed or absent."""
    rows = []
    for q, lab in truth.items():
        lin = {} if lab == "other" else parse_lineage(lab)
        res = ref.resolve_lineage(lin) if lin else None
        for r in RANKS:
            if not lin:
                st, ts = "unscorable", ""
            elif is_nocall(lin.get(r, "")):
                st, ts = "novel", ""
            elif res[r].call == "named":
                st, ts = "known", "|".join(res[r].taxon_ids)
            elif res[r].call == "group":
                st, ts = "unresolvable", "|".join(res[r].taxon_ids)
            else:
                st, ts = "novel", ""
            rows.append(dict(query_id=q, ref_id=ref_id, rank=r, truth_set=ts, status=st))
    return pd.DataFrame(rows)


def queries_table(qids, clusters=None, split="eval", tier="TX", dataset="d", region="V4"):
    clusters = clusters or {q: q for q in qids}
    return pd.DataFrame([dict(query_id=q, tier=tier, dataset=dataset, region=region,
                              split=split if isinstance(split, str) else split[q],
                              cluster_id=clusters[q]) for q in qids])


def preds_table(ref, calls: dict, run_id="r1", conf=None):
    """calls: {query_id: prefixed lineage string ('' = no call)}."""
    parsed = []
    for q, lab in calls.items():
        lin = {r: t for r, t in parse_lineage(lab).items()} if lab else {}
        c = {r: (conf or {}).get(q, 0.9) for r in lin}
        parsed.append(Parsed(q, lin, c, lab))
    return Adapter().to_predictions(parsed, ref, run_id, query_ids=list(calls))


def aligned(ref, truth: dict, calls: dict, **kw):
    t = truth_table(ref, truth)
    q = queries_table(list(truth), **kw)
    p = preds_table(ref, calls)
    return align(t, q, p, ref, list(truth))
