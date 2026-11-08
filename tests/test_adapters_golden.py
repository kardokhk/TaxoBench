"""Golden parser tests: raw outputs of the real tools (biome_digestive_system
V5-V7 of the paper benchmark, all eight classifiers) -> contract predictions.

The expected tables are regenerated only with TAXOBENCH_REGEN_GOLDEN=1; any
other difference fails. A change in an adapter's output needs a new adapter
version and a reviewed regeneration.
"""
import os
from pathlib import Path

import pandas as pd
import pytest

from taxobench.adapters import get_adapter
from taxobench.taxonomy import load_reference_lineages

F = Path(__file__).parent / "fixtures" / "adapters" / "bds_V5-V7"
RAW = {
    "taxonova": "raw/taxonova_raw.tsv", "idtaxa": "raw/idtaxa_raw.tsv",
    "spingo": "raw/spingo_raw.tsv", "mapseq": "raw/mapseq_raw.mseq",
    "kraken2": {"output": "raw/kraken2_classification.txt", "report": "raw/kraken2_report.txt"},
    "blast": "raw/blast_classification.tsv", "vsearch": "raw/vsearch_classification.tsv",
    "naive-bayes": "raw/naive-bayes_classification.tsv",
}
PARAMS = {"mapseq": {"cutoff": 0.4}}


@pytest.fixture(scope="module")
def ref():
    return load_reference_lineages(F / "reference_taxonomy.tsv", "fx/V5-V7")


def produce(name, ref):
    raw = RAW[name]
    raw = {k: F / v for k, v in raw.items()} if isinstance(raw, dict) else F / raw
    ad = get_adapter(name)
    q = [x.strip() for x in open(F / "query_ids.txt") if x.strip()]
    return ad.to_predictions(ad.parse(raw, **PARAMS.get(name, {})), ref, f"fx:{name}", query_ids=q)


@pytest.mark.parametrize("name", sorted(RAW))
def test_golden(name, ref):
    got = produce(name, ref)
    fp = F / "expected" / f"{name}.tsv"
    if os.environ.get("TAXOBENCH_REGEN_GOLDEN") == "1":
        fp.parent.mkdir(exist_ok=True)
        got.to_csv(fp, sep="\t", index=False)
    want = pd.read_csv(fp, sep="\t", dtype=str, keep_default_na=False)
    got = pd.read_csv(pd.io.common.StringIO(got.to_csv(sep="\t", index=False)), sep="\t",
                      dtype=str, keep_default_na=False)
    pd.testing.assert_frame_equal(got, want)


def _row(df, qid, rank):
    r = df[(df.query_id == qid) & (df["rank"] == rank)]
    assert len(r) == 1
    return r.iloc[0]


def test_spingo_space_names_repaired(ref):
    p = produce("spingo", ref)
    r = _row(p, "38cfbca9a08f8ca7aa7550eee2a4f362", "species")
    assert r.raw_label == "s__Blautia A 142742 liquoris"
    assert r.call == "named"
    assert ref.name[r.taxon_ids] == "Blautia_A_142742 liquoris"
    # upper ranks completed from the reference path of the genus
    assert _row(p, "38cfbca9a08f8ca7aa7550eee2a4f362", "family").call == "named"


def test_kraken2_suffixless_genus(ref):
    p = produce("kraken2", ref)
    g = p[p["rank"] == "genus"]
    oov = g[g.call == "oov"]
    assert len(oov) and oov.note.str.contains("ambiguous").all()
    # every repaired genus is a unique reference name
    g = g[g.call == "named"]
    named = g[g.raw_label.str[3:] != g.taxon_ids.map(lambda t: ref.name[t])]
    assert len(named) > 0  # at least one suffix repair happened
    for _, r in named.iterrows():
        assert ref.name[r.taxon_ids].startswith(r.raw_label[3:])
    assert get_adapter("kraken2").confidence_kind == "run_share"


def test_mapseq_is_hierarchical(ref):
    p = produce("mapseq", ref)
    w = p.pivot(index="query_id", columns="rank", values="call")
    order = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    for qid, row in w.iterrows():
        seen_nocall = False
        for rk in order:
            if row[rk] == "nocall":
                seen_nocall = True
            else:
                assert not seen_nocall, (qid, rk)
