"""`taxobench import-mocks`: mock builder output -> tier tables and query FASTA,
then a run on that FASTA ingested and scored."""
import pandas as pd
import pytest

from helpers import LINEAGES
from taxobench import cli
from taxobench.contract import read_table
from taxobench.mocks import REF_IDS, query_id

TRUTH = {"asv1": LINEAGES[0], "asv2": LINEAGES[2],
         "asv3": "d__B; p__P1; c__C1; o__O1; f__F1; g__G9; s__G9 z",   # genus absent from the reference
         "asv4": "other"}


@pytest.fixture()
def mock_tree(tmp_path):
    cell = tmp_path / "data" / "mock_community" / "m1" / "V4"
    (cell / "source" / "single").mkdir(parents=True)
    (cell / "trueish-taxonomies.tsv").write_text("id\ttaxonomy\n" + "".join(f"{k}\t{v}\n" for k, v in TRUTH.items()))
    (cell / "source" / "single" / "rep_seqs_m1.fasta").write_text("".join(f">{k}\nACGTACGT\n" for k in TRUTH))
    ref = tmp_path / "ref" / "V4"
    ref.mkdir(parents=True)
    (ref / "taxonomy.tsv").write_text("Feature ID\tTaxon\n" + "".join(f"r{i}\t{x}\n" for i, x in enumerate(LINEAGES)))
    return tmp_path


def _status(truth, qid, rank):
    row = truth[(truth.query_id == query_id("m1", "V4", qid)) & (truth["rank"] == rank)]
    return row.status.iloc[0]


def test_import_writes_tier_and_fasta(mock_tree):
    out = mock_tree / "tier"
    assert cli.main(["import-mocks", "--data-dir", str(mock_tree / "data"), "--ref-dir", str(mock_tree / "ref"),
                     "--out", str(out)]) == 0
    assert cli.main(["validate", "--tier", str(out)]) == 0
    q = read_table(out, "queries")
    assert set(q.query_id) == {query_id("m1", "V4", k) for k in TRUTH}
    assert set(q.tier) == {"mocks"} and set(q.split) == {"eval"}
    t = read_table(out, "truth")
    assert set(t.ref_id) == {REF_IDS["V4"]}
    assert _status(t, "asv1", "species") == "known"
    assert _status(t, "asv3", "family") == "known" and _status(t, "asv3", "genus") == "novel"
    assert _status(t, "asv4", "domain") == "unscorable"
    fasta = (out / "queries" / "m1__V4.fasta").read_text()
    assert fasta.count(">") == 4 and f">{query_id('m1', 'V4', 'asv2')}\n" in fasta


def test_imported_tier_scores_a_run(mock_tree):
    out = mock_tree / "tier"
    cli.main(["import-mocks", "--data-dir", str(mock_tree / "data"), "--ref-dir", str(mock_tree / "ref"),
              "--out", str(out)])
    raw = mock_tree / "raw.tsv"
    raw.write_text("Feature ID\tTaxon\tConfidence\n" + "".join(
        f"{query_id('m1', 'V4', k)}\t{v if v != 'other' else 'Unassigned'}\t0.9\n" for k, v in TRUTH.items()))
    assert cli.main(["ingest", "--adapter", "taxonomy-tsv", "--version", "1", "--ref-id", REF_IDS["V4"],
                     "--taxa", str(out), "--query-fasta", str(out / "queries" / "m1__V4.fasta"),
                     "--raw", str(raw), "--run-id", "R", "--out", str(mock_tree / "runs" / "R")]) == 0
    assert cli.main(["score", "--tier", str(out), "--runs", str(mock_tree / "runs" / "R"),
                     "--out", str(mock_tree / "sc"), "--n-boot", "50"]) == 0
    s = pd.read_csv(mock_tree / "sc" / "scores.tsv", sep="\t")
    f1 = s[(s["rank"] == "species") & (s.rule == "G1") & (s.metric == "f1")].value.iloc[0]
    assert f1 > 0


def test_ids_must_match(mock_tree):
    fa = mock_tree / "data" / "mock_community" / "m1" / "V4" / "source" / "single" / "rep_seqs_m1.fasta"
    fa.write_text(">asv1\nACGT\n")
    assert cli.main(["import-mocks", "--data-dir", str(mock_tree / "data"), "--ref-dir", str(mock_tree / "ref"),
                     "--out", str(mock_tree / "tier")]) == 1
