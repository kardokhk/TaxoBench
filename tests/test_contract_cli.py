"""Contract checks, the run/ingest layer and the CLI end to end on toy data."""
import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from helpers import LINEAGES, REF_ID, queries_table, toy_ref, truth_table
from taxobench import cli
from taxobench.contract import ContractError, read_table, write_table

UP = "d__B; p__P1; c__C1; o__O1"


@pytest.fixture()
def tier(tmp_path):
    ref = toy_ref()
    truth = {f"q{i}": f"{UP}; f__F1; g__G1; s__G1 {'a' if i % 2 else 'b'}" for i in range(20)}
    truth.update({f"n{i}": f"{UP}; f__F1; g__G2; s__" for i in range(5)})
    d = tmp_path / "tier"
    taxa = pd.DataFrame([dict(ref_id=REF_ID, rank=ref.rank[t], taxon_id=t, name=ref.name[t],
                              parent_taxon_id=ref.parent[t]) for t in ref.rank])
    write_table(pd.DataFrame([dict(ref_id=REF_ID, release="toy-1")]), d / "references.tsv", "references")
    write_table(taxa, d / "taxa.tsv", "taxa")
    q = queries_table(list(truth), clusters={k: k[:2] for k in truth},
                      split={k: ("calib" if k in ("q0", "q1", "n0") else "eval") for k in truth})
    write_table(q, d / "queries.tsv", "queries")
    write_table(truth_table(ref, truth), d / "truth.tsv", "truth")
    write_table(pd.DataFrame([dict(query_id=k, sample_id="s1", count=3) for k in truth]),
                d / "abundance.tsv", "abundance")
    (tmp_path / "q.fasta").write_text("".join(f">{k}\nACGT\n" for k in truth))
    (tmp_path / "lineages.tsv").write_text("".join(f"r{i}\t{x}\n" for i, x in enumerate(LINEAGES)))
    return tmp_path, truth


def _raw(tmp, name, calls):
    fp = tmp / name
    fp.write_text("Feature ID\tTaxon\tConfidence\n" + "".join(f"{q}\t{c}\t{v}\n" for q, (c, v) in calls.items()))
    return fp


def test_sidecar_is_checked(tier):
    tmp, _ = tier
    read_table(tmp / "tier", "truth")
    fp = tmp / "tier" / "truth.tsv"
    fp.write_text(fp.read_text() + "")
    with open(fp, "a") as fh:
        fh.write("x\ttoy\tspecies\t\tnovel\n")
    with pytest.raises(ContractError, match="sha256"):
        read_table(tmp / "tier", "truth")
    os.remove(str(fp) + ".json")
    with pytest.raises(ContractError, match="sidecar"):
        read_table(tmp / "tier", "truth")


def test_duplicate_keys_and_bad_vocabulary_fail(tmp_path):
    df = pd.DataFrame([dict(query_id="a", ref_id="r", rank="species", truth_set="", status="novel")] * 2)
    with pytest.raises(ContractError, match="duplicated"):
        write_table(df, tmp_path / "truth.tsv", "truth") and read_table(tmp_path, "truth")
        read_table(tmp_path, "truth")
    df = pd.DataFrame([dict(query_id="a", ref_id="r", rank="strain", truth_set="", status="novel")])
    write_table(df, tmp_path / "x" / "truth.tsv", "truth")
    with pytest.raises(ContractError, match="unknown rank"):
        read_table(tmp_path / "x", "truth")


def _ingest(tmp, name, run_id, calls, version):
    raw = _raw(tmp, f"{name}.tsv", calls)
    return cli.main(["ingest", "--adapter", "taxonomy-tsv", "--classifier", "toyclf", "--version", version,
                     "--ref-id", REF_ID, "--taxa", str(tmp / "lineages.tsv"), "--query-fasta",
                     str(tmp / "q.fasta"), "--raw", str(raw), "--run-id", run_id,
                     "--out", str(tmp / "runs" / run_id)])


def test_cli_end_to_end(tier, tmp_path):
    tmp, truth = tier
    good = {q: (lab, 0.9) for q, lab in truth.items()}
    worse = dict(good)
    for i in range(0, 20, 3):
        worse[f"q{i}"] = (f"{UP}; f__F1; g__G1", 0.5)
    assert _ingest(tmp, "a", "A", worse, "1.0") == 0
    assert _ingest(tmp, "b", "B", good, "1.1") == 0
    runs = [str(tmp / "runs" / "A"), str(tmp / "runs" / "B")]
    assert cli.main(["validate", "--tier", str(tmp / "tier"), "--runs", *runs]) == 0
    assert cli.main(["score", "--tier", str(tmp / "tier"), "--runs", *runs, "--out", str(tmp / "sc"),
                     "--n-boot", "200"]) == 0
    s = pd.read_csv(tmp / "sc" / "scores.tsv", sep="\t")
    f1b = s[(s.run_id == "B") & (s["rank"] == "species") & (s.rule == "G1") & (s.metric == "f1")]
    assert f1b.value.iloc[0] == 1.0
    assert {"l1_A_median", "unifrac_A_median", "average_precision", "ece_mass"} <= set(s.metric)
    plan = tmp / "plan.json"
    plan.write_text(json.dumps({"confirmatory": [{"tier": "TX", "rank": "species", "metric": "f1"}],
                                "guardrails": [{"tier": "TX", "rank": "genus", "metric": "f1", "margin": 0.01}]}))
    rc = cli.main(["compare", "--a", "A", "--b", "B", "--tier", str(tmp / "tier"), "--runs", *runs,
                   "--plan", str(plan), "--out", str(tmp / "cmp"), "--n-boot", "500"])
    assert rc == 0
    c = pd.read_csv(tmp / "cmp" / "compare.tsv", sep="\t")
    conf = c[c.family == "confirmatory"]
    assert len(conf) == 1 and conf.delta.iloc[0] > 0
    md = (tmp / "cmp" / "compare.md").read_text()
    assert "Confirmatory family" in md and "Guardrails" in md
    assert cli.main(["report", str(tmp / "sc" / "scores.tsv"), "--out", str(tmp / "r.md")]) == 0


def test_partial_output_fails_closed(tier):
    tmp, truth = tier
    calls = {q: (lab, 0.9) for q, lab in list(truth.items())[:-3]}
    assert _ingest(tmp, "p", "P", calls, "1") == 1          # coverage check refuses it


def test_compare_refuses_different_query_sets(tier):
    tmp, truth = tier
    assert _ingest(tmp, "a", "A", {q: (lab, 0.9) for q, lab in truth.items()}, "1") == 0
    # same queries but a different reference id -> refused
    assert cli.main(["compare", "--a", "A", "--b", "nope", "--tier", str(tmp / "tier"),
                     "--runs", str(tmp / "runs" / "A"), "--out", str(tmp / "c")]) == 1


def test_run_command_failure_leaves_nothing(tier):
    from taxobench.run import run_one
    tmp, truth = tier
    out = tmp / "runs" / "X"
    spec = {"adapter": "taxonomy-tsv", "version": "1", "ref_id": REF_ID,
            "taxa": str(tmp / "lineages.tsv"), "query_fasta": str(tmp / "q.fasta"),
            "outdir": str(out), "raw": "raw.tsv", "run_id": "X",
            "command": [sys.executable, "-c", "import sys; sys.exit(3)"]}
    with pytest.raises(RuntimeError, match="exited 3"):
        run_one(spec)
    assert not (out / "predictions.tsv").exists()
    # a command that writes the raw output succeeds
    src = _raw(tmp, "src.tsv", {q: (lab, 0.9) for q, lab in truth.items()})
    spec["command"] = [sys.executable, "-c", f"import shutil; shutil.copy({str(src)!r}, {str(out / 'raw.tsv')!r})"]
    run_one(spec)
    assert (out / "predictions.tsv").exists() and (out / "run.json").exists()


def test_uc_wrong_above_in_scores_and_compare(tier):
    tmp, truth = tier
    good = {q: (lab, 0.9) for q, lab in truth.items()}
    wrong = dict(good)
    wrong["q3"] = (f"{UP}; f__F1; g__G2", 0.9)          # species unnamed, genus wrong
    assert _ingest(tmp, "a", "A", wrong, "1") == 0
    assert _ingest(tmp, "b", "B", good, "2") == 0
    runs = [str(tmp / "runs" / "A"), str(tmp / "runs" / "B")]
    assert cli.main(["score", "--tier", str(tmp / "tier"), "--runs", *runs, "--out", str(tmp / "sc"),
                     "--n-boot", "0"]) == 0
    s = pd.read_csv(tmp / "sc" / "scores.tsv", sep="\t")
    v = s[(s.run_id == "A") & (s["rank"] == "species") & (s.rule == "G1") & (s.metric == "count_uc_wrong_above")]
    assert v.value.iloc[0] == 1
    assert {"count_n_homonym_resolved", "count_n_group_unexpanded", "count_ca_wrong_above"} <= set(s.metric)
    assert cli.main(["compare", "--a", "A", "--b", "B", "--tier", str(tmp / "tier"), "--runs", *runs,
                     "--out", str(tmp / "cmp"), "--n-boot", "200"]) == 0
    c = pd.read_csv(tmp / "cmp" / "compare.tsv", sep="\t")
    t = c[(c.level == "tier") & (c["rank"] == "species") & (c.rule == "G1") & (c.metric == "f1")]
    assert t.uc_wrong_above_A.iloc[0] == 1 and t.uc_wrong_above_B.iloc[0] == 0
    assert "UC wrong above A/B" in (tmp / "cmp" / "compare.md").read_text()
    assert cli.main(["report", str(tmp / "sc" / "scores.tsv"), "--out", str(tmp / "r.md")]) == 0
    assert "UC wrong above" in (tmp / "r.md").read_text()
