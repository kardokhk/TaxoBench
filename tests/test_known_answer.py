"""Known-answer tests of the scoring core.

Hand-computed micro cell, 10 queries, covering every outcome:

| q   | truth            | call                         | species G1 | G2  | G3          |
|-----|------------------|------------------------------|------------|-----|-------------|
| q1  | G1 a             | G1 a                         | TP         | TP  | TP          |
| q2  | G1 a             | G1 b                         | MC         | MC  | MC          |
| q3  | G1 b             | genus G1 only                | UC         | UC  | UC          |
| q4  | G2, species none | G2 c                         | OC         | OC  | OC          |
| q5  | G2, species none | genus G2 only                | CA         | CA  | CA          |
| q6  | G1 a             | s__G1 zz (not in reference)  | MC (oov)   | MC  | MC          |
| q7  | G1 a             | [G1 a|G1 b]                  | TP         | UC  | 0.5 TP/0.5 MC |
| q8  | G2 c             | genus G1, [G1 a|G1 b]        | MC         | UC* | MC          |
| q9  | other            | G1 a                         | excluded   |     |             |
| q10 | G3 d             | genus G1, species G3 d       | MC (wrong above) | MC | MC     |

(*) uc_wrong_above. K = 7, L = 2.
"""
import math

import numpy as np
import pandas as pd
import pytest

from helpers import aligned, preds_table, queries_table, toy_ref, truth_table
from taxobench import metrics as M
from taxobench.contract import ContractError
from taxobench.metrics import align

UP = "d__B; p__P1; c__C1; o__O1"
TRUTH = {
    "q1": f"{UP}; f__F1; g__G1; s__G1 a", "q2": f"{UP}; f__F1; g__G1; s__G1 a",
    "q3": f"{UP}; f__F1; g__G1; s__G1 b", "q4": f"{UP}; f__F1; g__G2; s__",
    "q5": f"{UP}; f__F1; g__G2; s__", "q6": f"{UP}; f__F1; g__G1; s__G1 a",
    "q7": f"{UP}; f__F1; g__G1; s__G1 a", "q8": f"{UP}; f__F1; g__G2; s__G2 c",
    "q9": "other", "q10": f"{UP}; f__F2; g__G3; s__G3 d",
}
CALLS = {
    "q1": f"{UP}; f__F1; g__G1; s__G1 a", "q2": f"{UP}; f__F1; g__G1; s__G1 b",
    "q3": f"{UP}; f__F1; g__G1", "q4": f"{UP}; f__F1; g__G2; s__G2 c",
    "q5": f"{UP}; f__F1; g__G2", "q6": f"{UP}; f__F1; g__G1; s__G1 zz",
    "q7": f"{UP}; f__F1; g__G1; s__[G1 a|G1 b]", "q8": f"{UP}; f__F1; g__G1; s__[G1 a|G1 b]",
    "q9": f"{UP}; f__F1; g__G1; s__G1 a", "q10": f"{UP}; f__F1; g__G1; s__G3 d",
}


@pytest.fixture(scope="module")
def al():
    return aligned(toy_ref(), TRUTH, CALLS)


def test_outcomes_g1(al):
    df = M.outcomes(al, "species", "G1")
    got = dict(zip(df.query_id, df.outcome))
    assert got == {"q1": "TP", "q2": "MC", "q3": "UC", "q4": "OC", "q5": "CA", "q6": "MC",
                   "q7": "TP", "q8": "MC", "q9": "excluded", "q10": "MC"}
    c = M.totals(df)
    assert (c["K"], c["L"], c["TP"], c["MC"], c["UC"], c["OC"], c["CA"]) == (7, 2, 2, 4, 1, 1, 1)
    assert c["n_oov"] == 1 and c["n_group"] == 2
    r = M.rates(c)
    assert r["precision"] == pytest.approx(2 / 7)
    assert r["recall"] == pytest.approx(2 / 7)
    assert r["f1"] == pytest.approx(2 / 7)
    assert r["edgar_accuracy"] == pytest.approx(2 / 8)
    assert r["mcr"] == pytest.approx(4 / 7)
    assert r["ucr"] == pytest.approx(1 / 7)
    assert r["ocr"] == pytest.approx(1 / 2)
    assert r["coverage"] == pytest.approx(7 / 9)


def test_outcomes_g2_strict(al):
    df = M.outcomes(al, "species", "G2")
    got = dict(zip(df.query_id, df.outcome))
    assert got["q7"] == "UC" and got["q8"] == "UC"
    c = M.totals(df)
    assert (c["TP"], c["MC"], c["UC"], c["OC"], c["CA"], c["uc_wrong_above"]) == (1, 3, 3, 1, 1, 1)
    r = M.rates(c)
    assert r["precision"] == pytest.approx(1 / 5)
    assert r["recall"] == pytest.approx(1 / 7)
    assert r["f1"] == pytest.approx(1 / 6)
    assert r["coverage"] == pytest.approx(5 / 9)


def test_outcomes_g3_size_weighted(al):
    c = M.totals(M.outcomes(al, "species", "G3"))
    assert c["TP"] == pytest.approx(1.5) and c["MC"] == pytest.approx(4.5)
    assert M.rates(c)["precision"] == pytest.approx(1.5 / 7)


def test_genus_rank(al):
    c = M.totals(M.outcomes(al, "genus", "G1"))
    assert (c["K"], c["TP"], c["MC"], c["UC"], c["L"]) == (9, 7, 2, 0, 0)


def test_f1_from_unrounded_counts():
    # R8: P = 0.9532.., R = 0.9404..; F must come from the counts, not rounded P and R
    c = dict.fromkeys(M.COUNT_COLS, 0.0)
    c.update(K=1000, TP=940.4, MC=46.2, L=0, OC=0)
    r = M.rates(c)
    p, rr = 940.4 / 986.6, 940.4 / 1000
    assert r["f1"] == pytest.approx(2 * p * rr / (p + rr), abs=1e-12)


def test_all_abstain_and_perfect():
    ref = toy_ref()
    truth = {k: v for k, v in TRUTH.items() if k != "q9"}
    a = aligned(ref, truth, {q: "" for q in truth})
    r = M.rates(M.totals(M.outcomes(a, "species")))
    assert r["recall"] == 0 and math.isnan(r["precision"])        # R9
    perfect = {q: lab for q, lab in truth.items()}
    a = aligned(ref, truth, perfect)
    r = M.rates(M.totals(M.outcomes(a, "species")))
    assert r["precision"] == 1 and r["recall"] == 1 and r["f1"] == 1
    assert r["ocr"] == 0


def test_invariant_row_permutation(al):
    ref = toy_ref()
    keys = list(TRUTH)[::-1]
    a2 = aligned(ref, {k: TRUTH[k] for k in keys}, {k: CALLS[k] for k in keys})
    for rule in ("G1", "G2", "G3"):
        assert M.rates(M.totals(M.outcomes(al, "species", rule))) == \
            pytest.approx(M.rates(M.totals(M.outcomes(a2, "species", rule))), nan_ok=True)


def test_invariant_consistent_renaming():
    ren = lambda s: (s.replace("G1", "Hx").replace("G2", "Hy").replace("F1", "Fz"))
    lins = [ren(x) for x in __import__("helpers").LINEAGES]
    ref2 = toy_ref(lins)
    a2 = aligned(ref2, {k: ren(v) for k, v in TRUTH.items()}, {k: ren(v) for k, v in CALLS.items()})
    a1 = aligned(toy_ref(), TRUTH, CALLS)
    for rule in ("G1", "G2", "G3"):
        assert M.rates(M.totals(M.outcomes(a1, "species", rule))) == \
            pytest.approx(M.rates(M.totals(M.outcomes(a2, "species", rule))), nan_ok=True)


# --- vocabulary (item 4) -------------------------------------------------------

def test_wrong_name_in_reference_is_mc_and_absent_is_oov():
    ref = toy_ref()
    p = preds_table(ref, {"x": f"{UP}; f__F1; g__G1; s__G1 b", "y": f"{UP}; f__F1; g__G1; s__G1 zz"})
    s = p[p["rank"] == "species"].set_index("query_id")
    assert s.loc["x", "call"] == "named" and s.loc["y", "call"] == "oov"


def test_unknown_rank_prefix_fails():
    # R7: k__Bacteria was a silent domain misclassification in v1/v2
    from taxobench.taxonomy import parse_lineage
    with pytest.raises(ContractError):
        parse_lineage("k__Bacteria; p__P1")


def test_homonym_resolved_by_parent():
    # R5: identity is the reference path, not the bare name
    ref = toy_ref(["d__B; p__P1; c__C1; o__O1; f__F1; g__Same",
                   "d__B; p__P1; c__C1; o__O1; f__F2; g__Same"])
    res = ref.resolve_lineage({"domain": "B", "phylum": "P1", "class": "C1", "order": "O1",
                               "family": "F2", "genus": "Same"})
    assert res["genus"].call == "named" and "f__F2" in res["genus"].taxon_ids[0]
    res = ref.resolve_lineage({"genus": "Same"})
    assert res["genus"].call == "oov" and "ambiguous" in res["genus"].note


# --- group handling (item 5) ----------------------------------------------------

def test_group_parser_cases():
    from taxobench.taxonomy import split_group
    assert split_group("[[Ruminococcus] gnavus group|Blautia]") == ["[Ruminococcus] gnavus group", "Blautia"]
    assert split_group("[Ruminococcus] gnavus group") == ["[Ruminococcus] gnavus group"]
    assert split_group("[A|B|C]") == ["A", "B", "C"]
    assert split_group("[A]") == ["A"]                      # R10: one member
    assert split_group("Plain name") == ["Plain name"]


def test_silva_bracket_name_is_one_name():
    ref = toy_ref(["d__B; p__P1; c__C1; o__O1; f__Lachno; g__[Ruminococcus] gnavus group",
                   "d__B; p__P1; c__C1; o__O1; f__Lachno; g__Blautia"])
    res = ref.resolve_lineage({"genus": "[Ruminococcus] gnavus group"})
    assert res["genus"].call == "named"
    res = ref.resolve_lineage({"genus": "[[Ruminococcus] gnavus group|Blautia]"})
    assert res["genus"].call == "group" and len(res["genus"].taxon_ids) == 2


def test_cross_genus_group_g1_correct_g2_nocall():
    ref = toy_ref()
    truth = {"q": f"{UP}; f__F1; g__G2; s__G2 c"}
    a = aligned(ref, truth, {"q": f"{UP}; f__F1; g__[G1|G2]; s__[G1 a|G2 c]"})
    assert M.outcomes(a, "species", "G1").outcome[0] == "TP"
    assert M.outcomes(a, "species", "G2").outcome[0] == "UC"


def test_unresolvable_truth_matches_any_member():
    ref = toy_ref()
    t = truth_table(ref, {"q": f"{UP}; f__F1; g__G1; s__[G1 a|G1 b]"})
    assert t[t["rank"] == "species"].status.iloc[0] == "unresolvable"
    a = align(t, queries_table(["q"]), preds_table(ref, {"q": f"{UP}; f__F1; g__G1; s__G1 b"}), ref, ["q"])
    assert M.outcomes(a, "species").outcome[0] == "TP"


# --- truth status -------------------------------------------------

def test_genus_only_truth_is_novel_not_tp():
    # R3: an identical genus-stopping call is a correct abstention, never a species TP
    ref = toy_ref()
    a = aligned(ref, {"q": f"{UP}; f__F1; g__G1; s__"}, {"q": f"{UP}; f__F1; g__G1"})
    df = M.outcomes(a, "species")
    assert df.outcome[0] == "CA" and df.TP.sum() == 0
    assert M.rates(M.totals(df))["coverage"] == 0          # D5: abstention visible


def test_unscorable_truth_never_matches():
    # R1/R2/D11: "other" or no named rank: excluded, never a match
    ref = toy_ref()
    a = aligned(ref, {"q": "other"}, {"q": f"{UP}; f__F1; g__G1; s__G1 a"})
    assert M.outcomes(a, "species").outcome[0] == "excluded"


# --- profiles (item 6) ----------------------------------------------------------

def _ab(rows):
    return pd.DataFrame(rows, columns=["query_id", "sample_id", "count"])


def test_l1_is_twice_bray_curtis_convention_a(al):
    ab = _ab([(q, "s1", n) for q, n in zip(TRUTH, [5, 3, 2, 4, 1, 7, 2, 3, 1, 6])])
    pr = M.profiles(al, toy_ref(), ab, "species")
    assert pr.l1_A.iloc[0] == pytest.approx(2 * pr.bray_curtis_A.iloc[0])


def test_unifrac_three_leaf_hand_value():
    # tree: root - B - P1 - ... - F1 - G1 - {G1 a, G1 b}; all mass on G1 a vs G1 b:
    # the paths differ on two unit edges (to each leaf), so the distance is 2.
    ref = toy_ref()
    a = "d__B;p__P1;c__C1;o__O1;f__F1;g__G1;s__G1 a"
    b = "d__B;p__P1;c__C1;o__O1;f__F1;g__G1;s__G1 b"
    assert M.weighted_unifrac(ref, {a: 1.0}, {b: 1.0}) == pytest.approx(2.0)
    assert M.weighted_unifrac(ref, {a: 1.0}, {a: 1.0}) == 0
    # half the mass at the genus node (a group's LCA): one edge differs for that half
    g = "d__B;p__P1;c__C1;o__O1;f__F1;g__G1"
    assert M.weighted_unifrac(ref, {a: 1.0}, {a: 0.5, g: 0.5}) == pytest.approx(0.5)


def test_group_mass_at_lca_independent_of_truth():
    # R12: the member choice never looks at the truth
    ref = toy_ref()
    call = {"q": f"{UP}; f__F1; g__G1; s__[G1 a|G1 b]"}
    ab = _ab([("q", "s", 10)])
    r1 = M.profiles(aligned(ref, {"q": f"{UP}; f__F1; g__G1; s__G1 a"}, call), ref, ab, "species")
    r2 = M.profiles(aligned(ref, {"q": f"{UP}; f__F1; g__G1; s__G1 b"}, call), ref, ab, "species")
    assert r1.unifrac_A.iloc[0] == pytest.approx(r2.unifrac_A.iloc[0]) == pytest.approx(1.0)


def test_perfect_profile_zero_and_duplicate_sample_changes_nothing():
    ref = toy_ref()
    truth = {k: v for k, v in TRUTH.items() if k != "q9"}
    a = aligned(ref, truth, dict(truth))
    ab = _ab([(q, "s1", i + 1) for i, q in enumerate(truth)])
    pr = M.profiles(a, ref, ab, "species")
    assert pr.l1_A.iloc[0] == 0 and pr.unifrac_A.iloc[0] == 0 and pr.l1_B.iloc[0] == 0
    ab2 = pd.concat([ab, ab.assign(sample_id="s2")])
    pr2 = M.profiles(aligned(ref, TRUTH, CALLS), ref, ab2, "species")
    cols = [c for c in pr2.columns if c != "sample_id"]
    assert np.allclose(pr2[cols].iloc[0].to_numpy(float), pr2[cols].iloc[1].to_numpy(float), equal_nan=True)


def test_unifrac_depends_on_the_reference_passed():
    # R13: no cache keyed only by the tree; a different taxonomy gives a different value
    a = "d__B;p__P1;c__C1;o__O1;f__F1;g__G1;s__G1 a"
    b = "d__B;p__P1;c__C1;o__O1;f__F1;g__G1;s__G1 b"
    deep = toy_ref()
    assert M.weighted_unifrac(deep, {a: 1}, {b: 1}) == pytest.approx(2.0)


# --- M8 and M9 -------------------------------------------------------------------

def test_curve_ap_and_threshold():
    df = pd.DataFrame({"outcome": ["TP", "TP", "MC", "TP", "UC"], "K": [1] * 5, "L": [0] * 5,
                       "TP": [1, 1, 0, 1, 0], "MC": [0, 0, 1, 0, 0], "OC": [0] * 5, "UC": [0, 0, 0, 0, 1],
                       "CA": [0] * 5, "confidence": [0.9, 0.8, 0.7, 0.6, np.nan],
                       **{c: [0] * 5 for c in ["uc_wrong_above", "n_oov", "n_group", "group_size_sum"]}})
    cv = M.curve(df)
    assert list(cv.precision.round(4)) == [1.0, 1.0, round(2 / 3, 4), 0.75]
    # AP = 0.2*1 + 0.2*1 + 0*(2/3) + 0.2*0.75
    assert M.average_precision(cv) == pytest.approx(0.55)
    assert M.fit_threshold(df) == pytest.approx(0.6)


def test_calibration_perfect_and_off():
    n = 100
    df = pd.DataFrame({"outcome": ["TP"] * n, "confidence": [1.0] * n})
    cal = M.calibration(df)
    assert cal["ece_mass"] == 0 and cal["brier"] == 0
    df = pd.DataFrame({"outcome": ["TP"] * 50 + ["MC"] * 50, "confidence": [0.9] * n})
    cal = M.calibration(df)
    assert cal["ece_mass"] == pytest.approx(0.4) and cal["ece_width"] == pytest.approx(0.4)


# --- homonyms and group members ---------------------------

HOMO = [
    "d__B; p__P1; c__C1; o__O1; f__Fa; g__Pantoea; s__Pantoea x",
    "d__B; p__P1; c__C1; o__O1; f__[Fa|Fb]; g__Pantoea; s__Pantoea y",
    "d__B; p__P1; c__C1; o__O1; f__Fb; g__Erwinia; s__Erwinia z",
    "d__B; p__P1; c__C1; o__O1; f__Fc; g__[Pantoea|Erwinia]",
    "d__B; p__P1; c__C1; o__O1; f__Fc; g__[Pantoea|Ghost]",
]


def test_group_member_on_homonymous_paths_resolved_by_name():
    # decision 2: g__Pantoea exists under f__Fa and under f__[Fa|Fb]; the group's
    # parent f__Fc narrows to neither, so the member stands for both paths
    ref = toy_ref(HOMO)
    gid = ref.index[("genus", "[Pantoea|Erwinia]")][0]
    mem = ref.group_members(gid)
    assert len(mem) == 2 and sorted(len(m) for m in mem) == [1, 2]
    assert gid in ref.homonym_groups
    truth = {"q": "d__B; p__P1; c__C1; o__O1; f__Fa; g__Pantoea; s__Pantoea x"}
    a = aligned(ref, truth, {"q": "d__B; p__P1; c__C1; o__O1; g__[Pantoea|Erwinia]"})
    df = M.outcomes(a, "genus", "G1")
    assert df.outcome[0] == "TP" and df.n_homonym_resolved[0] == 1
    assert df.group_size_sum[0] == 2                    # two members, not three paths
    assert M.outcomes(a, "genus", "G3").TP[0] == pytest.approx(0.5)


def test_group_with_absent_member_stays_unexpanded_and_is_counted():
    ref = toy_ref(HOMO)
    gid = ref.index[("genus", "[Pantoea|Ghost]")][0]
    assert ref.group_members(gid) == (frozenset([gid]),) and gid in ref.unexpanded_groups
    truth = {"q": "d__B; p__P1; c__C1; o__O1; f__Fa; g__Pantoea; s__Pantoea x"}
    a = aligned(ref, truth, {"q": "d__B; p__P1; c__C1; o__O1; g__[Pantoea|Ghost]"})
    df = M.outcomes(a, "genus", "G1")
    assert df.outcome[0] == "MC" and df.n_group_unexpanded[0] == 1 and df.n_homonym_resolved[0] == 0


def test_unnamed_but_wrong_above_is_uc_and_ca():
    # decision 1: not named at r, wrong higher up -> UC (K) or CA (L), with sub-counts;
    # precision at r judges only calls named at r
    ref = toy_ref()
    truth = {"k": f"{UP}; f__F1; g__G2; s__G2 c", "l": f"{UP}; f__F1; g__G2; s__"}
    a = aligned(ref, truth, {"k": f"{UP}; f__F1; g__G1", "l": f"{UP}; f__F1; g__G1"})
    df = M.outcomes(a, "species", "G1").set_index("query_id")
    assert df.loc["k", "outcome"] == "UC" and df.loc["k", "uc_wrong_above"] == 1
    assert df.loc["l", "outcome"] == "CA" and df.loc["l", "ca_wrong_above"] == 1
    r = M.rates(M.totals(df))
    assert r["mcr"] == 0 and r["ucr"] == 1 and r["ocr"] == 0


def test_placeholder_member_taxon_resolved_by_name():
    # a tier may write bracketed members as taxa with in_reference = 0 under the group's
    # parent; a member that is also a reference taxon elsewhere resolves by name
    import pandas as pd
    from taxobench.taxonomy import Reference
    rows = [("r", "domain", "d", "B", "", 0, 1), ("r", "family", "fA", "Fa", "d", 0, 1),
            ("r", "family", "fB", "Fb", "d", 0, 1), ("r", "family", "fG", "[Fa|Fb]", "d", 1, 1),
            ("r", "genus", "gP", "Pantoea", "fA", 0, 1), ("r", "genus", "gE", "Erwinia", "fB", 0, 1),
            ("r", "genus", "gPx", "Pantoea", "fG", 0, 0), ("r", "genus", "gEx", "Erwinia", "fG", 0, 0),
            ("r", "genus", "gG", "[Pantoea|Erwinia]", "fG", 1, 1)]
    taxa = pd.DataFrame(rows, columns=["ref_id", "rank", "taxon_id", "name", "parent_taxon_id",
                                       "is_group", "in_reference"]).astype(str)
    ref = Reference.from_taxa(taxa, "r")
    assert ref.expand(["gG"]) == frozenset(["gP", "gPx", "gE", "gEx"])
    assert "gG" in ref.homonym_groups
    assert ref.expand(["gG"]) & ref.expand(["gP"])          # truth Pantoea is a member
    assert len(ref.members(["gG"])) == 2                     # still two members
    assert ref.expand(["gPx"]) == frozenset(["gP", "gPx"])   # a placeholder named alone


def _placeholder_ref(rows):
    import pandas as pd
    from taxobench.taxonomy import Reference
    taxa = pd.DataFrame(rows, columns=["ref_id", "rank", "taxon_id", "name", "parent_taxon_id",
                                       "is_group", "in_reference"]).astype(str)
    return Reference.from_taxa(taxa, "r")


def test_guard_allows_one_lineage_through_parent_group():
    # Pantoea under Erwiniaceae and under [Erwiniaceae|Yersiniaceae] is one lineage: allowed
    ref = _placeholder_ref([
        ("r", "domain", "d", "B", "", 0, 1),
        ("r", "family", "fE", "Erwiniaceae", "d", 0, 1), ("r", "family", "fY", "Yersiniaceae", "d", 0, 1),
        ("r", "family", "fG", "[Erwiniaceae|Yersiniaceae]", "d", 1, 1),
        ("r", "genus", "gP", "Pantoea", "fE", 0, 1), ("r", "genus", "gS", "Serratia", "fY", 0, 1),
        ("r", "genus", "gPx", "Pantoea", "fG", 0, 0), ("r", "genus", "gSx", "Serratia", "fG", 0, 0),
        ("r", "genus", "gG", "[Pantoea|Serratia]", "fG", 1, 1)])
    assert ref.expand(["gG"]) == frozenset(["gP", "gPx", "gS", "gSx"])
    assert "gG" in ref.homonym_groups and "gG" not in ref.unexpanded_groups


def test_guard_blocks_distinct_lineages_and_non_specific_names():
    # endosymbionts under two unrelated families: distinct lineages, group stays unexpanded
    ref = _placeholder_ref([
        ("r", "domain", "d", "B", "", 0, 1),
        ("r", "family", "f1", "Rickettsiaceae", "d", 0, 1), ("r", "family", "f2", "Enterobacteriaceae", "d", 0, 1),
        ("r", "family", "f3", "Anaplasmataceae", "d", 0, 1),
        ("r", "family", "fG", "[Rickettsiaceae|Anaplasmataceae]", "d", 1, 1),
        ("r", "genus", "g1", "endosymbionts", "f1", 0, 1), ("r", "genus", "g2", "endosymbionts", "f2", 0, 1),
        ("r", "genus", "gW", "Wolbachia", "f3", 0, 1),
        ("r", "genus", "gEx", "endosymbionts", "fG", 0, 0), ("r", "genus", "gWx", "Wolbachia", "fG", 0, 0),
        ("r", "genus", "gG", "[endosymbionts|Wolbachia]", "fG", 1, 1),
        # a specific name whose candidates sit under unrelated families is guarded too
        ("r", "genus", "h1", "Same", "f1", 0, 1), ("r", "genus", "h2", "Same", "f2", 0, 1),
        ("r", "genus", "hx", "Same", "fG", 0, 0)])
    assert ref.group_members("gG") == (frozenset(["gG"]),) and "gG" in ref.unexpanded_groups
    assert ref.expand(["hx"]) == frozenset(["hx"]) and "hx" in ref.unexpanded_groups
    from taxobench.taxonomy import is_non_specific
    assert is_non_specific("uncultured") and is_non_specific("endosymbionts")
    assert is_non_specific("metagenome") and not is_non_specific("Pantoea")
