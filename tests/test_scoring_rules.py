"""Known-answer tests for the v2 scoring rules (abstention, group labels, ties
from other tools) and for the v1 behaviour they replace."""
import unittest

import _paths  # noqa: F401
import scoring_v2 as S
import evaluate_results as ER

ev = _paths.load_evaluator()

LIN = "d__Bacteria; p__P1; c__C1; o__O1; f__F1; g__G1"
A = LIN + "; s__G1 a"
B = LIN + "; s__G1 b"
C = LIN + "; s__G1 c"
GEN_ONLY = LIN + "; s__"                      # Taxonova / kraken2 / mapseq abstention
UNC = LIN + "; s__unclassified_g__G1"         # IDTAXA abstention
NB_GENUS = LIN                                # naive Bayes / BLAST / VSEARCH truncation
SPINGO_A = "d__; p__; c__; o__; f__; g__G1; s__G1 a"
SPINGO_AMB = "d__; p__; c__; o__; f__; g__G1; s__"
SP = 6


def pairs(*rows):
    return {f"id{i}": r for i, r in enumerate(rows)}


class TestAbstention(unittest.TestCase):
    def test_each_abstention_format_is_fn_not_fp(self):
        for call in (GEN_ONLY, UNC, NB_GENUS, SPINGO_AMB, "Unassigned"):
            res = S.score_rows(pairs((A, A), (B, call)), SP, "v2", universe_size=100)
            self.assertEqual((res["tp"], res["fp"], res["fn"]), (1, 0, 1), call)
            self.assertEqual(res["precision"], 1.0)
            self.assertEqual(res["recall"], 0.5)
            self.assertEqual(res["fp_taxa"], 0, call)          # never a FP taxon
            self.assertEqual(res["unassigned_fraction"], 0.5)
            # universe 100, union {A, B} -> TN 98, FP 0
            self.assertEqual(res["specificity"], 1.0)
            # Jaccard over named: E={A,B}, O={A}
            self.assertEqual(res["jaccard"], 0.5)
            # L1 over named calls: exp A .5 B .5 ; obs A 1.0 -> |.5-1| + |.5-0| = 1.0
            self.assertEqual(res["l1_error"], 1.0)

    def test_v1_counts_placeholder_as_fp_taxon(self):
        ev._UNIVERSE_CACHE.clear()
        ev._UNIVERSE_CACHE.update({("Rx", 6): 100})
        row = ev._v1_level_row(pairs((A, A), (B, GEN_ONLY)), "Rx", 6,
                               ER.precision_recall_fscore, drop_abstentions=False)
        self.assertEqual(row["fp_taxa"], 1)                   # the "...; g__G1; s__" taxon
        self.assertEqual(row["specificity"], round(97 / 98, 6))
        self.assertEqual((row["precision"], row["recall"]), (1.0, 0.5))  # P/R/F unaffected

    def test_abstention_at_genus_rank_also_no_call_at_species(self):
        call = "d__Bacteria; p__P1; c__C1; o__O1; f__F1; g__; s__"
        res = S.score_rows(pairs((A, call)), 5, "v2", universe_size=50)
        self.assertEqual((res["tp"], res["fp"], res["fn"]), (0, 0, 1))
        self.assertEqual(res["fp_taxa"], 0)

    def test_wrong_higher_rank_is_error_in_prf_but_not_species_taxon(self):
        call = "d__Bacteria; p__P1; c__C1; o__O1; f__F2"
        res = S.score_rows(pairs((A, call)), SP, "v2", universe_size=50)
        self.assertEqual((res["tp"], res["fp"], res["fn"]), (0, 1, 1))
        self.assertEqual(res["fp_taxa"], 0)
        res_f = S.score_rows(pairs((A, call)), 4, "v2", universe_size=50)
        self.assertEqual(res_f["fp_taxa"], 1)

    def test_l1_l2_over_named_and_unassigned_fraction(self):
        p = pairs((A, A), (A, A), (B, B), (B, GEN_ONLY))
        res = S.score_rows(p, SP, "v2", universe_size=100)
        # obs over 3 named calls: A 2/3, B 1/3 ; exp A 1/2, B 1/2
        self.assertAlmostEqual(res["l1_error"], round(1 / 6 + 1 / 6, 6))
        self.assertAlmostEqual(res["l2_error"], round(2 * (1 / 6) ** 2, 6))
        self.assertEqual(res["unassigned_fraction"], 0.25)


class TestGroupLabels(unittest.TestCase):
    G_IN = LIN + "; s__[G1 a|G1 b]"
    G_OUT = LIN + "; s__[G1 b|G1 c]"

    def test_truth_in_group_is_match(self):
        res = S.score_rows(pairs((A, self.G_IN)), SP, "v2", universe_size=10)
        self.assertEqual((res["tp"], res["fp"], res["fn"]), (1, 0, 0))
        self.assertEqual(res["fp_taxa"], 0)

    def test_truth_not_in_group_is_wrong(self):
        res = S.score_rows(pairs((A, self.G_OUT)), SP, "v2", universe_size=10)
        self.assertEqual((res["tp"], res["fp"], res["fn"]), (0, 1, 1))
        self.assertEqual(res["fp_taxa"], 1)                  # one opaque group taxon

    def test_v1_same_two_cases(self):
        self.assertEqual(ER.evaluate_classification(self.G_IN, A), "match")
        self.assertEqual(ER.evaluate_classification(self.G_OUT, A), "misclassification")

    def test_strict_group_is_no_call(self):
        res = S.score_rows(pairs((A, self.G_IN), (B, B)), SP, "v2-strict", universe_size=10)
        self.assertEqual((res["tp"], res["fp"], res["fn"]), (1, 0, 1))
        self.assertEqual(res["fp_taxa"], 0)
        self.assertEqual(res["unassigned_fraction"], 0.5)

    def test_cross_genus_group(self):
        # Taxonova vaginal example shape: members from two orders.
        g = ("d__Bacteria; p__P1; c__C1; o__[O1|O2]; f__[F1|F2]; g__[G1|G2]; "
             "s__[G1 a|G2 z]")
        self.assertEqual(S.classify(S.parse_lineage(g), S.parse_lineage(A), SP), "match")
        z_truth = "d__Bacteria; p__P1; c__C1; o__O2; f__F2; g__G2; s__G2 z"
        self.assertEqual(S.classify(S.parse_lineage(g), S.parse_lineage(z_truth), SP), "match")
        y_truth = "d__Bacteria; p__P1; c__C1; o__O2; f__F2; g__G2; s__G2 y"
        self.assertEqual(S.classify(S.parse_lineage(g), S.parse_lineage(y_truth), SP),
                         "misclassification")
        # strict: no call from order down -> under-classification at species
        self.assertEqual(S.classify(S.apply_strict(S.parse_lineage(g)), S.parse_lineage(A), SP),
                         "underclassification")
        # at class rank the group row is a plain correct call in both
        self.assertEqual(S.classify(S.apply_strict(S.parse_lineage(g)), S.parse_lineage(A), 2),
                         "match")

    def test_group_with_wrong_genus_is_wrong_even_if_species_member(self):
        g = "d__Bacteria; p__P1; c__C1; o__O1; f__F1; g__G9; s__[G1 a|G9 q]"
        self.assertEqual(S.classify(S.parse_lineage(g), S.parse_lineage(A), SP),
                         "misclassification")

    def test_real_bracketed_name_is_not_a_group(self):
        t = LIN + "; s__[Ruminococcus] gnavus"
        self.assertEqual(S.parse_lineage(t)[SP], frozenset(["[Ruminococcus] gnavus"]))

    def test_group_stats(self):
        p = pairs((A, self.G_IN), (A, self.G_OUT), (B, B), (C, GEN_ONLY),
                  (A, LIN + "; s__[G1 a|G1 b|G1 c|G1 d]"))
        g = S.group_stats(p)
        self.assertEqual(g["n_named_calls"], 4)
        self.assertEqual(g["n_group_calls"], 3)
        self.assertEqual(g["group_share"], 0.75)
        self.assertEqual(g["group_size_mean"], round((2 + 2 + 4) / 3, 3))
        self.assertEqual(g["group_size_max"], 4)
        self.assertEqual((g["group_truth_in"], g["group_truth_out"]), (2, 1))


class TestOtherToolFormats(unittest.TestCase):
    def test_spingo_blank_upper_ranks(self):
        # Correct SPINGO call: a match in P/R/F (v1 and v2) ...
        self.assertEqual(ER.evaluate_classification(SPINGO_A, A), "match")
        res = S.score_rows(pairs((A, SPINGO_A), (B, B)), SP, "v2", universe_size=100)
        self.assertEqual((res["tp"], res["fp"], res["fn"]), (2, 0, 0))
        # ... and the true taxon in the v2 taxon set
        self.assertEqual(res["fp_taxa"], 0)
        self.assertEqual(res["detection_rate"], 1.0)
        # v1 taxon set: the blank-rank string is a different taxon (FP + miss)
        ev._UNIVERSE_CACHE.clear()
        ev._UNIVERSE_CACHE.update({("Rx", 6): 100})
        row = ev._v1_level_row(pairs((A, SPINGO_A), (B, B)), "Rx", 6,
                               ER.precision_recall_fscore, drop_abstentions=False)
        self.assertEqual(row["fp_taxa"], 1)
        self.assertEqual(row["detection_rate"], 0.5)

    def test_blast_tie_reconstruction(self):
        hits = [A, B, B]
        lab = S.tie_group_label(NB_GENUS, hits)
        self.assertEqual(lab, LIN + "; s__[G1 a|G1 b]")
        self.assertEqual(S.classify(S.parse_lineage(lab), S.parse_lineage(A), SP), "match")
        self.assertEqual(S.classify(S.parse_lineage(lab), S.parse_lineage(C), SP),
                         "misclassification")
        # a single-species call is left alone; so is a tie within one species
        self.assertEqual(S.tie_group_label(A, hits), A)
        self.assertEqual(S.tie_group_label(NB_GENUS, [A, A]), NB_GENUS)

    def test_rule_is_identical_across_classifiers(self):
        # same label string -> same outcome whatever the method name: score_rows
        # has no method argument at all
        import inspect
        self.assertNotIn("method", inspect.signature(S.score_rows).parameters)


if __name__ == "__main__":
    unittest.main()
