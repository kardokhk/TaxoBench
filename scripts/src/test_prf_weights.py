#!/usr/bin/env python3
from evaluate_results import precision_recall_fscore

def test_weighted_vs_unweighted():
    exp = ['A','B','B','C']
    obs = ['A','B','B_wrong','Unclassified']
    weights = [1,2,3,4]
    _, _, p_w, r_w, f_w, _ = precision_recall_fscore(exp, obs, sample_weight=weights)
    _, _, p_u, r_u, f_u, _ = precision_recall_fscore(exp, obs)
    assert (p_w, r_w, round(f_w,3)) == (0.5, 0.3, 0.375)
    assert (p_u, r_u, round(f_u,3)) == (0.667, 0.5, 0.572)
