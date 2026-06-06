#!/usr/bin/env python

import os
import glob
import csv
import shutil
import re
import skbio
import qiime2
import pandas
import seaborn as sns
import matplotlib
import numpy
import pandas as pd
from skbio import DNA
from qiime2 import Artifact

def single_character_matches(char, query, subject):
    query = ((query == char) + (query == 'N')).astype(int)
    subject = ((subject == char) + (subject == 'N')).astype(int)
    return numpy.convolve(query[::-1], subject)


def gapless_mismatches(query, subject):
    matches = single_character_matches('A', query, subject)
    matches += single_character_matches('C', query, subject)
    matches += single_character_matches('G', query, subject)
    matches += single_character_matches('T', query, subject)
    return len(query) - matches.max()


def compare_observed_and_expected(observed, expected):
    ids = [obs.metadata['id'] for obs in observed]
    rev_observed = [numpy.array(list(str(obs.reverse_complement())))
                    for obs in observed]
    observed = [numpy.array(list(str(obs))) for obs in observed]
    keys = [' '.join([exp.metadata['id'], exp.metadata['description']]).strip()
            for exp in expected]
    expected = [numpy.array(list(str(exp))) for exp in expected]
    comparison = {}
    for _id, obs, rev in zip(ids, observed, rev_observed):
        matches = [(gapless_mismatches(obs, exp), key) for
                   exp, key in zip(expected, keys)]
        matches += [(gapless_mismatches(rev, exp), key) for
                    exp, key in zip(expected, keys)]
        matches.sort()
        comparison[_id] = [(e, k) for e, k in matches if e == matches[0][0]]
    return comparison



def generate_trueish_taxonomies(data_dir, feature_table,
                                expected_sequences, expected_taxonomy,
                                observed_sequences, trueish_taxonomy,
                                max_mismatches, force=False):
    # load the taxonomy map
    with open(expected_taxonomy) as tf_fh:
        reader = csv.reader(tf_fh, delimiter='\t')
        next(reader)
        tax_map = {seq: tax for seq, tax in reader}

    # load the expected sequences
    expected = skbio.io.read(
        expected_sequences, format='fasta', constructor=skbio.DNA)
    expected = list(expected)

    # compare with the observed sequences
    observed = skbio.io.read(
        observed_sequences, format='fasta', constructor=skbio.DNA)
    observed = list(observed)
    mismatches = compare_observed_and_expected(observed, expected)

    # save the taxonomies of the matched reads
    with open(trueish_taxonomy, 'w') as result_fh:
        writer = csv.writer(result_fh, delimiter='\t')
        writer.writerow(['id', 'taxonomy'])
        for rep_id, mismatches_for_id in mismatches.items():
            if mismatches_for_id[0][0] > max_mismatches:
                writer.writerow([rep_id, 'other'])
                continue
            matched_taxonomies = set()
            for m, e in mismatches_for_id:
                if e in tax_map:
                    matched_taxonomies.add(tax_map[e])
                else:
                    matched_taxonomies.add('NO-EXPECTED-TAXONOMY')
                    print('WARNING: no expected taxonomy for ' + e)
            matched_taxonomies = list(matched_taxonomies)
            if len(matched_taxonomies) > 1:
                print('WARNING: ' + rep_id + ' matches')
                print(' and\n'.join(matched_taxonomies))
                print('with ' + str(mismatches_for_id[0][0]) + ' mismatches')
            for matched_taxonomy in matched_taxonomies:
                writer.writerow([rep_id, matched_taxonomy])

    features = Artifact.load(feature_table)
    features = features.view(pd.DataFrame)
    best_mismatches = list(mm[0][0] for mm in mismatches.values())
    weighted_mismatches = []
    found = 0
    total = 0
    for obs in mismatches:
        try:
            count = int(features[obs])
        except TypeError:
            count = int(sum(features[obs]))
        num_mismatches = mismatches[obs][0][0]
        weighted_mismatches.extend([num_mismatches]*count)
        if num_mismatches <= max_mismatches:
            found += count
        total += count
    print('Guessed taxonomy for %d of %d reads (%.1f%%)' %
        (found, total, 100*found/total))
    sns.distplot(best_mismatches, kde=False, bins=50,
        axlabel='Number of Mismatches')
    matplotlib.pyplot.show()
    sns.distplot(weighted_mismatches, kde=False, bins=50,
        axlabel='Number of Mismatches')
    matplotlib.pyplot.show()
