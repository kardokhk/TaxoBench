"""Taxobench scoring rules, v2 (2026-10-29), with the published v1 rules kept.

Every rule here is applied to every classifier's ``taxonomy_final.tsv`` in the
same way; no rule looks at which classifier produced a label.

Rule sets (``RULES``):

``v1``
    The published evaluator. Per-row P/R/F from
    ``evaluate_results.evaluate_classification`` (group labels credited rank by
    rank). Taxon-set metrics and L1/L2 over the truncated label strings, so an
    abstention (``...; g__X; s__``) is a taxon in its own right (false positive).
    Implemented by the untouched v1 code path in ``04_evaluate.py``.

``v2``
    (a) Abstention. A rank with no name (empty, ``s__``, ``unclassified_*``,
        ``Unassigned``, ``Unclassified``, ``AMBIGUOUS``, ``No blast hit``) and
        every rank below it is "no call". In P/R/F a no-call row at the scored
        rank is a missed taxon (FN), never an FP. It is excluded from the
        taxon-set metrics (detection, sensitivity, specificity, Jaccard) and
        from L1/L2, which are computed over named calls only; the share of rows
        with no call is reported as ``unassigned_fraction``.
        Ranks above the deepest named rank that are blank (SPINGO writes
        ``d__; p__; c__; o__; f__; g__X; s__X y``) are not abstentions: they are
        not compared, as in v1.
    (b) Group labels. A label ``<rank>__[A|B|...]`` is correct at that rank if
        the truth taxon is one of the members and wrong otherwise. A row is a
        match at the scored rank only if this holds at every rank where both
        truth and call are named (as in v1). Groups from any classifier get the
        same rule.
    Taxon identity for the taxon-set metrics: a row that matches at the scored
    rank contributes the truth taxon; any other row named at that rank
    contributes rank prefix + name (a group not containing the truth is one
    opaque taxon). In v1 the whole label string was the identity, so a correct
    call whose upper ranks are blank (SPINGO writes ``d__; p__; ...``) never
    equalled the truth taxon and counted as a false-positive taxon plus a miss.

``v2-strict``
    As v2, except that a multi-member group is not credited: the group rank and
    every rank below it are treated as no call (the classifier did not name a
    single taxon there). This is the "strict single-species" view.

Group statistics (per classifier, mock, region; at species rank): number of
named species calls, the number and share that are groups, mean and maximum
group size, groups containing / not containing the truth species.
"""

import re
from collections import Counter

RANK_PREFIXES = ('d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__')
NOCALL_TOKENS = {'', 'AMBIGUOUS', 'Unclassified', 'Unassigned', 'No blast hit',
                 'unassigned', 'unclassified', 'other_nocall'}
RULES = ('v1', 'v2', 'v2-strict')
_GROUP = re.compile(r"^\[(.*)\]$")


def _strip_prefix(part):
    p = part.strip()
    if p[:3] in RANK_PREFIXES:
        return p[:3], p[3:].strip()
    return '', p


def is_nocall_content(content):
    return (content in NOCALL_TOKENS) or content.startswith('unclassified_')


def members(content):
    """Member set of one rank's content. ``[A|B]`` -> {A, B}; ``[A]`` -> {A};
    plain ``A`` -> {A}; no call -> empty set."""
    if content is None:
        return frozenset()
    s = content.strip()
    if is_nocall_content(s):
        return frozenset()
    m = _GROUP.match(s)
    if m:
        inner = [x.strip() for x in m.group(1).split('|') if x.strip()]
        # A bracket that is part of a real name, e.g. "[Ruminococcus] gnavus",
        # does not match ^\[...\]$ and is kept as a plain name.
        return frozenset(inner) if inner else frozenset()
    return frozenset([s])


def parse_lineage(tax_str):
    """Return a list of 7 member sets (domain..species). Strings without rank
    prefixes ('Unassigned', 'other') are handled: 'other' (truth for ASVs that
    match no expected sequence within 3 mismatches) is kept as a domain-rank
    name so it can never be matched, as in v1."""
    out = [frozenset()] * 7
    if tax_str is None:
        return out
    t = tax_str.strip()
    if t in NOCALL_TOKENS:
        return out
    parts = [x.strip() for x in t.split(';')]
    for i, part in enumerate(parts[:7]):
        _, content = _strip_prefix(part)
        out = out[:i] + [members(content)] + out[i + 1:]
    return out


def deepest_named(sets):
    for i in range(6, -1, -1):
        if sets[i]:
            return i
    return -1


def apply_strict(sets):
    """v2-strict: the first multi-member rank and every rank below it -> no call."""
    out = list(sets)
    for i, s in enumerate(out):
        if len(s) > 1:
            for j in range(i, 7):
                out[j] = frozenset()
            break
    return out


def classify(obs_sets, exp_sets, level):
    """Per-row outcome at ``level`` (0..6): 'match', 'underclassification',
    'overclassification' or 'misclassification'. Truth is single-valued (LCA)."""
    o = list(obs_sets[:level + 1])
    e = list(exp_sets[:level + 1])
    de = deepest_named(e + [frozenset()] * (6 - level))
    if de < 0:
        return 'match'
    do = deepest_named(o + [frozenset()] * (6 - level))
    if do < 0:
        return 'underclassification'
    for i in range(min(de, do) + 1):
        if e[i] and o[i] and not (e[i] & o[i]):
            return 'misclassification'
    if do < de:
        return 'underclassification'
    if do > de:
        return 'overclassification'
    return 'match'


def prf_from_outcomes(outcomes, rounding=3):
    """v1 counting: match -> TP; under -> FN; over/mis -> FP + FN.
    Returns (tp, fp, fn, p, r, f) with v1 rounding (3 decimals)."""
    tp = fp = fn = 0
    for c in outcomes:
        if c == 'match':
            tp += 1
        elif c == 'underclassification':
            fn += 1
        else:
            fp += 1
            fn += 1
    if tp == 0:
        return tp, fp, fn, 0, 0, 0
    p = round(tp / (tp + fp), rounding)
    r = round(tp / (tp + fn), rounding)
    f = round(2 * p * r / (p + r), rounding)
    return tp, fp, fn, p, r, f


def taxon_key(sets, level, raw):
    """Identity of a taxon at ``level`` for the taxon-set metrics: the rank
    prefix plus the name at that rank (names are unique within a rank in the
    GTDB-style Greengenes2 and GROND taxonomies), so a correct name is the same
    taxon whatever the classifier writes at the ranks above it. A group label
    is one opaque taxon. A truth label not named at ``level`` keeps its
    truncated string, as in v1."""
    s = sets[level]
    if len(s) == 1:
        return RANK_PREFIXES[level] + next(iter(s))
    if len(s) > 1:
        return 'GROUP|' + RANK_PREFIXES[level] + '[' + '|'.join(sorted(s)) + ']'
    return 'TRUNC|' + truncate(raw, level)


def truncate(tax_str, level):
    parts = [x.strip() for x in tax_str.split(';')]
    return ';'.join(parts[:level + 1])


def score_rows(pairs, level, rule, universe_size=None):
    """Score one (mock, region, method, param) at ``level`` under ``rule``
    ('v2' or 'v2-strict'). ``pairs``: dict id -> (truth_str, call_str).
    Returns a dict of metrics plus counts."""
    assert rule in ('v2', 'v2-strict')
    outcomes = []
    exp_keys, obs_keys = [], []
    n = 0
    n_nocall = 0
    for _id, (te, to) in pairs.items():
        n += 1
        es = parse_lineage(te)
        os_ = parse_lineage(to)
        if rule == 'v2-strict':
            os_ = apply_strict(os_)
        c = classify(os_, es, level)
        outcomes.append(c)
        ekey = taxon_key(es, level, te)
        exp_keys.append(ekey)
        if c == 'match':
            # credited (includes the case where truth is not named at this rank
            # and the call stops at the same depth)
            obs_keys.append(ekey)
        elif os_[level]:
            # a named (possibly wrong, possibly group) call at this rank
            obs_keys.append(taxon_key(os_, level, to))
        else:
            # no call at this rank: not a taxon (it may still be an error at a
            # higher rank, which P/R/F count and the higher rank's taxon set sees)
            n_nocall += 1
    tp, fp, fn, p, r, f = prf_from_outcomes(outcomes)
    exp_set, obs_set = set(exp_keys), set(obs_keys)
    tp_t = len(exp_set & obs_set)
    fp_t = len(obs_set - exp_set)
    fn_t = len(exp_set - obs_set)
    det = round(tp_t / max(1, len(exp_set)), 3)
    sens = round(tp_t / (tp_t + fn_t), 4) if (tp_t + fn_t) else 0.0
    union = exp_set | obs_set
    jac = round(tp_t / len(union), 4) if union else float('nan')
    if universe_size is None or universe_size <= 0:
        spec = float('nan')
    else:
        tn_t = max(0, universe_size - len(union))
        spec = round(tn_t / (tn_t + fp_t), 6) if (tn_t + fp_t) else 0.0
    l1, l2 = abundance_norms(exp_keys, obs_keys)
    return {
        'precision': p, 'recall': r, 'f-score': f,
        'detection_rate': det, 'sensitivity': sens, 'specificity': spec,
        'jaccard': jac, 'l1_error': l1, 'l2_error': l2,
        'unassigned_fraction': round(n_nocall / n, 6) if n else float('nan'),
        'tp': tp, 'fp': fp, 'fn': fn, 'n': n,
        'fp_taxa': fp_t, 'tp_taxa': tp_t, 'fn_taxa': fn_t,
    }


def abundance_norms(exp_taxa, obs_taxa):
    """v1 definition (L1 = sum |p_exp - p_obs|, L2 = sum (p_exp - p_obs)^2), with
    p_obs over the named calls only (the caller has dropped no-call rows)."""
    n_exp = max(1, len(exp_taxa))
    n_obs = max(1, len(obs_taxa))
    cexp, cobs = Counter(exp_taxa), Counter(obs_taxa)
    l1 = l2 = 0.0
    for t in set(cexp) | set(cobs):
        d = cexp.get(t, 0) / n_exp - cobs.get(t, 0) / n_obs
        l1 += abs(d)
        l2 += d * d
    return round(l1, 6), round(l2, 6)


def group_stats(pairs, level=6):
    """Group-label statistics at ``level`` (default species)."""
    n_named = n_group = n_in = n_out = 0
    sizes = []
    for _id, (te, to) in pairs.items():
        os_ = parse_lineage(to)
        es = parse_lineage(te)
        s = os_[level]
        if not s:
            continue
        n_named += 1
        if len(s) > 1:
            n_group += 1
            sizes.append(len(s))
            if es[level] and (es[level] & s):
                n_in += 1
            else:
                n_out += 1
    return {
        'n_named_calls': n_named,
        'n_group_calls': n_group,
        'group_share': round(n_group / n_named, 4) if n_named else 0.0,
        'group_size_mean': round(sum(sizes) / len(sizes), 3) if sizes else 0.0,
        'group_size_max': max(sizes) if sizes else 0,
        'group_truth_in': n_in,
        'group_truth_out': n_out,
    }


# ---------------------------------------------------------------------------
# Tie reconstruction for search-based classifiers (sensitivity option)
# ---------------------------------------------------------------------------

def tie_group_label(call_str, hit_taxa):
    """BLASTn / VSEARCH (QIIME 2 consensus) resolve a tie between top hits by
    truncating the call. Given the taxonomies of the hits tied at the best
    score, rebuild a group label in Taxonova's format: ranks where all tied hits
    agree keep the name, deeper ranks become ``<rank>__[A|B|...]``. Applied only
    when the emitted call has no species and the tied hits name at least two
    species; otherwise the call is returned unchanged."""
    if not hit_taxa:
        return call_str
    call = parse_lineage(call_str)
    if call[6]:
        return call_str
    lin = [parse_lineage(t) for t in hit_taxa]
    sp = set()
    for l in lin:
        sp |= set(l[6])
    if len(sp) < 2:
        return call_str
    parts = []
    for i, pref in enumerate(RANK_PREFIXES):
        names = sorted(set().union(*[l[i] for l in lin]))
        if not names:
            parts.append(pref)
        elif len(names) == 1:
            parts.append(pref + names[0])
        else:
            parts.append(pref + '[' + '|'.join(names) + ']')
    return '; '.join(parts)
