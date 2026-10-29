"""Gapless best-offset matching of ASVs to expected (source) sequences.

Factored out of ``process_mock.py`` (taxobench v1: ``single_character_matches``,
``gapless_mismatches``, ``compare_observed_and_expected``, ``compute_lca``) so the
truth builder can be run and tested without QIIME 2. The v1 functions in
``process_mock.py`` are kept and delegate here.

Semantics are those of v1, including two details that matter for exact
reproduction of the published truth files:

* ``numpy.convolve`` in full mode, so partial overlaps (query hanging off either
  end of the subject) are scored, with the overhang counted as mismatches;
* an ``N`` in the *query* scores one match against any of A/C/G/T.

The single change is the treatment of ``N`` in the *subject* (expected/source)
sequence. v1 (``n_wildcard=True``) scores it as a match for every query base,
so a source sequence with an N run at least as long as the ASV ties with every
ASV at 0 mismatches. v2 (``n_wildcard=False``, the default) scores it as a
mismatch.

``best_hits`` is an exact, faster equivalent of the v1 all-pairs search for the
question the truth builder asks (which expected sequences reach the minimum
mismatch count, when that minimum is <= max_mismatches). It is exact by the
pigeonhole argument documented in the function and is checked against the
brute-force v1 search in ``tests/test_truth_match.py``.
"""

import numpy

BASES = "ACGT"
_RC = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb")


def revcomp(seq):
    return seq.translate(_RC)[::-1]


# ---------------------------------------------------------------------------
# v1-compatible brute force (numpy.convolve over the full overlap range)
# ---------------------------------------------------------------------------

def single_character_matches(char, query, subject, n_wildcard=False):
    query = ((query == char) + (query == 'N')).astype(int)
    if n_wildcard:
        subject = ((subject == char) + (subject == 'N')).astype(int)
    else:
        subject = (subject == char).astype(int)
    return numpy.convolve(query[::-1], subject)


def gapless_mismatches(query, subject, n_wildcard=False):
    """``query`` and ``subject`` are numpy arrays of single characters."""
    matches = single_character_matches('A', query, subject, n_wildcard)
    for c in "CGT":
        matches += single_character_matches(c, query, subject, n_wildcard)
    return len(query) - matches.max()


def compute_lca(taxonomies, delimiter='; '):
    split_taxa = [t.split(delimiter) for t in taxonomies]
    lca = []
    for rank_group in zip(*split_taxa):
        if all(rank == rank_group[0] for rank in rank_group):
            lca.append(rank_group[0])
        else:
            break
    return delimiter.join(lca) if lca else 'unassigned'


def brute_force_min(query, subject, n_wildcard=False):
    """v1 gapless_mismatches on plain strings."""
    q = numpy.array(list(query))
    s = numpy.array(list(subject))
    return int(gapless_mismatches(q, s, n_wildcard))


# ---------------------------------------------------------------------------
# Exact fast search
# ---------------------------------------------------------------------------

def _score_at(query, subject, offset, n_wildcard):
    """Number of mismatches (v1 scoring) with query[0] aligned to subject[offset].
    Overhang positions count as mismatches. Returns len(query) - matches."""
    L = len(query)
    lo = max(0, -offset)
    hi = min(L, len(subject) - offset)
    if hi <= lo:
        return L
    q = numpy.frombuffer(query[lo:hi].encode(), numpy.uint8)
    s = numpy.frombuffer(subject[offset + lo:offset + hi].encode(), numpy.uint8)
    N = ord('N')
    qb = q != N
    sb = s != N
    # per-position match count, mirroring the sum over A/C/G/T convolutions
    acgt = numpy.frombuffer(BASES.encode(), numpy.uint8)
    q_is_base = numpy.isin(q, acgt)
    s_is_base = numpy.isin(s, acgt)
    m = ((q == s) & q_is_base).astype(numpy.int64)            # base == base
    m += ((~qb) & s_is_base)                                   # query N vs base
    if n_wildcard:
        m += (q_is_base & (~sb))                               # base vs subject N
        m += 4 * ((~qb) & (~sb))                               # N vs N (v1 quirk)
    return L - int(m.sum())


def _n_runs(seq):
    """List of (start, end) half-open N runs."""
    runs, i, n = [], 0, len(seq)
    while i < n:
        if seq[i] == 'N':
            j = i
            while j < n and seq[j] == 'N':
                j += 1
            runs.append((i, j))
            i = j
        else:
            i += 1
    return runs


def _candidate_offsets(query, subject, max_mismatches, n_wildcard, subj_n_runs):
    """Offsets (query[0] -> subject[offset]) that can reach <= max_mismatches.

    Pigeonhole: split the query into max_mismatches + 1 disjoint segments. Any
    alignment with <= max_mismatches mismatches leaves at least one segment
    with zero mismatches and no overhang. A zero-mismatch segment is an exact
    substring occurrence, except where it touches an N (query N scores as a
    match to any base; with n_wildcard a subject N does too). So the candidates
    are (i) offsets implied by exact occurrences of any segment, plus (ii) when
    the query holds an N, every offset (fallback), plus (iii) with n_wildcard,
    every offset whose window overlaps a subject N run."""
    L, S = len(query), len(subject)
    lo_off, hi_off = -max_mismatches, S - L + max_mismatches
    if 'N' in query:
        return range(lo_off, hi_off + 1)
    k = max_mismatches + 1
    bounds = [round(i * L / k) for i in range(k + 1)]
    cands = set()
    for a, b in zip(bounds[:-1], bounds[1:]):
        seg = query[a:b]
        if not seg:
            continue
        start = subject.find(seg)
        while start != -1:
            off = start - a
            if lo_off <= off <= hi_off:
                cands.add(off)
            start = subject.find(seg, start + 1)
    if n_wildcard:
        for (rs, re_) in subj_n_runs:
            # window [off, off+L) overlaps [rs, re_)
            for off in range(max(lo_off, rs - L + 1), min(hi_off, re_ - 1) + 1):
                cands.add(off)
    return cands


def _min_fft(query, subject, n_wildcard):
    """Minimum v1 mismatch count over every offset (full overlap range), via
    FFT cross-correlation of the per-base indicator vectors. Exact after
    rounding (counts are small integers)."""
    from scipy.signal import fftconvolve
    L = len(query)
    q = numpy.frombuffer(query.encode(), numpy.uint8)
    s = numpy.frombuffer(subject.encode(), numpy.uint8)
    qN = (q == ord('N')).astype(float)
    sN = (s == ord('N')).astype(float)
    tot = None
    for c in BASES:
        qc = (q == ord(c)).astype(float) + qN
        sc = (s == ord(c)).astype(float)
        if n_wildcard:
            sc = sc + sN
        m = fftconvolve(qc[::-1], sc, mode='full')
        tot = m if tot is None else tot + m
    return L - int(numpy.rint(tot.max()))


FFT_THRESHOLD = 300   # candidate offsets above which the FFT path is used


def best_hits(query, subjects, max_mismatches=3, n_wildcard=False,
              include_reverse_complement=True, n_runs=None):
    """Return (min_mismatches, [keys at the minimum]) over ``subjects``
    (dict key -> sequence string, upper case), both strands.

    If no alignment reaches <= max_mismatches, returns (None, []): v1 then
    writes 'other' and does not use the tie set. Exact with respect to the v1
    brute force whenever the minimum is <= max_mismatches."""
    if n_runs is None:
        n_runs = {k: _n_runs(s) for k, s in subjects.items()} if n_wildcard else {}
    queries = [query]
    if include_reverse_complement:
        queries.append(revcomp(query))
    best = {}
    for key, subj in subjects.items():
        runs = n_runs.get(key, [])
        m_best = None
        for q in queries:
            cands = _candidate_offsets(q, subj, max_mismatches, n_wildcard, runs)
            if len(cands) > FFT_THRESHOLD:
                mm = _min_fft(q, subj, n_wildcard)
                if m_best is None or mm < m_best:
                    m_best = mm
                continue
            for off in cands:
                mm = _score_at(q, subj, off, n_wildcard)
                if m_best is None or mm < m_best:
                    m_best = mm
        if m_best is not None and m_best <= max_mismatches:
            best[key] = m_best
    if not best:
        return None, []
    m = min(best.values())
    return m, sorted(k for k, v in best.items() if v == m)


def trueish_label(min_mm, tied_keys, tax_map, max_mismatches=3):
    """v1 labelling rule (lca_mode=True): 'other' above max_mismatches,
    otherwise the LCA of the set of tied expected taxonomies."""
    if min_mm is None or min_mm > max_mismatches:
        return 'other'
    taxa = set()
    for k in tied_keys:
        taxa.add(tax_map.get(k, 'NO-EXPECTED-TAXONOMY'))
    return compute_lca(list(taxa))
