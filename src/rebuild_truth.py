"""Rebuild a trueish taxonomy with the taxobench v1 rule (gapless best offset,
both strands, all expected sequences at the minimum, LCA, 'other' above 3
mismatches), with N in the expected sequence as a wildcard (v1) or as a
mismatch (v2). Uses truth_match.best_hits (exact; no QIIME 2 needed).

Keys follow process_mock.compare_observed_and_expected: the expected key is
'<id> <description>'.strip(), looked up in the taxonomy map by that string.
"""
import csv
from collections import OrderedDict

import truth_match as T


def read_fasta(fp, keep_description=True):
    out = OrderedDict()
    key, buf = None, []
    with open(fp) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if key is not None:
                    out[key] = "".join(buf).upper()
                h = line[1:].strip()
                key = h if keep_description else h.split()[0]
                buf = []
            elif key is not None:
                buf.append(line.strip())
    if key is not None:
        out[key] = "".join(buf).upper()
    return out


def read_taxmap(fp):
    with open(fp) as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)
        return {row[0]: row[1] for row in reader if len(row) >= 2}


def read_truth(fp):
    out = OrderedDict()
    with open(fp) as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader)
        for row in reader:
            if len(row) >= 2:
                out[row[0]] = row[1]
    return out


def rebuild(observed_fp, expected_fp, taxonomy_fp, n_wildcard=False, max_mismatches=3,
            details=None):
    """Return OrderedDict asv_id -> label. If ``details`` is a list, append
    (asv, min_mm, tied_keys) per ASV."""
    obs = read_fasta(observed_fp, keep_description=False)
    exp = read_fasta(expected_fp, keep_description=True)
    tax = read_taxmap(taxonomy_fp)
    runs = {k: T._n_runs(s) for k, s in exp.items()} if n_wildcard else {}
    out = OrderedDict()
    for asv, q in obs.items():
        m, tied = T.best_hits(q, exp, max_mismatches=max_mismatches,
                              n_wildcard=n_wildcard, n_runs=runs)
        out[asv] = T.trueish_label(m, tied, tax, max_mismatches)
        if details is not None:
            details.append((asv, m, tied))
    return out


def write_truth(labels, fp):
    with open(fp, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["id", "taxonomy"])
        for k, v in labels.items():
            w.writerow([k, v])


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description=(
        "Rebuild one trueish taxonomy from the observed ASVs and the source "
        "sequences and taxonomy of a mock (no QIIME 2 needed). Default: N in a "
        "source sequence is a mismatch; --n-wildcard restores the published v1 rule."))
    ap.add_argument("--observed", required=True, help="rep_seqs_<mock>.fasta")
    ap.add_argument("--expected", required=True, help="input_data/<mock>_seqs.fasta")
    ap.add_argument("--taxonomy", required=True, help="input_data/<mock>_taxa.txt")
    ap.add_argument("--out", required=True, help="trueish taxonomy TSV to write")
    ap.add_argument("--n-wildcard", action="store_true")
    ap.add_argument("--max-mismatches", type=int, default=3)
    a = ap.parse_args(argv)
    write_truth(rebuild(a.observed, a.expected, a.taxonomy, n_wildcard=a.n_wildcard,
                        max_mismatches=a.max_mismatches), a.out)


if __name__ == "__main__":
    main()
