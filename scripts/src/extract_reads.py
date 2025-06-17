# ----------------------------------------------------------------------------
# Copyright (c) 2016-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
#
# Modified script from Qiime2 plugin q2-feature-classifier
# Original script available at https://github.com/qiime2/q2-feature-classifier/blob/dev/q2_feature_classifier/_cutter.py
# ----------------------------------------------------------------------------

import os
import argparse
import subprocess
import skbio
from typing import Union
from numpy import ceil
from itertools import islice
from joblib import Parallel, delayed, effective_n_jobs

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=DeprecationWarning)


def _chunks(reads, chunk_size):
    reads = iter(reads)
    while True:
        chunk = list(islice(reads, chunk_size))
        if len(chunk) == 0:
            break
        yield chunk


def _autotune_reads_per_batch(reads, n_jobs):
    # detect effective jobs. Will raise error if n_jobs == 0
    if n_jobs == 0:
        raise ValueError("Value other than zero must be specified as number "
                         "of jobs to run.")
    else:
        n_jobs = effective_n_jobs(n_jobs)

    # we really only want to calculate this if running in parallel
    if n_jobs != 1:
        seq_count = subprocess.run(
            ['grep', '-c', '^>', str(reads)], check=True,
            stdout=subprocess.PIPE)
        # set a max value to avoid blowing up memory
        return min(int(ceil(int(seq_count.stdout.decode('utf-8')) / n_jobs)),
                   20000)
    # otherwise reads_per_batch = 20000, which has a modest memory overhead
    else:
        return 20000


def _seq_to_regex(seq):
    """Build a regex out of a IUPAC consensus sequence"""
    result = []
    for base in str(seq):
        if base in skbio.DNA.degenerate_chars:
            result.append('[{0}]'.format(
                ''.join(sorted(skbio.DNA.degenerate_map[base]))))
        else:
            result.append(base)

    return ''.join(result)


def _primers_to_regex(f_primer, r_primer):
    return '({0}.*{1})'.format(_seq_to_regex(f_primer),
                               _seq_to_regex(r_primer.reverse_complement()))


def _local_aln(primer, sequence):
    best_score = None
    for one_primer in sorted([str(s) for s in primer.expand_degenerates()]):
        # `sequence` may contain degenerates. These will usually be N
        # characters, which SSW will score as zero. Although undocumented, SSW
        # will treat other degenerate characters as a mismatch. We acknowledge
        # that this approach is a heuristic to finding an optimal alignment and
        # may be revisited in the future if there's an aligner that explicitly
        # handles degenerates.
        this_aln = \
            skbio.alignment.local_pairwise_align_ssw(skbio.DNA(one_primer),
                                                     sequence)
        score = this_aln[1]
        if best_score is None or score > best_score:
            best_score = score
            best_aln = this_aln
    return best_aln


def _semisemiglobal(primer, sequence, reverse=False):
    if reverse:
        primer = primer.reverse_complement()

    # locally align the primer
    (aln_prim, aln_seq), score, (prim_pos, seq_pos) = \
        _local_aln(primer, sequence)
    amplicon_pos = seq_pos[1]+len(primer)-prim_pos[1]

    # naively extend the alignment to be semi-global
    bits = [primer[:prim_pos[0]], aln_prim, primer[prim_pos[1]+1:]]
    aln_prim = ''.join(map(str, bits))
    bits = ['-'*(prim_pos[0]-seq_pos[0]),
            sequence[max(seq_pos[0]-prim_pos[0], 0):seq_pos[0]],
            aln_seq,
            sequence[seq_pos[1]+1:amplicon_pos],
            '-'*(amplicon_pos-len(sequence))]
    aln_seq = ''.join(map(str, bits))

    # count the matches
    matches = sum(s in skbio.DNA.degenerate_map.get(p, {p})
                  for p, s in zip(aln_prim, aln_seq))

    if reverse:
        amplicon_pos = max(seq_pos[0]-prim_pos[0], 0)

    return amplicon_pos, matches, len(aln_prim)


def _exact_match(seq, f_primer, r_primer):
    try:
        regex = _primers_to_regex(f_primer, r_primer)
        match = next(seq.find_with_regex(regex))
        beg, end = match.start + len(f_primer), match.stop - len(r_primer)
        return seq[beg:end]
    except StopIteration:
        return None


def _approx_match(seq, f_primer, r_primer, identity):
    beg, b_matches, b_length = _semisemiglobal(f_primer, seq)
    end, e_matches, e_length = _semisemiglobal(r_primer, seq, reverse=True)
    if (b_matches + e_matches) / (b_length + e_length) >= identity:
        return seq[beg:end]
    return None


def _gen_reads(sequence, f_primer, r_primer, trim_right, trunc_len, trim_left,
               identity, min_length, max_length, read_orientation):
    f_primer = skbio.DNA(f_primer)
    r_primer = skbio.DNA(r_primer)
    amp = None
    if read_orientation in ['forward', 'both']:
        amp = _exact_match(sequence, f_primer, r_primer)
    if not amp and read_orientation in ['reverse', 'both']:
        amp = _exact_match(sequence.reverse_complement(), f_primer, r_primer)
    if not amp and read_orientation in ['forward', 'both']:
        amp = _approx_match(sequence, f_primer, r_primer, identity)
    if not amp and read_orientation in ['reverse', 'both']:
        amp = _approx_match(
            sequence.reverse_complement(), f_primer, r_primer, identity)
    if not amp:
        return
    # we want to filter by max length before trimming
    if max_length > 0 and len(amp) > max_length:
        return
    if trim_right > 0:
        amp = amp[:-trim_right]
    if trunc_len > 0:
        amp = amp[:trunc_len]
    if trim_left > 0:
        amp = amp[trim_left:]
    if min_length > 0 and len(amp) < min_length:
        return
    if not amp:
        return
    return amp


def cutter(sequences: str, f_primer: str,
           r_primer: str, trim_right: int = 0,
           trunc_len: int = 0, trim_left: int = 0,
           min_length: int = 0, max_length: int = 0,
           identity: float = 0.8, n_jobs: int = 1,
           batch_size: Union[int, str] = 'auto', read_orientation: str = 'both',
           output_path: str = "extracted_reads.fasta") \
           -> str:
    """Extract the read selected by a primer or primer pair. Only sequences
    which match the primers at greater than the specified identity are returned

    Parameters
    ----------
    sequences : str
        Path to input FASTA or sequence file.
    f_primer : skbio.sequence.DNA
        Forward primer sequence
    r_primer : skbio.sequence.DNA
        Reverse primer sequence
    trim_right : int, optional
        `trim_right` nucleotides are removed from the 3' end if trim_right is
        positive. Applied before trunc_len.
    trunc_len : int, optional
        Read is cut to trunc_len if trunc_len is positive. Applied after
        trim_right.
    trim_left : int, optional
        `trim_left` nucleotides are removed from the 5' end if trim_left is
        positive. Applied after trim_right and trunc_len.
    min_length: int, optional
        Minimum amplicon length. Shorter amplicons are discarded. Default: 0
    max_length: int, optional
        Maximum amplicon length. Longer amplicons are discarded.
    identity : float, optional
        Minimum combined primer match identity threshold. Default: 0.8
    n_jobs: int, optional
        Number of seperate processes to break the task into.
    batch_size: int, optional
        Number of samples to be processed in one batch.
    read_orientation: str, optional
        'Orientation of primers relative to the sequences: "forward" searches '
        'for primer hits in the forward direction, "reverse" searches the '
        'reverse-complement, and "both" searches both directions.'
    Returns
    -------
    - str: Path to output FASTA file with extracted reads.
    """
    if min_length > trunc_len - (trim_left + trim_right) and trunc_len > 0:
        raise ValueError('The minimum length setting is greater than the '
                         'length of the truncated sequences. This will cause '
                         'all sequences to be removed from the dataset. To '
                         'proceed, set '
                         'min_length ≤ trunc_len - (trim_left  + '
                         'trim_right).')

    n_jobs = effective_n_jobs(n_jobs)
    if batch_size == 'auto':
        batch_size = _autotune_reads_per_batch(sequences, n_jobs)

    sequences = skbio.read(sequences, format='fasta', constructor=skbio.DNA)

    with open(output_path, 'w') as fh:
        with Parallel(n_jobs) as parallel:
            for chunk in _chunks(sequences, batch_size):
                amplicons = parallel(delayed(_gen_reads)(sequence, f_primer,
                                                         r_primer,
                                                         trim_right,
                                                         trunc_len,
                                                         trim_left,
                                                         identity,
                                                         min_length,
                                                         max_length,
                                                         read_orientation)
                                     for sequence in chunk)
                for amplicon in amplicons:
                    if amplicon is not None:
                        skbio.write(amplicon, format='fasta', into=fh)

    if os.stat(output_path).st_size == 0:
        raise RuntimeError("No matching amplicons were found.")

    return output_path




def main():
    parser = argparse.ArgumentParser(
        description='Extract amplicon reads from a reference database using in-silico PCR.'
    )

    # Input fasta sequence
    parser.add_argument('--sequences', required=True,
                        help='Path to the input reference sequences in fasta format.')

    # Output fasta sequence
    parser.add_argument('--output', required=True,
                        help='Path to write the extracted reads in fasta format.')

    # Primers
    parser.add_argument('--f_primer', required=True,
                        help='Forward primer sequence (5\' -> 3\').')
    parser.add_argument('--r_primer', required=True,
                        help='Reverse primer sequence (5\' -> 3\'). Do not use reverse-complemented primer.')

    # Trimming and truncation
    parser.add_argument('--trim_right', type=int, default=0,
                        help='Number of bases to trim from the 3\' end before truncation.')
    parser.add_argument('--trunc_len', type=int, default=0,
                        help='Truncate read to this length after trim_right.')
    parser.add_argument('--trim_left', type=int, default=0,
                        help='Number of bases to trim from the 5\' end after trunc_len.')

    # Length filtering
    parser.add_argument('--min_length', type=int, default=0,
                        help='Minimum read length to retain after all trimming steps.')
    parser.add_argument('--max_length', type=int, default=0,
                        help='Maximum read length allowed before trimming.')

    # Identity threshold
    parser.add_argument('--identity', type=float, default=0.8,
                        help='Minimum combined primer match identity threshold.')

    # Parallelization
    parser.add_argument('--n_jobs', type=int, default=1,
                        help='Number of parallel processes to run (>=1).')
    parser.add_argument('--batch_size', default='auto',
                        help='Batch size for processing. Can be an integer or "auto".')

    # Read orientation
    parser.add_argument('--read_orientation', choices=['both', 'forward', 'reverse'],
                        default='both',
                        help='Orientation of primers relative to the reference sequences.')

    args = parser.parse_args()

    # Example usage of parsed args
    print("Running amplicon extraction with parameters:")
    for arg, value in vars(args).items():
        print(f"  {arg}: {value}")

    # run
    cutter(
        sequences=args.sequences,
        f_primer=args.f_primer,
        r_primer=args.r_primer,
        trim_right=args.trim_right,
        trunc_len=args.trunc_len,
        trim_left=args.trim_left,
        min_length=args.min_length,
        max_length=args.max_length,
        identity=args.identity,
        n_jobs=args.n_jobs,
        batch_size=args.batch_size,
        read_orientation=args.read_orientation,
        output_path=args.output
    )

if __name__ == "__main__":
    main()
