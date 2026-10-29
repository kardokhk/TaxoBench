
from pathlib import Path
from qiime2 import Artifact
import qiime2.plugins.demux.actions as demux_actions
import qiime2.plugins.dada2.actions as dada2_actions
import qiime2.plugins.taxa.actions as taxa_actions
from q2_types.feature_data import DNAIterator

import numpy
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from copy import copy, deepcopy
import pandas as pd
import os
import gzip
import biom
import skbio
from skbio import DNA
from skbio.io import write
from os.path import join, basename, dirname, exists, splitext
import random
import subprocess
import gzip
import shutil
import tempfile
import glob
from joblib import Parallel, delayed

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=DeprecationWarning)


def art_simulation(input_fasta, output_prefix, dir_out, seq_method = "MSv1", length = 200, count = 1,
                   n_samples = 1, mlen = 400, paired = False, seed=123, minQ=25, parallel_jobs=1,
                   force=False):
    """
    Simulate Miseq Illumina Reads with ART simulation tools [1]_.
    The function of amplicon simulation is used. More information `here <https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm>`_ .


    Arguments
    ---------
        input_fasta: str
            Reference fasta file
        output_prefix: str
            Output file prefix
        dir_out: str
            Output directory
        seq_method: str, default = 'MSv1'
            Sequencing method to simulate
        length: int, default = 200
            Length of the reads to simulate
        count: int, default = 1
            Number of reads per sequence in the fasta file.
        n_samples: int, default = 1
            Number of samples to generate
        mlen: int, default = 400
            Minimum total read length, only applicable if paired = True
        paired: bool, default = False
            If False, simulate single end reads
            If True, simulate paired end reads
        minQ: int, default = 25
            Minimal read quality to simulate


    Returns
    -------
    dict
        Dictionary of fasta paths for each sample {sample:fasta_path}

    References
    ----------
    .. [1] Weichun Huang, Leping Li, Jason R Myers, and Gabor T Marth. ART: a next-generation sequencing read simulator, Bioinformatics (2012) 28 (4): 593-594

    See Also
    --------
    eval_ref()

    """

    def _with_logging(cmd, idx, paired, log_dir=dir_out):
        log_dir = Path(log_dir)
        if paired:
            subdir = log_dir / "paired"
        else:
            subdir = log_dir / "single"
        subdir.mkdir(exist_ok=True, parents=True)
        log_file = subdir / f"art_sim_{idx}.log"
        with open(log_file, "a") as f:
            return subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)

    def _compress_and_rename(src, dst):
        with open(src, 'rb') as f_in, gzip.open(dst, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(src)

    # Create out_path
    if paired:
        dir_out_prefix = join(dir_out, "paired")
    else:
        dir_out_prefix= join(dir_out, "single")
    if not os.path.exists(dir_out_prefix):
        os.makedirs(dir_out_prefix)

    # Simulation
    commands = []
    outfiles_art = []

    for rep in range(1, n_samples + 1):
        sample_id = f"{output_prefix}-{rep}_S001"
        base_outname = join(dir_out_prefix, f"{sample_id}_L001")

        if paired:
            raw_out1 = base_outname + "1.fq"
            raw_out2 = base_outname + "2.fq"
            final_out1 = base_outname + "_R1_001.fq"
            final_out2 = base_outname + "_R2_001.fq"
            outname = base_outname  # passed to ART
            art_outname = [final_out1, final_out2]
        else:
            raw_out1 = base_outname + ".fq"
            final_out1 = base_outname + "_R1_001.fq"
            outname = base_outname
            art_outname = [final_out1]

        if not all(exists(f) for f in art_outname) or force:
            cmd = f"art_illumina -ss {seq_method} -amp -na -nf 0 -i {input_fasta} -l {length} -c {count} -rs {seed} -o {outname} --minQ {minQ}"
            if paired:
                cmd += f" -p -m {mlen}"
            commands.append((cmd.split(), rep))

        seed += 1
        outfiles_art.extend(art_outname)

    # Run ART simulations
    if commands:
        Parallel(n_jobs=parallel_jobs)(
            delayed(_with_logging)(cmd, idx, paired) for cmd, idx in commands
        )

    # Rename and compress output files to Illumina-style with .fq.gz
    for rep in range(1, n_samples + 1):
        sample_id = f"{output_prefix}-{rep}_S001"
        base_outname = join(dir_out_prefix, f"{sample_id}_L001")

        if paired:
            raw1 = base_outname + "1.fq"
            raw2 = base_outname + "2.fq"
            new1 = base_outname + "_R1_001.fastq.gz"
            new2 = base_outname + "_R2_001.fastq.gz"
            if exists(raw1):
                _compress_and_rename(raw1, new1)
            if exists(raw2):
                _compress_and_rename(raw2, new2)
        else:
            raw = base_outname + ".fq"
            new = base_outname + "_R1_001.fastq.gz"
            if exists(raw):
                _compress_and_rename(raw, new)

    return outfiles_art

def qiime_demultiplex(out_dir, fq_out_dir, paired=False, force=False):
    """
    Merges samples generated with art simulation tools using qiime2.
    More information in `here <https://docs.qiime2.org/2022.8/tutorials/moving-pictures-usage/>`_.


    Arguments
    ---------
        out_dir: str
            Output directory
        fq_out_dir : str
            Input directory, where the fastq files are located
        paired: bool, default = False
            If reads are single end, set to False. If reads are paired end, set to True
        force: bool, default = False
            If True, overwrites already existing results.

    Returns
    -------
        demux: qiime2 artifact
            Qiime2 artifact with the demultiplexed sequences
        demux_viz: qiime visualization
            Qiime2 visualization of the demux results
        qiime_out: str
            Output directory where the resulting files are located

    """

    # Outputs and parameter
    if paired:
        qiime_out = join(out_dir, "paired/")
        demux_out = join(qiime_out, "demux-paired-end.qza")
        end_type = 'SampleData[PairedEndSequencesWithQuality]'
    else:
        qiime_out = join(out_dir, "single/")
        demux_out = join(qiime_out, "demux-single-end.qza")
        end_type = 'SampleData[SequencesWithQuality]'
    if not os.path.exists(qiime_out):
        os.makedirs(qiime_out)

    # remove log files
    for filename in os.listdir(fq_out_dir):
        if filename.endswith(".log"):
            os.remove(os.path.join(fq_out_dir, filename))
            
    if not exists(demux_out) or force:
        # Demultiplex
        Artifact.import_data(end_type, fq_out_dir, 'CasavaOneEightSingleLanePerSampleDirFmt').save(demux_out)
        # summarize
        demux = Artifact.load(demux_out)
        demux_viz, = demux_actions.summarize(data=demux)
        demux_viz.save(splitext(demux_out)[0])

    return demux, demux_viz, qiime_out

def denoise_dada2(qiime_out,mockname, demux, paired = False,
                  trim_left = 0, trunc_len=0, trunc_len_f = 0, trunc_len_r = 0,
                  trim_left_f = 0, trim_left_r = 0, n_reads_learn = 1000000,
                  threads = 1, force = False):


    """ Dada2 denoising of Illumina reads using qiime2 implementation.
        More information in `here <https://docs.qiime2.org/2022.8/tutorials/moving-pictures-usage/>`_.

    Arguments
    ---------
        qiime_out: str
            Output directory
        mockname: str
            Label for the output files
        demux: qiime2 artifact
            Qiime2 artifact with demultiplexed sequences
        paired: bool, default = False
            If reads are single end, set to False. If reads are paired end, set to True
        trim_left: int, default = 0
            Position at which sequences should be trimmed due to
            low quality. This trims the 5' end of the of the
            input sequences, which will be the bases that were
            sequenced in the first cycles.
            Only applies if paired = False
        trunc_len: int, default = 0
            Position at which sequences should be truncated due
            to decrease in quality. This truncates the 3' end of
            the of the input sequences, which will be the bases
            that were sequenced in the last cycles. Reads that
            are shorter than this value will be discarded. If 0
            is provided, no truncation or length filtering will
            be performed.
            Only applies if paired = False
        trunc_len_f: int, default = 0
            Like trunc_len for the forward reads, only applies if paired = True
        trunc_len_r: int, default = 0
            Like trunc_len for the reverse reads, only applies if paired = True
        trim_left_f: int, default = 0
            Like trim_left for the forward reads, only applies if paired = True
        trim_left_r: int, default = 0
            Like trim_left for the reverse reads, only applies if paired = True
        n_reads_learn: int, default = 1000000
            The number of reads to use when training the error
            model. Smaller numbers will result in a shorter run
            time but a less reliable error model.
        threads: int, default = 1
            Number of processes to use
        force: bool, default = False
            If True, overwrite existing results.

    Returns
    --------
        table: qiime2 artifact
            Feature table
        rep_seqs: qiime2 artifact
            Representative sequences
        stats: qiime2 artifact
            Stats of the denoising process


    See also
    --------
    qiime_demultiplex()

    """

    # Outputs
    table_out = join(qiime_out, f"table_{mockname}")
    rep_seqs_out = join(qiime_out, f"rep_seqs_{mockname}")
    stats_out =  join(qiime_out, f"stats_{mockname}")

    if not os.path.exists(f"{table_out}.qza") or not os.path.exists(f"{rep_seqs_out}.qza") or not os.path.exists(f"{stats_out}.qza") or force:
        if paired:
            table, rep_seqs, stats = dada2_actions.denoise_paired(
                demultiplexed_seqs=demux,
                trunc_len_f = trunc_len_f,
                trunc_len_r=trunc_len_r,
                trim_left_f = trim_left_f,
                trim_left_r = trim_left_r,
                n_threads=threads,
                n_reads_learn = n_reads_learn
            )
        else:
            table, rep_seqs, stats = dada2_actions.denoise_single(
                demultiplexed_seqs=demux,
                trim_left=trim_left,
                trunc_len=trunc_len,
                n_threads=threads
            )
        ### Save files
        ## FEATURE-TABLE
        table.save(table_out)
        # to biom table
        biom_table = table.view(biom.Table)
        with biom.util.biom_open(f"{table_out}.biom", "w") as f:
            biom_table.to_hdf5(f, "Generated by QIIME 2 Python API")
        # to tsv table
        with open(f"{table_out}.tsv", "w") as f:
            biom_table.to_tsv(direct_io=f)
        ## REP-SEQS
        rep_seqs.save(rep_seqs_out)
        # to fasta
        seqs = rep_seqs.view(DNAIterator)
        with open(f"{rep_seqs_out}.fasta", "w") as f:
            for seq in seqs:
                write(seq, format='fasta', into=f)
        ## STATS
        stats.save(stats_out)
        # to tsv
        stats.export_data(join(qiime_out, "tmp"))
        shutil.copy(os.path.join(qiime_out, "tmp", "stats.tsv"), f"{stats_out}.tsv")
        shutil.rmtree(join(qiime_out, "tmp"))
    else:
        table = Artifact.load(f"{table_out}.qza")
        rep_seqs = Artifact.load(f"{rep_seqs_out}.qza")
        stats = Artifact.load(f"{stats_out}.qza")

    return table, rep_seqs, stats





# N handling in the expected (source) sequence. v1 (published) treated N as a
# wildcard (N_WILDCARD = True). v2 default: N is a mismatch. The matching code
# lives in truth_match.py (no QIIME 2 dependency) and is re-exported here.
N_WILDCARD = False

def single_character_matches(char, query, subject, n_wildcard=None):
    import truth_match
    nw = N_WILDCARD if n_wildcard is None else n_wildcard
    return truth_match.single_character_matches(char, query, subject, nw)

def gapless_mismatches(query, subject, n_wildcard=None):
    import truth_match
    nw = N_WILDCARD if n_wildcard is None else n_wildcard
    return truth_match.gapless_mismatches(query, subject, nw)

def compute_lca(taxonomies, delimiter='; '):
    split_taxa = [t.split(delimiter) for t in taxonomies]
    lca = []
    for rank_group in zip(*split_taxa):
        if all(rank == rank_group[0] for rank in rank_group):
            lca.append(rank_group[0])
        else:
            break
    return delimiter.join(lca) if lca else 'unassigned'

def compare_observed_and_expected(observed, expected, include_reverse_complement=True, return_strand=False,
                                  n_wildcard=None):
    ids = [obs.metadata['id'] for obs in observed]
    observed_arr = [numpy.array(list(str(obs))) for obs in observed]

    if include_reverse_complement:
        rev_observed_arr = [numpy.array(list(str(obs.reverse_complement()))) for obs in observed]

    keys = [' '.join([exp.metadata['id'], exp.metadata.get('description', '')]).strip()
            for exp in expected]
    expected_arr = [numpy.array(list(str(exp))) for exp in expected]

    comparison = {}
    for i, _id in enumerate(ids):
        obs = observed_arr[i]
        matches = [(gapless_mismatches(obs, exp, n_wildcard), key, 'forward')
                   for exp, key in zip(expected_arr, keys)]

        if include_reverse_complement:
            rev = rev_observed_arr[i]
            matches += [(gapless_mismatches(rev, exp, n_wildcard), key, 'reverse')
                        for exp, key in zip(expected_arr, keys)]

        matches.sort(key=lambda x: x[0])  # sort by mismatch count
        min_mismatches = matches[0][0]

        if return_strand:
            comparison[_id] = [(e, k, strand) for e, k, strand in matches if e == min_mismatches]
        else:
            comparison[_id] = [(e, k) for e, k, _ in matches if e == min_mismatches]

    return comparison

def generate_trueish_taxonomies_master(data_dir, feature_table,
                                expected_sequences, expected_taxonomy,
                                observed_sequences, trueish_taxonomy,
                                reverse_complement=True,
                                max_mismatches=3,
                                force=False,
                                n_wildcard=None):
    # run only if not exist
    if not os.path.exists(trueish_taxonomy) or force:
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
        mismatches = compare_observed_and_expected(observed, expected, 
                                                   include_reverse_complement=reverse_complement, 
                                                   return_strand=False,
                                                   n_wildcard=n_wildcard)

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


def generate_trueish_taxonomies(data_dir, feature_table,
                                expected_sequences, expected_taxonomy,
                                observed_sequences, trueish_taxonomy,
                                lca_mode=True,
                                reverse_complement=True,
                                max_mismatches=3,
                                force=False,
                                n_wildcard=None):
    # run only if not exist
    if not os.path.exists(trueish_taxonomy) or force:
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
        mismatches = compare_observed_and_expected(observed, expected, 
                                                   include_reverse_complement=reverse_complement, 
                                                   return_strand=False,
                                                   n_wildcard=n_wildcard)

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
                if lca_mode:
                    lca_taxonomy = compute_lca(matched_taxonomies)
                    writer.writerow([rep_id, lca_taxonomy])
                else:
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


def mutate_sequence(seq_dict, mutation_rate):
    """
    Mutate `mutation_rate * len(seq)` bases in the sequence stored in `seq_dict['seq']`.
    Returns a new dictionary with the same keys, but the sequence mutated.
    """
    seq_str = str(seq_dict['seq'])
    seq_len = len(seq_str)
    num_mutations = round(mutation_rate * seq_len)

    if num_mutations == 0:
        return seq_dict

    bases = ['A', 'T', 'C', 'G']
    seq_list = list(seq_str)

    # Randomly select positions to mutate
    positions = random.sample(range(seq_len), num_mutations)

    for pos in positions:
        original = seq_list[pos]
        if original not in bases:
            continue  # skip ambiguous or non-standard bases
        alternatives = [b for b in bases if b != original]
        seq_list[pos] = random.choice(alternatives)

    mutated_seq = Seq(''.join(seq_list))

    # Return a new dict with same taxa, mutated Seq
    return {
        'taxa': seq_dict['taxa'],
        'seq': mutated_seq
    }
