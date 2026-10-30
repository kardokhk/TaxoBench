import os
import subprocess
from os.path import join, exists


def _query_artifact(rep_seq_base, classifier):
    """Path of the imported query artifact for one classifier.

    Derived from the (input, classifier) pair, not from the input alone. All
    three QIIME classifiers previously imported to the same "<base>.qza": run
    concurrently they truncated each other mid-write, which surfaced as
    "Archive does not contain a correctly formatted VERSION file" and read like
    a broken installation.
    """
    return "%s.%s.qza" % (rep_seq_base, classifier)


def _stale(artifact, rep_seq_base):
    """True when the artifact is missing or older than the FASTA it came from.

    The previous guard was `if not os.path.exists(...)`, which never noticed
    that the source had changed. data/real_world_input/rep_seqs_zuo.qza held
    1,769 sequences against a 3,218-sequence FASTA written five hours later, so
    three classifiers were scored on a different query set from every other
    program, and because feature ids are positional the evaluator paired counts
    with the wrong sequences' taxonomy rather than failing.
    """
    if not os.path.exists(artifact):
        return True
    src = rep_seq_base + ".fasta"
    return os.path.exists(src) and os.path.getmtime(src) > os.path.getmtime(artifact)

# QIIME 2 is imported only by the builders that need it, so that this module
# (and the dry-run and --help paths of the scripts) loads without QIIME 2.
def _qiime2():
    from qiime2 import Artifact
    from qiime2.plugins.feature_classifier.methods import makeblastdb
    from qiime2.plugins.rescript.pipelines import evaluate_fit_classifier
    return Artifact, makeblastdb, evaluate_fit_classifier


#: Ranks above species, in order, with their taxonova threshold flag names.
TAXONOVA_PARENT_RANKS = ("domain", "phylum", "class", "order", "family", "genus")

#: Defaults for the "taxonova" block of a classifier config. Used when a config
#: has no such block, which keeps older configs working unchanged.
TAXONOVA_DEFAULTS = {
    "cli": "classification",
    "models": {},
    "parent_gates": None,
}


def taxonova_settings(classifiers_cfg):
    """The "taxonova" block of a classifier config, with defaults filled in."""
    out = dict(TAXONOVA_DEFAULTS)
    out.update((classifiers_cfg or {}).get("taxonova", {}))
    return out


def resolve_taxonova_model(region, region_dir, settings, ref_root=None):
    """Model for one region: a model id (for example "v4"), a file path, or
    None when nothing is configured or found.

    settings["models"][region] may be
      - a model id such as "v4" or "16s", passed to taxonova as is (the
        package downloads and verifies it on first use);
      - a path to a checkpoint file, absolute or relative to the reference
        database root (for example "V4/taxonova/gg2_v4_mrgroups_model.pt");
      - a dict {"path": ..., "md5": ..., "url": ...}; only "path" is used here.
    Without an entry, the first *model.pt under <region_dir>/taxonova/ (sorted)
    is used, as in earlier versions of this benchmark.
    """
    import glob
    entry = (settings.get("models") or {}).get(region)
    spec = entry if isinstance(entry, dict) else {}
    if isinstance(entry, dict):
        entry = entry.get("path")
    if entry:
        looks_like_path = entry.endswith(".pt") or os.sep in entry
        if not looks_like_path:
            return entry
        if not os.path.isabs(entry) and ref_root:
            entry = join(ref_root, entry)
        _check_pinned_model(entry, spec)
        return entry
    models = sorted(glob.glob(join(region_dir, "taxonova", "*model.pt")))
    return models[0] if models else None


_VERIFIED_MODELS = {}


def _check_pinned_model(path, spec):
    """Verify a pinned checkpoint against its md5 (once per process).
    A missing file is reported with its download URL; a file with the wrong
    md5 raises ValueError, because it would not reproduce the pinned runs."""
    want = spec.get("md5")
    if not want:
        return
    if not os.path.exists(path):
        url = spec.get("url")
        print(f"[taxonova] pinned model not found: {path}"
              + (f" (download: {url})" if url else ""))
        return
    if path not in _VERIFIED_MODELS:
        import hashlib
        h = hashlib.md5()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
        _VERIFIED_MODELS[path] = h.hexdigest()
    if _VERIFIED_MODELS[path] != want:
        raise ValueError(f"{path}: md5 {_VERIFIED_MODELS[path]} does not match "
                         f"the pinned {want}")


def build_taxonova_command(
    rep_seq_base, outdir, param_dict, reference_dict,
    scripts_dir, conda_env, dry_run, settings=None
):
    """Command line for one taxonova run.

    settings is the "taxonova" block of the classifier config
    (see taxonova_settings):
      cli           "classify" (taxonova >= 26.10.1: model ids, package
                    default thresholds) or "classification" (the original
                    CLI, used by the published runs; works in every release)
      parent_gates  None to keep the package defaults for the ranks above
                    species, or a dict {rank: threshold} (for example 0.5 at
                    every parent rank, as in the published runs)
    The species threshold comes from param_dict["confidence_threshold"], so
    it names the results directory (confidence_threshold-<value>).
    """
    settings = taxonova_settings({"taxonova": settings or {}})
    cli = settings["cli"]
    if cli not in ("classify", "classification"):
        raise ValueError(f"taxonova cli must be 'classify' or 'classification', not {cli!r}")

    input_sequences = f"{rep_seq_base}.fasta"
    model = reference_dict.get('model')
    taxonomy_raw = join(outdir, 'taxonomy_raw.tsv')
    taxonomy_final = join(outdir, 'taxonomy_final.tsv')

    # No configured model and no checkpoint on disk: train one on the region's
    # reference (slow; a GPU is recommended).
    if model is None:
        reference_dir = join(os.path.dirname(reference_dict['seqs']), 'taxonova')
        model = join(reference_dir, 'taxonova_model.pt')
        if not dry_run:
            os.makedirs(reference_dir, exist_ok=True)
            reference_seqs = reference_dict.get('merged') or reference_dict['seqs']
            subprocess.run(["taxonova", "training", "--verbose", "--smart_full_data",
                            "--train_db", reference_seqs,
                            "--model_output", model], check=True)

    # Thresholds. The parent gates are stated explicitly when the config asks
    # for them, rather than left to the package default: the published runs
    # used a flat 0.5 at every parent rank, while taxonova >= 26.9 scales a
    # graded profile (0.5 domain to family, 0.7 genus, 0.8 species) with the
    # species threshold. Keeping the gates in the "taxonova" block rather than
    # in param_dict keeps them out of the param_id that names the results
    # directory.
    gates = {}
    for rank, value in (settings.get("parent_gates") or {}).items():
        if rank not in TAXONOVA_PARENT_RANKS:
            raise ValueError(f"unknown taxonova parent rank {rank!r}")
        gates[rank] = value
    for rank in TAXONOVA_PARENT_RANKS:          # a param_dict entry wins
        if f"{rank}_threshold" in param_dict:
            gates[rank] = param_dict[f"{rank}_threshold"]
    species = param_dict.get("confidence_threshold")

    if cli == "classify":
        threshold_args = f"--confidence {species} " if species is not None else ""
        threshold_args += "".join(f"--{r}_threshold {gates[r]} "
                                  for r in TAXONOVA_PARENT_RANKS if r in gates)
        run = (f"taxonova classify {input_sequences} "
               f"--model {model} "
               f"{threshold_args}"
               f"-o {taxonomy_raw}")
    else:
        threshold_args = "".join(f"--{r}_threshold {gates[r]} "
                                 for r in TAXONOVA_PARENT_RANKS if r in gates)
        if species is not None:
            threshold_args += f"--confidence_threshold {species} "
        run = (f"taxonova classification "
               f"{threshold_args}"
               f"--input {input_sequences} "
               f"--model_input {model} "
               f"--output {taxonomy_raw}")

    cmd = (
        f"conda run -p {conda_env} bash -c \""
        f"mkdir -p {outdir} && "
        f"{run} && "
        f"python {scripts_dir}/process_taxonova_output.py "
            f"{taxonomy_raw} {taxonomy_final} False \""
    )
    return cmd



def build_vsearch_command(
    rep_seq_base, outdir, param_dict,
    reference_dict, conda_env, n_threads, dry_run
):
    # reference database dir
    reference_dir = join(os.path.dirname(reference_dict['seqs']), 'vsearch')
    if not exists(reference_dir):
        if not dry_run:
            os.makedirs(reference_dir)

    # build reference/input/output paths
    input_sequences = _query_artifact(rep_seq_base, "vsearch")
    reference_seqs = join(reference_dir, "sequences.qza")
    reference_taxa = join(reference_dir, "taxonomy.qza")
    classification_out = join(outdir, 'taxonomy_classification_out.qza')
    search_out = join(outdir, 'taxonomy_search_out.qza')
    taxonomy_final = join(outdir, 'taxonomy_final.tsv')

    # Import reference sequences and taxonomy
    if not all(os.path.exists(f) for f in [reference_seqs, reference_taxa]):
        if not dry_run:
            # Import seqs
            Artifact, makeblastdb, evaluate_fit_classifier = _qiime2()
            ref_seqs = Artifact.import_data(
                'FeatureData[Sequence]',
                reference_dict['seqs']
            )
            ref_seqs.save(reference_seqs)

            # Import taxonomy
            ref_taxonomy = Artifact.import_data(
                'FeatureData[Taxonomy]',
                reference_dict['taxa'],
                'HeaderlessTSVTaxonomyFormat'
            )
            ref_taxonomy.save(reference_taxa)

    # Check if input QZA file exists, if not create it
    if _stale(input_sequences, rep_seq_base):
        if not dry_run:
            # Use helper script to import sequences
            import subprocess
            helper_script = join(os.path.dirname(__file__), 'qiime_helper.py')
            cmd_import = f'conda run -p {conda_env} python {helper_script} import_fasta_to_qza {rep_seq_base}.fasta {input_sequences}'
            subprocess.run(cmd_import, shell=True, check=True)

    # assemble command using helper script for conversion
    helper_script = join(os.path.dirname(__file__), 'qiime_helper.py')
    cmd = (
        f'conda run -p {conda_env} bash -c '
        f'"mkdir -p {outdir} && '
        f'qiime feature-classifier classify-consensus-vsearch '
        f'--i-query {input_sequences} '
        f'--i-reference-reads {reference_seqs} '
        f'--i-reference-taxonomy {reference_taxa} '
        f'--p-threads {n_threads} '
        f'--p-strand both '
        f'--p-query-cov 0.80 '
        f'--p-maxaccepts {param_dict["p-maxaccepts"]} '
        f'--p-perc-identity {param_dict["p-perc-identity"]} '
        f'--p-min-consensus {param_dict["p-min-consensus"]} '
        f'--o-classification {classification_out} '
        f'--o-search-results {search_out} && '
        f'python {helper_script} convert_qza_to_tsv {classification_out} {taxonomy_final}"'
    )
        
    return cmd



def build_blast_command(
    rep_seq_base, outdir, param_dict,
    reference_dict, conda_env, n_threads, dry_run
):
    # reference database dir
    reference_dir = join(os.path.dirname(reference_dict['seqs']), 'blast')
    if not exists(reference_dir):
        if not dry_run:
            os.makedirs(reference_dir)

    # build reference/input/output paths
    input_sequences = _query_artifact(rep_seq_base, "blast")
    reference_seqs = join(reference_dir, "sequences.qza")
    reference_taxa = join(reference_dir, "taxonomy.qza")
    reference_db = join(reference_dir, "blastdb.qza")
    classification_out = join(outdir, 'taxonomy_classification_out.qza')
    search_out = join(outdir, 'taxonomy_search_out.qza')
    taxonomy_final = join(outdir, 'taxonomy_final.tsv')

    # Import reference sequences and taxonomy
    if not all(os.path.exists(f) for f in [reference_seqs, reference_taxa, reference_db]):
        if not dry_run:
            # Import seqs
            Artifact, makeblastdb, evaluate_fit_classifier = _qiime2()
            ref_seqs = Artifact.import_data(
                'FeatureData[Sequence]',
                reference_dict['seqs']
            )
            ref_seqs.save(reference_seqs)

            # Import taxonomy
            ref_taxonomy = Artifact.import_data(
                'FeatureData[Taxonomy]',
                reference_dict['taxa'],
                'HeaderlessTSVTaxonomyFormat'
            )
            ref_taxonomy.save(reference_taxa)

            # Make blastdb
            ref_db = makeblastdb(ref_seqs)
            ref_db.database.save(reference_db)

    # Check if input QZA file exists, if not create it
    if _stale(input_sequences, rep_seq_base):
        if not dry_run:
            # Use helper script to import sequences
            import subprocess
            helper_script = join(os.path.dirname(__file__), 'qiime_helper.py')
            cmd_import = f'conda run -p {conda_env} python {helper_script} import_fasta_to_qza {rep_seq_base}.fasta {input_sequences}'
            subprocess.run(cmd_import, shell=True, check=True)

    # assemble command using helper script for conversion
    helper_script = join(os.path.dirname(__file__), 'qiime_helper.py')
    cmd = (
        f'conda run -p {conda_env} bash -c '
        f'"mkdir -p {outdir} && '
        f'qiime feature-classifier classify-consensus-blast '
        f'--i-query {input_sequences} '
        f'--i-blastdb {reference_db} '
        f'--i-reference-taxonomy {reference_taxa} '
        f'--p-num-threads {n_threads} '
        f'--p-strand both '
        f'--p-query-cov 0.80 '
        f'--p-evalue 0.001 '
        f'--p-maxaccepts {param_dict["p-maxaccepts"]} '
        f'--p-perc-identity {param_dict["p-perc-identity"]} '
        f'--p-min-consensus {param_dict["p-min-consensus"]} '
        f'--o-classification {classification_out} '
        f'--o-search-results {search_out} && '
        f'python {helper_script} convert_qza_to_tsv {classification_out} {taxonomy_final}"'
    )
        
    return cmd



def build_naive_bayes_command(
    rep_seq_base, outdir, param_dict,
    reference_dict, conda_env, n_threads, dry_run
):
    # reference database dir
    reference_dir = join(os.path.dirname(reference_dict['seqs']), 'naive-bayes')
    if not exists(reference_dir):
        if not dry_run:
            os.makedirs(reference_dir)

    # build reference/input/output paths
    input_sequences = _query_artifact(rep_seq_base, "naive_bayes")
    model_file = reference_dict['model']
    evaluation_file = join(reference_dir, "evaluation.qza")
    prediction_file = join(reference_dir, "predictions.qza")
    classification_out = join(outdir, 'taxonomy_classification_out.qza')
    taxonomy_final = join(outdir, 'taxonomy_final.tsv')

    # Import reference sequences and taxonomy and build classifier
    if model_file is None:
        if not dry_run:
            # Import seqs
            Artifact, makeblastdb, evaluate_fit_classifier = _qiime2()
            ref_seqs = Artifact.import_data(
                'FeatureData[Sequence]',
                reference_dict['seqs']
            )

            # Import taxonomy
            ref_taxonomy = Artifact.import_data(
                'FeatureData[Taxonomy]',
                reference_dict['taxa'],
                'HeaderlessTSVTaxonomyFormat'
            )

            # Train classifier
            classifier, evaluation, observed_taxonomy = (
            evaluate_fit_classifier(
                sequences = ref_seqs, 
                taxonomy = ref_taxonomy,
                reads_per_batch = 'auto',
                n_jobs = n_threads,
                confidence = 0.7
                )
            )

            # Save outputs
            classifier.save(model_file)
            evaluation.save(evaluation_file)
            observed_taxonomy.save(prediction_file)

    # Check if input QZA file exists, if not create it
    if _stale(input_sequences, rep_seq_base):
        if not dry_run:
            # Use helper script to import sequences
            import subprocess
            helper_script = join(os.path.dirname(__file__), 'qiime_helper.py')
            cmd_import = f'conda run -p {conda_env} python {helper_script} import_fasta_to_qza {rep_seq_base}.fasta {input_sequences}'
            subprocess.run(cmd_import, shell=True, check=True)

    # assemble command using helper script for conversion
    helper_script = join(os.path.dirname(__file__), 'qiime_helper.py')
    cmd = (
        f'conda run -p {conda_env} bash -c '
        f'"mkdir -p {outdir} && '
        f'qiime feature-classifier classify-sklearn '
        f'--i-reads {input_sequences} '
        f'--i-classifier {model_file} '
        f'--p-n-jobs {n_threads} '
        f'--p-pre-dispatch "2*n_jobs" '
        f'--p-confidence {param_dict["p-confidence"]} '
        f'--o-classification {classification_out} && '
        f'python {helper_script} convert_qza_to_tsv {classification_out} {taxonomy_final}"'
    )
        
    return cmd


def build_idtaxa_command(
    rep_seq_base, outdir, param_dict, reference_dict, 
    scripts_dir, conda_env, dry_run
):
    # reference database dir
    reference_dir = join(os.path.dirname(reference_dict['seqs']), 'idtaxa')
    if not exists(reference_dir):
        if not dry_run:
            os.makedirs(reference_dir)

    # build reference/input/output paths
    input_sequences = f"{rep_seq_base}.fasta"
    model_file = reference_dict['model']
    taxonomy_raw = join(outdir, 'taxonomy_raw.tsv')
    taxonomy_final = join(outdir, 'taxonomy_final.tsv')
    
    # Get confidence threshold from parameters
    confidence_threshold = param_dict.get('p-confidence', 50)

    # assemble command - first create raw output, then process to final
    cmd = (
        f"conda run -p {conda_env} bash -c \""
        f"mkdir -p {outdir} && "
        f"Rscript {scripts_dir}/run_idtaxa_classification.R "
            f"{input_sequences} "
            f"{model_file} "
            f"{taxonomy_raw} "
            f"{confidence_threshold} && "
        f"python {scripts_dir}/process_idtaxa_output.py "
            f"{taxonomy_raw} {taxonomy_final}\""
    )

    return cmd


def build_spingo_command(
    rep_seq_base, outdir, param_dict, reference_dict, 
    scripts_dir, conda_env, n_threads, dry_run
):
    """
    Build SPINGO classification command with dual environment support
    
    SPINGO parameters:
    --database (-d): Location of the fasta format database
    --input (-i): Location of the input file
    --kmersize (-k): Kmer size (default=8)
    --bootstrap (-b): Number of bootstrap samples (default=10)
    --subsample (-s): Proportion of kmers for bootstrap (default=kmersize)
    --processors (-p): Number of processor threads (default=1)
    
    Note: SPINGO requires a separate environment from the main processing environment
    """
    # reference database dir
    reference_dir = join(os.path.dirname(reference_dict['seqs']), 'spingo')
    if not exists(reference_dir):
        if not dry_run:
            os.makedirs(reference_dir)

    # build reference/input/output paths
    input_sequences = f"{rep_seq_base}.fasta"
    database_file = reference_dict['database']  # SPINGO database file
    taxonomy_raw = join(outdir, 'taxonomy_raw.tsv')
    taxonomy_final = join(outdir, 'taxonomy_final.tsv')
    
    # Get SPINGO parameters from param_dict
    kmersize = param_dict.get('kmersize', 8)
    bootstrap = param_dict.get('bootstrap', 10)
    subsample = param_dict.get('subsample', kmersize)

    # Handle environment configuration
    # Check if conda_env is a dictionary with separate environments for SPINGO
    if isinstance(conda_env, dict):
        spingo_env = conda_env['spingo_env']
        processing_env = conda_env['processing_env']
        
        # Build command with environment switching
        cmd = (
            f"bash -c \""
            f"mkdir -p {outdir} && "
            f"conda run -p {spingo_env} spingo "
                f"--database {database_file} "
                f"--input {input_sequences} "
                f"--processors {n_threads} "
                f"--kmersize {kmersize} "
                f"--bootstrap {bootstrap} "
                f"--subsample {subsample} "
                f"> {taxonomy_raw} && "
            f"conda run -p {processing_env} python {scripts_dir}/process_spingo_output.py "
                f"{taxonomy_raw} {taxonomy_final}\""
        )
    else:
        # Fallback to single environment (for backward compatibility)
        cmd = (
            f"conda run -p {conda_env} bash -c \""
            f"mkdir -p {outdir} && "
            f"spingo "
                f"--database {database_file} "
                f"--input {input_sequences} "
                f"--processors {n_threads} "
                f"--kmersize {kmersize} "
                f"--bootstrap {bootstrap} "
                f"--subsample {subsample} "
                f"> {taxonomy_raw} && "
            f"python {scripts_dir}/process_spingo_output.py "
                f"{taxonomy_raw} {taxonomy_final}\""
        )

    return cmd


def build_kraken2_command(
    rep_seq_base, outdir, param_dict,
    reference_dict, scripts_dir, conda_env, n_threads, dry_run
):
    """
    Build Kraken2 classification command.
    
    Args:
        rep_seq_base: Base path for representative sequences (without extension)
        outdir: Output directory for results
        param_dict: Dictionary of parameters (confidence, quick, etc.)
        reference_dict: Dictionary containing database path
        scripts_dir: Path to scripts directory
        conda_env: Conda environment path
        n_threads: Number of threads to use
        dry_run: Whether this is a dry run
        
    Returns:
        str: Command string to execute
    """
    # Build reference/input/output paths
    input_sequences = f"{rep_seq_base}.fasta"
    database_path = reference_dict['database']  # Path to kraken2 database directory
    
    # Output files
    classification_out = join(outdir, 'kraken2_classification.txt')
    report_out = join(outdir, 'kraken2_report.txt')
    unclassified_out = join(outdir, 'kraken2_unclassified.fasta')
    classified_out = join(outdir, 'kraken2_classified.fasta')
    taxonomy_final = join(outdir, 'taxonomy_final.tsv')
    
    # Extract parameters with defaults
    confidence = param_dict.get('confidence', 0.0)
    quick = param_dict.get('quick', False)
    use_names = param_dict.get('use_names', True)
    
    # Build Kraken2 command
    kraken2_cmd_parts = [
        "kraken2",
        f"--db {database_path}",
        f"--threads {n_threads}",
        f"--confidence {confidence}",
        f"--output {classification_out}",
        f"--report {report_out}",
        f"--unclassified-out {unclassified_out}",
        f"--classified-out {classified_out}"
    ]
    
    if quick:
        kraken2_cmd_parts.append("--quick")
    
    if use_names:
        kraken2_cmd_parts.append("--use-names")
    
    kraken2_cmd_parts.append(input_sequences)
    
    kraken2_cmd = " ".join(kraken2_cmd_parts)
    
    # Assemble full command
    cmd = (
        f"conda run -p {conda_env} bash -c \""
        f"mkdir -p {outdir} && "
        f"{kraken2_cmd} && "
        f"python {scripts_dir}/process_kraken2_output.py "
            f"{classification_out} {taxonomy_final} {confidence}\""
    )
    
    return cmd


def build_mapseq_command(
    rep_seq_base, outdir, param_dict, reference_dict, 
    scripts_dir, conda_env, n_threads, dry_run
):
    """
    Build MAPseq classification command
    
    MAPseq command format:
    mapseq [options] input.fa reference.fa reference.tax > output.mseq
    
    Parameters from param_dict:
    - nthreads: number of threads to use (default: 4)
    - tophits: number of reference sequences to include in alignment phase (default: 20)
    - minscore: minimum score cutoff to consider for classification (default: 30)
    - minid1: minimum number of shared kmers in second phase kmer search (default: 1)
    - minid2: minimum number of shared kmers in alignment phase (default: 1)
    """
    # reference database dir
    reference_dir = join(os.path.dirname(reference_dict['seqs']), 'mapseq')
    if not exists(reference_dir):
        if not dry_run:
            os.makedirs(reference_dir)

    # build reference/input/output paths
    input_sequences = f"{rep_seq_base}.fasta"
    reference_fasta = reference_dict['database_fasta']
    reference_tax = reference_dict['database_tax']
    taxonomy_raw = join(outdir, 'taxonomy_raw.mseq')
    taxonomy_final = join(outdir, 'taxonomy_final.tsv')

    # Check if reference database files exist
    if not reference_fasta or not reference_tax:
        raise FileNotFoundError(f"MAPseq database files not found in {reference_dir}")
    
    if not exists(reference_fasta) or not exists(reference_tax):
        raise FileNotFoundError(f"MAPseq database files missing: {reference_fasta}, {reference_tax}")

    # Extract parameters with defaults
    nthreads = param_dict.get('nthreads', 4)
    tophits = param_dict.get('tophits', 20)
    minscore = param_dict.get('minscore', 30)
    minid1 = param_dict.get('minid1', 1)
    minid2 = param_dict.get('minid2', 1)
    cutoff = param_dict.get('cutoff', 0.0)  # Add cutoff parameter support

    # Build MAPseq command
    mapseq_cmd = (
        f"mapseq "
        f"-nthreads {nthreads} "
        f"-tophits {tophits} "
        f"-minscore {minscore} "
        f"-minid1 {minid1} "
        f"-minid2 {minid2} "
        f"{input_sequences} "
        f"{reference_fasta} "
        f"{reference_tax} "
        f"> {taxonomy_raw}"
    )

    # Assemble full command
    cutoff_arg = f" {cutoff}" if cutoff > 0.0 else ""
    cmd = (
        f"conda run -p {conda_env} bash -c \""
        f"mkdir -p {outdir} && "
        f"{mapseq_cmd} && "
        f"python {scripts_dir}/process_mapseq_output.py "
            f"{taxonomy_raw} {taxonomy_final}{cutoff_arg}\""
    )
    
    return cmd


