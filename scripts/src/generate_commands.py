import os
import subprocess
from os.path import join, exists
from qiime2 import Artifact
from qiime2.plugins.feature_classifier.methods import makeblastdb
from qiime2.plugins.rescript.pipelines import evaluate_fit_classifier



def build_taxonova_command(
    rep_seq_base, outdir, param_dict, reference_dict, 
    scripts_dir, conda_env, dry_run
):
    # reference database dir
    reference_dir = join(os.path.dirname(reference_dict['seqs']), 'taxonova')
    if not exists(reference_dir):
        if not dry_run:
            os.makedirs(reference_dir)

    # build reference/input/output paths
    input_sequences = f"{rep_seq_base}.fasta"
    reference_seqs = reference_dict['merged']
    model_file = reference_dict['model']
    taxonomy_raw = join(outdir, 'taxonomy_raw.tsv')
    taxonomy_final = join(outdir, 'taxonomy_final.tsv')

    # train on reference sequence if model is not available
    if model_file is None:
        if not dry_run:
            # save final qza file
            cmd = ["taxonova", "training", "--verbose",
                   "--full_data", 
                   "--train_db", reference_seqs, 
                   "--model_output", model_file]
            subprocess.run(cmd, check=True)

    # Build threshold arguments from param_dict
    threshold_args = ""
    threshold_params = [
        'species_threshold', 'genus_threshold', 'family_threshold',
        'order_threshold', 'class_threshold', 'phylum_threshold'
    ]
    
    for param in threshold_params:
        if param in param_dict:
            threshold_args += f"--{param} {param_dict[param]} "

    # assemble command
    cmd = (
        f"conda run -p {conda_env} bash -c \""
        f"mkdir -p {outdir} && "
        f"taxonova classification "
            f"{threshold_args}"
            f"--input {input_sequences} "
            f"--model_input {model_file} "
            f"--output {taxonomy_raw} && "
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
    input_sequences = f"{rep_seq_base}.qza"
    reference_seqs = join(reference_dir, "sequences.qza")
    reference_taxa = join(reference_dir, "taxonomy.qza")
    classification_out = join(outdir, 'taxonomy_classification_out.qza')
    search_out = join(outdir, 'taxonomy_search_out.qza')
    taxonomy_final = join(outdir, 'taxonomy_final.tsv')

    # Import reference sequences and taxonomy
    if not all(os.path.exists(f) for f in [reference_seqs, reference_taxa]):
        if not dry_run:
            # Import seqs
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
    if not os.path.exists(input_sequences):
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
    input_sequences = f"{rep_seq_base}.qza"
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
    if not os.path.exists(input_sequences):
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
    input_sequences = f"{rep_seq_base}.qza"
    model_file = reference_dict['model']
    evaluation_file = join(reference_dir, "evaluation.qza")
    prediction_file = join(reference_dir, "predictions.qza")
    classification_out = join(outdir, 'taxonomy_classification_out.qza')
    taxonomy_final = join(outdir, 'taxonomy_final.tsv')

    # Import reference sequences and taxonomy and build classifier
    if model_file is None:
        if not dry_run:
            # Import seqs
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
    if not os.path.exists(input_sequences):
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
