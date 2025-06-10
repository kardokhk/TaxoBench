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
        # define model_file
        model_file = join(reference_dir, 'greengens2_taxonova_model.pt')
        if not dry_run:
            if not exists(model_file):
                # train a model
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

    # assemble command
    cmd = (
            f'conda run -p {conda_env} bash -c \"'
            f'mkdir -p {outdir} && '
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
            f'python -c \\"import pandas as pd; '
                f'from qiime2 import Artifact; '
                f'df = Artifact.load(\'{classification_out}\').view(pd.DataFrame); '
                f'df.to_csv(\'{taxonomy_final}\', sep=\'\\\\t\', index=True)\\"'
            f'"'
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

    # assemble command
    cmd = (
            f'conda run -p {conda_env} bash -c \"'
            f'mkdir -p {outdir} && '
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
            f'python -c \\"import pandas as pd; '
                f'from qiime2 import Artifact; '
                f'df = Artifact.load(\'{classification_out}\').view(pd.DataFrame); '
                f'df.to_csv(\'{taxonomy_final}\', sep=\'\\\\t\', index=True)\\"'
            f'"'
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
        # define model_file
        model_file = join(reference_dir, 'greengens2_sklearn_classifier.qza')
        if not dry_run:
            if not exists(model_file):
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
                    confidence = 0.70
                    )
                )

                # Save outputs
                evaluation.save(evaluation_file)
                observed_taxonomy.save(prediction_file)
                classifier.save(model_file)

    # assemble command
    cmd = (
            f'conda run -p {conda_env} bash -c \"'
            f'mkdir -p {outdir} && '
            f'qiime feature-classifier classify-sklearn '
                f'--i-reads {input_sequences} '
                f'--i-classifier {model_file} '
                f'--p-n-jobs {n_threads} '
                f'--p-pre-dispatch "2*n_jobs" '
                f'--p-confidence {param_dict["p-confidence"]} '
                f'--o-classification {classification_out} && '
            f'python -c \\"import pandas as pd; '
                f'from qiime2 import Artifact; '
                f'df = Artifact.load(\'{classification_out}\').view(pd.DataFrame); '
                f'df.to_csv(\'{taxonomy_final}\', sep=\'\\\\t\', index=True)\\"'
            f'"'
    )
        
    return cmd



def build_IdTaxa_command(
    rep_seq_base, outdir, param_dict, reference_dict, 
    scripts_dir, conda_env, n_threads, dry_run
):
    # reference database dir
    reference_dir = join(os.path.dirname(reference_dict['seqs']), 'IdTaxa')
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
        # define model_file
        model_file = join(reference_dir, 'greengens2_IdTaxa_classifier.RDS')
        if not dry_run:
            if not exists(model_file):
                # train a model
                cmd = (
                    f"conda run -p {conda_env} bash -c \""
                    f"Rscript {scripts_dir}/train_with_idTaxa.R "
                        f"{reference_seqs} {model_file} \""
                )
                subprocess.run(cmd, shell=True, check=True)

    # assemble command
    cmd = (
        f"conda run -p {conda_env} bash -c \""
        f"mkdir -p {outdir} && "
        f"Rscript {scripts_dir}/classify_with_idTaxa.R "
            f"{model_file} {input_sequences} {taxonomy_raw} "
            f"{param_dict['threshold']} {n_threads} && "
        f"python {scripts_dir}/process_taxonova_output.py "
            f"{taxonomy_raw} {taxonomy_final} False \""
    )

    return cmd