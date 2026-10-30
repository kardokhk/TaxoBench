#!/usr/bin/env python3
"""Runtime / memory benchmarking for the 8 supported classifiers.

Supported classifiers:
    taxonova, vsearch, blast, naive-bayes, idtaxa, spingo, kraken2, mapseq

CSV-only output. No plotting. No removed classifiers (rdp / mothur / emu / qiime-OTU).

This module is a library; orchestration lives in
`scripts/06_runtime.py` which builds a config and calls
`PerformanceBenchmarker`.
"""

import os
import sys
import time
import json
import glob
import random
import threading
import subprocess
from datetime import datetime

import psutil
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Make sibling modules importable
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_commands  # noqa: E402


SUPPORTED_CLASSIFIERS = (
    'taxonova', 'vsearch', 'blast', 'naive-bayes',
    'idtaxa', 'spingo', 'kraken2', 'mapseq',
)


class ResourceMonitor:
    """Sample peak RSS of a subprocess + its children every 100 ms."""

    def __init__(self):
        self.reset()

    def reset(self):
        self.max_memory_mb = 0
        self.start_time = None
        self.end_time = None
        self.monitoring = False
        self.process = None

    def start_monitoring(self, process):
        self.process = process
        self.start_time = time.time()
        self.monitoring = True
        self.max_memory_mb = 0
        self.monitor_thread = threading.Thread(target=self._monitor_resources)
        self.monitor_thread.daemon = True
        self.monitor_thread.start()

    def stop_monitoring(self):
        self.monitoring = False
        self.end_time = time.time()
        if hasattr(self, 'monitor_thread'):
            self.monitor_thread.join(timeout=1.0)
        return {
            'runtime_seconds': self.end_time - self.start_time,
            'max_memory_mb': self.max_memory_mb,
        }

    def _monitor_resources(self):
        while self.monitoring:
            try:
                if self.process and self.process.poll() is None:
                    try:
                        parent = psutil.Process(self.process.pid)
                        children = parent.children(recursive=True)
                        total = parent.memory_info().rss
                        for c in children:
                            try:
                                total += c.memory_info().rss
                            except (psutil.NoSuchProcess, psutil.AccessDenied):
                                pass
                        self.max_memory_mb = max(self.max_memory_mb, total / (1024 * 1024))
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
                time.sleep(0.1)
            except Exception as exc:
                print(f"[ResourceMonitor] error: {exc}")
                break


class PerformanceBenchmarker:
    """Build runtime/memory benchmarks for the 8 supported classifiers."""

    def __init__(self, config):
        if isinstance(config, str):
            with open(config) as fh:
                config = json.load(fh)
        self.config = config
        self._validate_config()
        self.results = []

    # ------------------------------------------------------------------
    def _validate_config(self):
        required = ('project_dir', 'reference_database_dir',
                    'output_dir', 'classifiers', 'sequence_sizes',
                    'v_regions', 'n_replicates', 'conda_envs', 'source_fasta')
        for key in required:
            if key not in self.config:
                raise ValueError(f"Missing config key: {key}")
        bad = [c for c in self.config['classifiers'] if c not in SUPPORTED_CLASSIFIERS]
        if bad:
            raise ValueError(f"Unsupported classifiers: {bad}. Allowed: {SUPPORTED_CLASSIFIERS}")
        os.makedirs(self.config['output_dir'], exist_ok=True)

    # ------------------------------------------------------------------
    def generate_test_sequences(self, source_fasta, output_fasta, n_sequences,
                                random_seed=42):
        random.seed(random_seed)
        seqs = list(SeqIO.parse(source_fasta, "fasta"))
        if len(seqs) < n_sequences:
            selected = random.choices(seqs, k=n_sequences)
        else:
            selected = random.sample(seqs, n_sequences)
        renamed = [SeqRecord(s.seq, id=f"test_seq_{i+1:06d}", description=f"Test sequence {i+1}")
                   for i, s in enumerate(selected)]
        SeqIO.write(renamed, output_fasta, "fasta")
        return len(renamed)

    # ------------------------------------------------------------------
    def setup_test_datasets(self):
        """Create one fasta per (region x size) under output_dir/test_datasets."""
        test_datasets = {}
        for v_region in self.config['v_regions']:
            test_datasets[v_region] = {}
            for size in self.config['sequence_sizes']:
                tdir = os.path.join(self.config['output_dir'], 'test_datasets',
                                    v_region, f'size_{size}')
                os.makedirs(tdir, exist_ok=True)
                fasta = os.path.join(tdir, f'test_sequences_{size}.fasta')
                if not os.path.exists(fasta):
                    self.generate_test_sequences(self.config['source_fasta'], fasta, size)
                test_datasets[v_region][size] = fasta
        return test_datasets

    # ------------------------------------------------------------------
    def get_reference_database(self):
        """Discover reference DB files for each (classifier x region)."""
        ref_root = self.config['reference_database_dir']
        common = {}
        for v_region in self.config['v_regions']:
            v_dir = os.path.join(ref_root, v_region)
            if not os.path.isdir(v_dir):
                continue
            if v_region == 'RRN':
                merged = os.path.join(v_dir, "grond_clean.fasta")
                seq    = os.path.join(v_dir, "grond_clean.fasta")
                taxa   = os.path.join(v_dir, "grond_clean.tsv")
            else:
                merged = os.path.join(v_dir, "merged.fasta")
                seq    = os.path.join(v_dir, "sequences.fasta")
                taxa   = os.path.join(v_dir, "taxonomy.tsv")
            if os.path.isfile(seq) and os.path.isfile(taxa):
                common[v_region] = {
                    'seqs': seq, 'taxa': taxa,
                    'merged': merged if os.path.isfile(merged) else None,
                }

        ref_db = {c: {} for c in SUPPORTED_CLASSIFIERS}
        for v_region, files in common.items():
            v_dir = os.path.join(ref_root, v_region)

            taxonova_model = generate_commands.resolve_taxonova_model(
                v_region, v_dir, self.config.get('taxonova_settings') or {},
                ref_root=ref_root)
            nb_model       = _first_glob(os.path.join(v_dir, 'naive-bayes', '*sklearn*qza'))
            idtaxa_model   = (_first_glob(os.path.join(v_dir, 'idtaxa', '*.rds'))
                              or _first_glob(os.path.join(v_dir, 'idtaxa', '*.RData')))
            spingo_db      = _first_glob(os.path.join(v_dir, 'spingo', '*.fasta'))

            kraken2_db = None
            kraken2_dir = os.path.join(v_dir, 'kraken2')
            if os.path.isdir(kraken2_dir) and all(
                    os.path.exists(os.path.join(kraken2_dir, f))
                    for f in ('hash.k2d', 'opts.k2d', 'taxo.k2d')):
                kraken2_db = kraken2_dir

            mapseq_fa = mapseq_tx = None
            mapseq_dir = os.path.join(v_dir, 'mapseq')
            if os.path.isdir(mapseq_dir):
                fa = os.path.join(mapseq_dir, 'reference.fasta')
                tx = os.path.join(mapseq_dir, 'reference.tax')
                if os.path.exists(fa) and os.path.exists(tx):
                    mapseq_fa, mapseq_tx = fa, tx

            ref_db['taxonova'][v_region] = {
                'seqs': files['seqs'], 'taxa': files['taxa'],
                'merged': files['merged'], 'model': taxonova_model,
            }
            ref_db['naive-bayes'][v_region] = {
                'seqs': files['seqs'], 'taxa': files['taxa'], 'model': nb_model,
            }
            ref_db['vsearch'][v_region] = {'seqs': files['seqs'], 'taxa': files['taxa']}
            ref_db['blast'][v_region]   = {'seqs': files['seqs'], 'taxa': files['taxa']}
            ref_db['idtaxa'][v_region]  = {'seqs': files['seqs'], 'taxa': files['taxa'],
                                           'model': idtaxa_model}
            ref_db['spingo'][v_region]  = {'seqs': files['seqs'], 'taxa': files['taxa'],
                                           'database': spingo_db}
            ref_db['kraken2'][v_region] = {'seqs': files['seqs'], 'taxa': files['taxa'],
                                           'database': kraken2_db}
            ref_db['mapseq'][v_region]  = {'seqs': files['seqs'], 'taxa': files['taxa'],
                                           'database_fasta': mapseq_fa,
                                           'database_tax': mapseq_tx}
        return ref_db

    # ------------------------------------------------------------------
    def _build_command(self, classifier, test_base, bench_dir, v_region,
                       reference_database, n_threads):
        params = self.config.get('classifier_params', {}).get(classifier, {})
        envs = self.config['conda_envs']
        scripts_dir = self.config.get('scripts_dir', os.path.dirname(os.path.abspath(__file__)))

        if classifier == 'taxonova':
            return generate_commands.build_taxonova_command(
                rep_seq_base=test_base, outdir=bench_dir, param_dict=params,
                reference_dict=reference_database['taxonova'][v_region],
                scripts_dir=scripts_dir, conda_env=envs['taxonova'], dry_run=False,
                settings=self.config.get('taxonova_settings'), n_threads=n_threads)
        if classifier == 'vsearch':
            return generate_commands.build_vsearch_command(
                rep_seq_base=test_base, outdir=bench_dir, param_dict=params,
                reference_dict=reference_database['vsearch'][v_region],
                conda_env=envs['vsearch'], n_threads=n_threads, dry_run=False)
        if classifier == 'blast':
            return generate_commands.build_blast_command(
                rep_seq_base=test_base, outdir=bench_dir, param_dict=params,
                reference_dict=reference_database['blast'][v_region],
                conda_env=envs['blast'], n_threads=n_threads, dry_run=False)
        if classifier == 'naive-bayes':
            return generate_commands.build_naive_bayes_command(
                rep_seq_base=test_base, outdir=bench_dir, param_dict=params,
                reference_dict=reference_database['naive-bayes'][v_region],
                conda_env=envs['naive-bayes'], n_threads=n_threads, dry_run=False)
        if classifier == 'idtaxa':
            return generate_commands.build_idtaxa_command(
                rep_seq_base=test_base, outdir=bench_dir, param_dict=params,
                reference_dict=reference_database['idtaxa'][v_region],
                scripts_dir=scripts_dir, conda_env=envs['idtaxa'], dry_run=False)
        if classifier == 'spingo':
            spingo_env = envs['spingo']
            if isinstance(spingo_env, str):
                spingo_env = {'spingo_env': spingo_env,
                              'processing_env': envs.get('taxonova', spingo_env)}
            return generate_commands.build_spingo_command(
                rep_seq_base=test_base, outdir=bench_dir, param_dict=params,
                reference_dict=reference_database['spingo'][v_region],
                scripts_dir=scripts_dir, conda_env=spingo_env,
                n_threads=n_threads, dry_run=False)
        if classifier == 'kraken2':
            return generate_commands.build_kraken2_command(
                rep_seq_base=test_base, outdir=bench_dir, param_dict=params,
                reference_dict=reference_database['kraken2'][v_region],
                scripts_dir=scripts_dir, conda_env=envs['kraken2'],
                n_threads=n_threads, dry_run=False)
        if classifier == 'mapseq':
            return generate_commands.build_mapseq_command(
                rep_seq_base=test_base, outdir=bench_dir, param_dict=params,
                reference_dict=reference_database['mapseq'][v_region],
                scripts_dir=scripts_dir, conda_env=envs['mapseq'],
                n_threads=n_threads, dry_run=False)
        raise ValueError(f"Unsupported classifier: {classifier}")

    # ------------------------------------------------------------------
    def run_classifier_benchmark(self, classifier, test_fasta, v_region,
                                 reference_database, replicate=1):
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        bench_dir = os.path.join(self.config['output_dir'], 'benchmark_runs',
                                 f"{classifier}_{v_region}_"
                                 f"{os.path.basename(test_fasta).replace('.fasta','')}_"
                                 f"r{replicate}_{ts}")
        os.makedirs(bench_dir, exist_ok=True)

        test_base = os.path.splitext(test_fasta)[0]
        n_threads = self.config.get('n_threads', max(1, (os.cpu_count() or 2) - 1))

        try:
            cmd = self._build_command(classifier, test_base, bench_dir, v_region,
                                      reference_database, n_threads)
        except Exception as exc:
            print(f"[build] {classifier} failed: {exc}")
            return None

        print(f"[run] {classifier} | {v_region} | rep {replicate}")
        monitor = ResourceMonitor()
        try:
            proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, universal_newlines=True)
            monitor.start_monitoring(proc)
            stdout, stderr = proc.communicate()
            metrics = monitor.stop_monitoring()
            success = proc.returncode == 0
            if not success:
                print(f"  return-code {proc.returncode}: {stderr[:200]}")
            return {
                'classifier': classifier, 'v_region': v_region,
                'test_file': test_fasta, 'replicate': replicate,
                'success': success,
                'runtime_seconds': metrics['runtime_seconds'],
                'max_memory_mb': metrics['max_memory_mb'],
                'output_dir': bench_dir,
                'timestamp': ts,
            }
        except Exception as exc:
            print(f"[run] {classifier} failed: {exc}")
            return None

    # ------------------------------------------------------------------
    def run_benchmarks(self):
        test_datasets = self.setup_test_datasets()
        ref_db = self.get_reference_database()
        for v_region, sizes in test_datasets.items():
            for size, fasta in sizes.items():
                n = len(list(SeqIO.parse(fasta, "fasta")))
                for classifier in self.config['classifiers']:
                    if v_region not in ref_db.get(classifier, {}):
                        print(f"  skip: {classifier} has no ref for {v_region}")
                        continue
                    for rep in range(1, self.config['n_replicates'] + 1):
                        r = self.run_classifier_benchmark(
                            classifier, fasta, v_region, ref_db, replicate=rep)
                        if r:
                            r.update({'sequence_count': n,
                                      'sequence_size_category': size})
                            self.results.append(r)
                        self.save_results()

    # ------------------------------------------------------------------
    def save_results(self):
        if not self.results:
            return None
        df = pd.DataFrame(self.results)
        out_dir = self.config['output_dir']
        # one CSV per classifier
        for classifier, sub in df.groupby('classifier'):
            sub_path = os.path.join(out_dir, f"{classifier}_runtime.csv")
            sub.rename(columns={
                'sequence_count': 'n_seqs',
                'runtime_seconds': 'wall_seconds',
                'max_memory_mb': 'peak_rss_mb',
            })[['n_seqs', 'replicate', 'wall_seconds', 'peak_rss_mb', 'v_region']
               ].rename(columns={'v_region': 'region'}
                        ).to_csv(sub_path, index=False)
        combined = os.path.join(out_dir, "all_runtime_results.csv")
        df.to_csv(combined, index=False)
        return combined


# ---------------------------------------------------------------------------
def _first_glob(pattern):
    matches = sorted(glob.glob(pattern))
    return matches[0] if matches else None


def main():
    import argparse
    ap = argparse.ArgumentParser(description="Runtime/memory benchmark for the 8 classifiers")
    ap.add_argument('--config', required=True, help="JSON config file")
    args = ap.parse_args()
    bench = PerformanceBenchmarker(args.config)
    bench.run_benchmarks()
    bench.save_results()


if __name__ == '__main__':
    main()
