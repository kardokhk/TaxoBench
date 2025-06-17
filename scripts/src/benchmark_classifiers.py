#!/usr/bin/env python3
"""
Performance Benchmarking Script for Taxonomic Classifiers

This script measures classification run times and memory usage for 6 different 
taxonomic classifiers (TaxoNova, VSsearch, BLAST, Naive Bayes, IdTaxa, SPINGO)
across different sequence set sizes and amplicon regions.

The script generates test datasets of varying sizes (100-10000 sequences) and
measures performance metrics for each classifier.
"""

import os
import sys
import time
import psutil
import subprocess
import threading
import random
import shutil
import glob
from datetime import datetime
from pathlib import Path
from collections import defaultdict
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import json

# Add scripts directory to path
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(script_dir)

import utils
import generate_commands

class ResourceMonitor:
    """Monitor system resource usage during command execution"""
    
    def __init__(self):
        self.reset()
    
    def reset(self):
        self.max_memory_mb = 0
        self.start_time = None
        self.end_time = None
        self.monitoring = False
        self.process = None
        
    def start_monitoring(self, process):
        """Start monitoring resources for a given process"""
        self.process = process
        self.start_time = time.time()
        self.monitoring = True
        self.max_memory_mb = 0
        
        # Start monitoring thread
        self.monitor_thread = threading.Thread(target=self._monitor_resources)
        self.monitor_thread.daemon = True
        self.monitor_thread.start()
    
    def stop_monitoring(self):
        """Stop monitoring and return results"""
        self.monitoring = False
        self.end_time = time.time()
        
        # Wait for monitoring thread to finish
        if hasattr(self, 'monitor_thread'):
            self.monitor_thread.join(timeout=1.0)
        
        return {
            'runtime_seconds': self.end_time - self.start_time,
            'max_memory_mb': self.max_memory_mb
        }
    
    def _monitor_resources(self):
        """Monitor memory usage in a separate thread"""
        while self.monitoring:
            try:
                if self.process and self.process.poll() is None:
                    # Get memory usage of main process and all children
                    try:
                        parent = psutil.Process(self.process.pid)
                        children = parent.children(recursive=True)
                        
                        total_memory = parent.memory_info().rss
                        for child in children:
                            try:
                                total_memory += child.memory_info().rss
                            except (psutil.NoSuchProcess, psutil.AccessDenied):
                                pass
                        
                        memory_mb = total_memory / (1024 * 1024)
                        self.max_memory_mb = max(self.max_memory_mb, memory_mb)
                        
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
                
                time.sleep(0.1)  # Monitor every 100ms
                
            except Exception as e:
                print(f"Error monitoring resources: {e}")
                break

class PerformanceBenchmarker:
    """Main class for benchmarking taxonomic classifiers"""
    
    def __init__(self, config_file=None):
        self.load_config(config_file)
        self.setup_paths()
        self.results = []
        
    def load_config(self, config_file):
        """Load configuration from file or use defaults"""
        if config_file and os.path.exists(config_file):
            with open(config_file, 'r') as f:
                config = json.load(f)
        else:
            config = {}
            
        # Default configuration
        self.config = {
            'project_dir': config.get('project_dir', '/home/kardokhk/taxonova_bench/taxobench_eight'),
            'sequence_sizes': config.get('sequence_sizes', [100, 500, 1000]),
            'v_regions': config.get('v_regions', ['V4', 'V3-V4', 'V1-V9']),
            'mock_names': config.get('mock_names', ['mockrobiota-16S']),
            'classifiers': config.get('classifiers', ['taxonova', 'vsearch', 'blast', 'naive-bayes', 'idtaxa', 'spingo']),
            'n_replicates': config.get('n_replicates', 3),
            'output_dir': config.get('output_dir', 'benchmark_results'),
            'conda_env': config.get('conda_env', '/home/kardokhk/miniforge3/envs/mock-taxonomy')
        }
        
        # Mock to v_region mapping
        self.mock_v_regions = {
            'gutmock': ['V4', 'V3-V4', 'V1-V9'],
            'vaginalmock': ['V4', 'V3-V4', 'V1-V9'],
            'mockrobiota-16S': ['V4', 'V3-V4', 'V1-V9'],
            'mockrobiota-ITS': ['ITS1', 'ITS2', 'ITS12']
        }
        
    def setup_paths(self):
        """Setup directory paths"""
        self.project_dir = self.config['project_dir']
        self.scripts_dir = os.path.join(self.project_dir, 'scripts', 'src')
        self.reference_database_dir = os.path.join(self.project_dir, 'benchmark', 'data', 'ref_dbs')
        self.mock_dir = os.path.join(self.project_dir, 'benchmark', 'data', 'mock-community')
        self.output_dir = os.path.join(self.project_dir, self.config['output_dir'])
        
        # Create output directory
        os.makedirs(self.output_dir, exist_ok=True)
        
    def generate_test_sequences(self, source_fasta, output_fasta, n_sequences, random_seed=42):
        """Generate test sequence files with specified number of sequences"""
        random.seed(random_seed)
        
        # Read all sequences
        sequences = list(SeqIO.parse(source_fasta, "fasta"))
        
        if len(sequences) < n_sequences:
            # If we need more sequences, sample with replacement
            selected_sequences = random.choices(sequences, k=n_sequences)
        else:
            # Sample without replacement
            selected_sequences = random.sample(sequences, n_sequences)
        
        # Rename sequences to avoid duplicates
        renamed_sequences = []
        for i, seq in enumerate(selected_sequences):
            new_id = f"test_seq_{i+1:06d}"
            new_seq = SeqRecord(seq.seq, id=new_id, description=f"Test sequence {i+1}")
            renamed_sequences.append(new_seq)
        
        # Write to output file
        SeqIO.write(renamed_sequences, output_fasta, "fasta")
        return len(renamed_sequences)
    
    def setup_test_datasets(self):
        """Generate test datasets for all combinations of parameters"""
        test_datasets = {}
        
        for mock_name in self.config['mock_names']:
            test_datasets[mock_name] = {}
            
            for v_region in self.config['v_regions']:
                # Skip if region not applicable to mock
                if v_region not in self.mock_v_regions.get(mock_name, []):
                    continue
                    
                test_datasets[mock_name][v_region] = {}
                
                # Find source sequences
                mockdir = os.path.join(self.mock_dir, mock_name, v_region, 'source')
                rep_seq_files = glob.glob(os.path.join(mockdir, "*", f'rep_seqs_{mock_name}.fasta'))
                
                if not rep_seq_files:
                    print(f"Warning: No rep_seqs file found for {mock_name}/{v_region}")
                    continue
                    
                source_fasta = rep_seq_files[0]
                
                for seq_size in self.config['sequence_sizes']:
                    # Create test dataset directory
                    test_dir = os.path.join(self.output_dir, 'test_datasets', mock_name, v_region, f'size_{seq_size}')
                    os.makedirs(test_dir, exist_ok=True)
                    
                    test_fasta = os.path.join(test_dir, f'test_sequences_{seq_size}.fasta')
                    
                    # Generate test sequences if not exists
                    if not os.path.exists(test_fasta):
                        actual_size = self.generate_test_sequences(source_fasta, test_fasta, seq_size)
                        print(f"Generated {actual_size} test sequences for {mock_name}/{v_region}/size_{seq_size}")
                    
                    test_datasets[mock_name][v_region][seq_size] = test_fasta
        
        return test_datasets
    
    def get_reference_database(self):
        """Load reference database configuration"""
        # This is adapted from the original assign_mock_community.ipynb
        v_regions = self.config['v_regions'] + ['ITS1', 'ITS2', 'ITS12']  # Include all possible regions
        
        common_database = {}
        for v_region in v_regions:
            v_region_dir = os.path.join(self.reference_database_dir, v_region)
            if os.path.isdir(v_region_dir):
                merged_file = os.path.join(v_region_dir, "merged.fasta")
                seq_file = os.path.join(v_region_dir, "sequences.fasta")
                taxa_file = os.path.join(v_region_dir, "taxonomy.tsv")

                if os.path.isfile(seq_file) and os.path.isfile(taxa_file):
                    common_database[v_region] = {
                        'seqs': seq_file,
                        'taxa': taxa_file,
                        'merged': merged_file if os.path.isfile(merged_file) else None
                    }

        # Tool-specific models
        reference_database = {
            'taxonova': {},
            'vsearch': {},
            'blast': {},
            'naive-bayes': {},
            'idtaxa': {},
            'spingo': {}
        }

        for v_region, files in common_database.items():
            v_region_dir = os.path.join(self.reference_database_dir, v_region)
            
            # Find tool-specific models
            taxonova_model = None
            nb_model = None
            idtaxa_model = None
            spingo_database = None

            # Search for models
            taxonova_dir = os.path.join(v_region_dir, 'taxonova')
            if os.path.isdir(taxonova_dir):
                models = glob.glob(os.path.join(taxonova_dir, '*model.pt'))
                if models:
                    taxonova_model = models[0]

            nb_dir = os.path.join(v_region_dir, 'naive-bayes')
            if os.path.isdir(nb_dir):
                nb_models = glob.glob(os.path.join(nb_dir, '*sklearn*qza'))
                if nb_models:
                    nb_model = nb_models[0]

            idtaxa_dir = os.path.join(v_region_dir, 'idtaxa')
            if os.path.isdir(idtaxa_dir):
                idtaxa_models = glob.glob(os.path.join(idtaxa_dir, '*.rds')) + glob.glob(os.path.join(idtaxa_dir, '*.RData'))
                if idtaxa_models:
                    idtaxa_model = idtaxa_models[0]

            spingo_dir = os.path.join(v_region_dir, 'spingo')
            if os.path.isdir(spingo_dir):
                spingo_db_files = glob.glob(os.path.join(spingo_dir, '*.fasta'))
                if spingo_db_files:
                    spingo_database = spingo_db_files[0]

            # Populate reference_database
            reference_database['taxonova'][v_region] = {
                'seqs': files['seqs'],
                'taxa': files['taxa'],
                'merged': files['merged'],
                'model': taxonova_model
            }
            reference_database['naive-bayes'][v_region] = {
                'seqs': files['seqs'],
                'taxa': files['taxa'],
                'model': nb_model
            }
            reference_database['vsearch'][v_region] = {
                'seqs': files['seqs'],
                'taxa': files['taxa']
            }
            reference_database['blast'][v_region] = {
                'seqs': files['seqs'],
                'taxa': files['taxa']
            }
            reference_database['idtaxa'][v_region] = {
                'seqs': files['seqs'],
                'taxa': files['taxa'],
                'model': idtaxa_model
            }
            reference_database['spingo'][v_region] = {
                'seqs': files['seqs'],
                'taxa': files['taxa'],
                'database': spingo_database
            }

        return reference_database
    
    def run_classifier_benchmark(self, classifier, test_fasta, v_region, reference_database, replicate=1):
        """Run benchmark for a specific classifier"""
        
        # Create unique output directory for this benchmark
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        bench_dir = os.path.join(self.output_dir, 'benchmark_runs', 
                                 f'{classifier}_{v_region}_{os.path.basename(test_fasta).replace(".fasta", "")}_{replicate}_{timestamp}')
        os.makedirs(bench_dir, exist_ok=True)
        
        # Default parameters for benchmarking (using single parameter set for speed)
        default_params = {
            'taxonova': {'species_threshold': 0.7},
            'vsearch': {'p-maxaccepts': 10, 'p-perc-identity': 0.97, 'p-min-consensus': 0.51},
            'blast': {'p-maxaccepts': 10, 'p-perc-identity': 0.97, 'p-min-consensus': 0.51},
            'naive-bayes': {'p-confidence': 0.7},
            'idtaxa': {'p-confidence': 50},
            'spingo': {'kmersize': 8, 'bootstrap': 50, 'subsample': 10}
        }
        
        param_dict = default_params.get(classifier, {})
        
        # Build command based on classifier
        cmd = None
        conda_env = self.config['conda_env']
        #n_cpus = os.cpu_count() - 1
        n_cpus = 30
        
        # Convert test_fasta to base name (without extension) for compatibility
        test_base = os.path.splitext(test_fasta)[0]
        
        try:
            if classifier == 'taxonova':
                cmd = generate_commands.build_taxonova_command(
                    rep_seq_base=test_base,
                    outdir=bench_dir,
                    param_dict=param_dict,
                    reference_dict=reference_database['taxonova'][v_region],
                    scripts_dir=self.scripts_dir,
                    conda_env=conda_env,
                    dry_run=False
                )
            elif classifier == 'vsearch':
                cmd = generate_commands.build_vsearch_command(
                    rep_seq_base=test_base,
                    outdir=bench_dir,
                    param_dict=param_dict,
                    reference_dict=reference_database['vsearch'][v_region],
                    conda_env=conda_env,
                    n_threads=n_cpus,
                    dry_run=False
                )
            elif classifier == 'blast':
                cmd = generate_commands.build_blast_command(
                    rep_seq_base=test_base,
                    outdir=bench_dir,
                    param_dict=param_dict,
                    reference_dict=reference_database['blast'][v_region],
                    conda_env=conda_env,
                    n_threads=n_cpus,
                    dry_run=False
                )
            elif classifier == 'naive-bayes':
                cmd = generate_commands.build_naive_bayes_command(
                    rep_seq_base=test_base,
                    outdir=bench_dir,
                    param_dict=param_dict,
                    reference_dict=reference_database['naive-bayes'][v_region],
                    conda_env=conda_env,
                    n_threads=n_cpus,
                    dry_run=False
                )
            elif classifier == 'idtaxa':
                cmd = generate_commands.build_idtaxa_command(
                    rep_seq_base=test_base,
                    outdir=bench_dir,
                    param_dict=param_dict,
                    reference_dict=reference_database['idtaxa'][v_region],
                    scripts_dir=self.scripts_dir,
                    conda_env=conda_env,
                    dry_run=False
                )
            elif classifier == 'spingo':
                spingo_conda_config = {
                    'spingo_env': '/home/kardokhk/miniforge3/envs/spingo-env',
                    'processing_env': conda_env
                }
                cmd = generate_commands.build_spingo_command(
                    rep_seq_base=test_base,
                    outdir=bench_dir,
                    param_dict=param_dict,
                    reference_dict=reference_database['spingo'][v_region],
                    scripts_dir=self.scripts_dir,
                    conda_env=spingo_conda_config,
                    n_threads=n_cpus,
                    dry_run=False
                )
                
        except Exception as e:
            print(f"Error building command for {classifier}: {e}")
            return None
            
        if not cmd:
            print(f"Failed to build command for {classifier}")
            return None
            
        # Run command with resource monitoring
        print(f"Running {classifier} benchmark...")
        print(f"Command: {cmd}")
        
        monitor = ResourceMonitor()
        
        try:
            # Start process
            process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, 
                                       stderr=subprocess.PIPE, universal_newlines=True)
            
            # Start monitoring
            monitor.start_monitoring(process)
            
            # Wait for completion
            stdout, stderr = process.communicate()
            
            # Stop monitoring and get results
            performance_metrics = monitor.stop_monitoring()
            
            # Check if successful
            success = process.returncode == 0
            
            if not success:
                print(f"Command failed with return code {process.returncode}")
                print(f"Error: {stderr}")
            
            return {
                'classifier': classifier,
                'v_region': v_region,
                'test_file': test_fasta,
                'replicate': replicate,
                'success': success,
                'runtime_seconds': performance_metrics['runtime_seconds'],
                'max_memory_mb': performance_metrics['max_memory_mb'],
                'output_dir': bench_dir,
                'stdout': stdout,
                'stderr': stderr,
                'timestamp': timestamp
            }
            
        except Exception as e:
            print(f"Error running benchmark for {classifier}: {e}")
            return None
    
    def run_benchmarks(self):
        """Run all benchmarking combinations"""
        # Setup test datasets
        print("Setting up test datasets...")
        test_datasets = self.setup_test_datasets()
        
        # Get reference databases
        print("Loading reference databases...")
        reference_database = self.get_reference_database()
        
        # Run benchmarks
        total_runs = sum(
            len(self.config['classifiers']) * len(sizes) * self.config['n_replicates']
            for mock_data in test_datasets.values()
            for sizes in mock_data.values()
        )
        
        current_run = 0
        
        for mock_name, mock_data in test_datasets.items():
            for v_region, region_data in mock_data.items():
                for seq_size, test_fasta in region_data.items():
                    
                    # Get number of sequences in test file
                    n_sequences = len(list(SeqIO.parse(test_fasta, "fasta")))
                    
                    for classifier in self.config['classifiers']:
                        # Check if classifier has reference data for this region
                        if v_region not in reference_database.get(classifier, {}):
                            print(f"Skipping {classifier} for {v_region} - no reference data")
                            continue
                            
                        for replicate in range(1, self.config['n_replicates'] + 1):
                            current_run += 1
                            print(f"\nProgress: {current_run}/{total_runs}")
                            print(f"Running: {classifier} | {mock_name} | {v_region} | {n_sequences} seqs | replicate {replicate}")
                            
                            result = self.run_classifier_benchmark(
                                classifier=classifier,
                                test_fasta=test_fasta,
                                v_region=v_region,
                                reference_database=reference_database,
                                replicate=replicate
                            )
                            
                            if result:
                                # Add additional metadata
                                result.update({
                                    'mock_name': mock_name,
                                    'sequence_count': n_sequences,
                                    'sequence_size_category': seq_size
                                })
                                self.results.append(result)
                            
                            # Save intermediate results
                            self.save_results()
    
    def save_results(self):
        """Save benchmark results to files"""
        if not self.results:
            return
            
        # Convert to DataFrame
        df = pd.DataFrame(self.results)
        
        # Save detailed results
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        detailed_file = os.path.join(self.output_dir, f'detailed_benchmark_results_{timestamp}.csv')
        df.to_csv(detailed_file, index=False)
        
        # Save summary results
        summary_df = df[df['success'] == True].groupby([
            'classifier', 'v_region', 'mock_name', 'sequence_count'
        ]).agg({
            'runtime_seconds': ['mean', 'std', 'min', 'max'],
            'max_memory_mb': ['mean', 'std', 'min', 'max']
        }).round(3)
        
        summary_df.columns = ['_'.join(col).strip() for col in summary_df.columns.values]
        summary_file = os.path.join(self.output_dir, f'summary_benchmark_results_{timestamp}.csv')
        summary_df.to_csv(summary_file)
        
        print(f"Results saved to:")
        print(f"  Detailed: {detailed_file}")
        print(f"  Summary: {summary_file}")
        
        return detailed_file, summary_file

def main():
    """Main function to run benchmarking"""
    # Initialize benchmarker
    benchmarker = PerformanceBenchmarker()
    
    print("Starting Performance Benchmarking...")
    print(f"Configuration:")
    print(f"  Project directory: {benchmarker.config['project_dir']}")
    print(f"  Sequence sizes: {benchmarker.config['sequence_sizes']}")
    print(f"  V-regions: {benchmarker.config['v_regions']}")
    print(f"  Classifiers: {benchmarker.config['classifiers']}")
    print(f"  Replicates: {benchmarker.config['n_replicates']}")
    print(f"  Output directory: {benchmarker.output_dir}")
    
    # Run benchmarks
    benchmarker.run_benchmarks()
    
    # Generate final report
    print("\nGenerating final report...")
    benchmarker.save_results()
    
    print("Benchmarking completed!")

if __name__ == "__main__":
    main()
