#!/usr/bin/env python3
"""
Command-line script to run the taxonomic classifier benchmarking

This script provides a simple command-line interface to the benchmarking
framework without requiring Jupyter notebook.

Usage:
    python run_benchmark.py --config benchmark_config.json
    python run_benchmark.py --quick-test  # Run with minimal test data
"""

import argparse
import sys
import os
from pathlib import Path

# Add src directory to path
script_dir = Path(__file__).parent
sys.path.append(str(script_dir / "src"))

from benchmark_classifiers import PerformanceBenchmarker
from analyze_benchmark_results import BenchmarkAnalyzer

def main():
    parser = argparse.ArgumentParser(description='Run taxonomic classifier benchmarking')
    parser.add_argument('--config', 
                       help='Path to configuration JSON file',
                       default='benchmark_config.json')
    parser.add_argument('--quick-test', 
                       action='store_true',
                       help='Run a quick test with minimal data')
    parser.add_argument('--analyze-only', 
                       help='Only analyze existing results (provide results directory)')
    parser.add_argument('--output-dir', 
                       help='Output directory for results',
                       default='benchmark_results')
    
    args = parser.parse_args()
    
    if args.analyze_only:
        print(f"Analyzing existing results from: {args.analyze_only}")
        analyzer = BenchmarkAnalyzer(args.analyze_only)
        analyzer.run_complete_analysis()
        return
    
    # Setup quick test configuration if requested
    if args.quick_test:
        print("Running quick test with minimal configuration...")
        config = {
            'project_dir': "/home/kardokhk/taxonova_bench/taxobench_eight",
            'sequence_sizes': [100, 500],  # Minimal sizes
            'v_regions': ['V4'],  # Single region
            'mock_names': ['mockrobiota-16S'],
            'classifiers': ['taxonova', 'vsearch'],  # Just 2 classifiers
            'n_replicates': 1,  # Single replicate
            'output_dir': args.output_dir,
            'conda_env': "/home/kardokhk/miniforge3/envs/mock-taxonomy"
        }
        
        # Save quick test config
        import json
        with open('quick_test_config.json', 'w') as f:
            json.dump(config, f, indent=2)
        args.config = 'quick_test_config.json'
    
    # Run benchmarking
    print(f"Starting benchmarking with config: {args.config}")
    
    try:
        benchmarker = PerformanceBenchmarker(args.config)
        print("Running benchmarks...")
        benchmarker.run_benchmarks()
        
        print("Generating analysis...")
        analyzer = BenchmarkAnalyzer(benchmarker.output_dir)
        analyzer.run_complete_analysis()
        
        print(f"Benchmarking completed! Results in: {benchmarker.output_dir}")
        
    except Exception as e:
        print(f"Error during benchmarking: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
