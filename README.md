# TaxoBench

A comprehensive benchmarking framework for taxonomic classification tools.

## Overview

TaxoBench provides tools and scripts for evaluating the performance of various taxonomic classification methods on mock microbial communities.

## Structure

- `benchmark/` - Benchmark datasets and test files
- `scripts/` - Main analysis scripts and notebooks
  - `src/` - Core Python and R modules for data processing and evaluation
  - `*.ipynb` - Jupyter notebooks for interactive analysis
  - `*.py` - Python scripts for automated benchmarking

## Features

- Mock community generation and evaluation
- Support for multiple taxonomic classifiers
- Comprehensive performance metrics
- Runtime and memory usage benchmarking

## Getting Started

1. Clone this repository
2. Install required dependencies
3. Run the notebooks in the `scripts/` directory to get started

## Scripts

- `generate_mock_community.ipynb` - Generate synthetic mock communities
- `prepare_mock_community.ipynb` - Prepare mock communities for analysis
- `assign_mock_community.ipynb` - Assign taxonomy to mock communities
- `evaluate_mock_community.ipynb` - Evaluate classification performance
- `runtime_usage_benchmark.py` - Performance benchmarking

## Contributing

Contributions are welcome! Please feel free to submit issues and pull requests.

## License

This project is open source.
