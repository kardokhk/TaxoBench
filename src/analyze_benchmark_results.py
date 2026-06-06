#!/usr/bin/env python3
"""Aggregate runtime/memory benchmark CSVs into summary tables.

CSV-only output (no plots). Reads either the per-classifier CSVs produced
by `benchmark_classifiers.PerformanceBenchmarker.save_results` or the
combined `all_runtime_results.csv`.
"""

import os
import glob
import json
import argparse

import pandas as pd


class BenchmarkAnalyzer:
    """Build CSV summaries from runtime benchmark results."""

    def __init__(self, results_dir):
        self.results_dir = results_dir
        self.output_dir = os.path.join(results_dir, 'analysis')
        os.makedirs(self.output_dir, exist_ok=True)

    # ------------------------------------------------------------------
    def load_results(self, results_file=None):
        if results_file is None:
            for pattern in ('all_runtime_results.csv',
                            '*_runtime.csv',
                            '*benchmark*.csv'):
                hits = sorted(glob.glob(os.path.join(self.results_dir, pattern)))
                if hits:
                    results_file = hits[-1]
                    break
        if not results_file or not os.path.isfile(results_file):
            raise FileNotFoundError(f"No benchmark CSV found in {self.results_dir}")
        print(f"[load] {results_file}")
        df = pd.read_csv(results_file)
        if 'success' in df.columns:
            df = df[df['success'].astype(bool)]
        if df.empty:
            raise ValueError("No successful benchmark rows.")
        # Normalise column names to runtime spec
        rename = {'runtime_seconds': 'wall_seconds',
                  'max_memory_mb':   'peak_rss_mb',
                  'sequence_count':  'n_seqs',
                  'v_region':        'region'}
        df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
        self.df = df
        return df

    # ------------------------------------------------------------------
    def summarize(self):
        if not hasattr(self, 'df'):
            raise ValueError("Call load_results() first.")
        df = self.df

        meta = {
            'total_runs': int(len(df)),
            'classifiers': sorted(df['classifier'].unique().tolist()),
            'regions': sorted(df['region'].unique().tolist()) if 'region' in df else [],
            'n_seqs': sorted(df['n_seqs'].unique().tolist()) if 'n_seqs' in df else [],
        }
        with open(os.path.join(self.output_dir, 'summary_meta.json'), 'w') as fh:
            json.dump(meta, fh, indent=2)

        agg_cols = {}
        if 'wall_seconds' in df:
            agg_cols['wall_seconds'] = ['mean', 'std', 'min', 'max']
        if 'peak_rss_mb' in df:
            agg_cols['peak_rss_mb'] = ['mean', 'std', 'min', 'max']

        group_cols = [c for c in ('classifier', 'region', 'n_seqs') if c in df.columns]
        summary = df.groupby(group_cols).agg(agg_cols).round(3)
        summary.columns = ['_'.join(c).strip() for c in summary.columns.values]
        out = os.path.join(self.output_dir, 'runtime_summary.csv')
        summary.to_csv(out)
        print(f"[write] {out}")
        return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--results-dir', required=True)
    ap.add_argument('--results-file', default=None)
    args = ap.parse_args()
    az = BenchmarkAnalyzer(args.results_dir)
    az.load_results(args.results_file)
    az.summarize()


if __name__ == '__main__':
    main()
