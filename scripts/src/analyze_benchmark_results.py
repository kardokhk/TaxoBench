#!/usr/bin/env python3
"""
Benchmark Results Analysis and Visualization

This script analyzes the benchmark results from benchmark_classifiers.py and
generates comprehensive reports including plots and tables.
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import glob
import json

class BenchmarkAnalyzer:
    """Analyze and visualize benchmark results"""
    
    def __init__(self, results_dir):
        self.results_dir = results_dir
        self.output_dir = os.path.join(results_dir, 'analysis')
        os.makedirs(self.output_dir, exist_ok=True)
        
        # Set plotting style
        plt.style.use('default')
        sns.set_palette("husl")
        
    def load_results(self, results_file=None):
        """Load benchmark results from CSV file"""
        if results_file is None:
            # Find the most recent results file
            csv_files = glob.glob(os.path.join(self.results_dir, 'detailed_benchmark_results_*.csv'))
            if not csv_files:
                csv_files = glob.glob(os.path.join(self.results_dir, '*benchmark*.csv'))
            
            if not csv_files:
                raise FileNotFoundError("No benchmark results files found")
                
            results_file = max(csv_files, key=os.path.getctime)
            print(f"Loading results from: {results_file}")
        
        self.df = pd.read_csv(results_file)
        
        # Filter successful runs only
        self.df_success = self.df[self.df['success'] == True].copy()
        
        if len(self.df_success) == 0:
            raise ValueError("No successful benchmark runs found")
            
        print(f"Loaded {len(self.df)} total runs, {len(self.df_success)} successful")
        
        return self.df_success
    
    def generate_summary_statistics(self):
        """Generate summary statistics"""
        if not hasattr(self, 'df_success'):
            raise ValueError("Load results first using load_results()")
        
        # Overall summary - convert numpy types to native Python types for JSON serialization
        summary_stats = {
            'total_runs': int(len(self.df_success)),
            'classifiers': [str(x) for x in self.df_success['classifier'].unique()],
            'v_regions': [str(x) for x in self.df_success['v_region'].unique()],
            'sequence_counts': [int(x) for x in sorted(self.df_success['sequence_count'].unique())],
            'date_range': {
                'start': str(self.df_success['timestamp'].min()),
                'end': str(self.df_success['timestamp'].max())
            }
        }
        
        # Performance summary by classifier
        perf_summary = self.df_success.groupby('classifier').agg({
            'runtime_seconds': ['mean', 'std', 'min', 'max'],
            'max_memory_mb': ['mean', 'std', 'min', 'max'],
            'sequence_count': ['min', 'max']
        }).round(3)
        
        perf_summary.columns = ['_'.join(col).strip() for col in perf_summary.columns.values]
        
        # Save summary statistics
        with open(os.path.join(self.output_dir, 'summary_statistics.json'), 'w') as f:
            json.dump(summary_stats, f, indent=2)
        
        perf_summary.to_csv(os.path.join(self.output_dir, 'performance_summary_by_classifier.csv'))
        
        return summary_stats, perf_summary
    
    def plot_runtime_comparison(self):
        """Generate runtime comparison plots"""
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Classification Runtime Analysis', fontsize=16, fontweight='bold')
        
        # 1. Runtime by classifier (box plot)
        ax1 = axes[0, 0]
        sns.boxplot(data=self.df_success, x='classifier', y='runtime_seconds', ax=ax1)
        ax1.set_title('Runtime Distribution by Classifier')
        ax1.set_xlabel('Classifier')
        ax1.set_ylabel('Runtime (seconds)')
        ax1.tick_params(axis='x', rotation=45)
        
        # 2. Runtime vs sequence count (line plot)
        ax2 = axes[0, 1]
        runtime_by_size = self.df_success.groupby(['classifier', 'sequence_count'])['runtime_seconds'].mean().reset_index()
        for classifier in runtime_by_size['classifier'].unique():
            data = runtime_by_size[runtime_by_size['classifier'] == classifier]
            ax2.plot(data['sequence_count'], data['runtime_seconds'], marker='o', label=classifier)
        
        ax2.set_title('Runtime vs Sequence Count')
        ax2.set_xlabel('Number of Sequences')
        ax2.set_ylabel('Average Runtime (seconds)')
        ax2.legend()
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        
        # 3. Runtime by V-region (bar plot)
        ax3 = axes[1, 0]
        runtime_by_region = self.df_success.groupby(['v_region', 'classifier'])['runtime_seconds'].mean().reset_index()
        runtime_pivot = runtime_by_region.pivot(index='v_region', columns='classifier', values='runtime_seconds')
        runtime_pivot.plot(kind='bar', ax=ax3)
        ax3.set_title('Average Runtime by V-Region')
        ax3.set_xlabel('V-Region')
        ax3.set_ylabel('Average Runtime (seconds)')
        ax3.legend(title='Classifier', bbox_to_anchor=(1.05, 1), loc='upper left')
        ax3.tick_params(axis='x', rotation=45)
        
        # 4. Runtime scaling efficiency (normalized by smallest dataset)
        ax4 = axes[1, 1]
        scaling_data = []
        for classifier in self.df_success['classifier'].unique():
            classifier_data = self.df_success[self.df_success['classifier'] == classifier]
            size_means = classifier_data.groupby('sequence_count')['runtime_seconds'].mean()
            if len(size_means) > 1:
                baseline = size_means.iloc[0]  # Smallest size
                normalized = size_means / baseline
                for size, norm_time in normalized.items():
                    scaling_data.append({
                        'classifier': classifier,
                        'sequence_count': size,
                        'normalized_runtime': norm_time
                    })
        
        if scaling_data:
            scaling_df = pd.DataFrame(scaling_data)
            for classifier in scaling_df['classifier'].unique():
                data = scaling_df[scaling_df['classifier'] == classifier]
                ax4.plot(data['sequence_count'], data['normalized_runtime'], marker='o', label=classifier)
        
        ax4.set_title('Runtime Scaling (Normalized to Smallest Dataset)')
        ax4.set_xlabel('Number of Sequences')
        ax4.set_ylabel('Normalized Runtime')
        ax4.legend()
        ax4.set_xscale('log')
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'runtime_analysis.png'), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(self.output_dir, 'runtime_analysis.pdf'), bbox_inches='tight')
        plt.close()
        
    def plot_memory_comparison(self):
        """Generate memory usage comparison plots"""
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Memory Usage Analysis', fontsize=16, fontweight='bold')
        
        # 1. Memory by classifier (box plot)
        ax1 = axes[0, 0]
        sns.boxplot(data=self.df_success, x='classifier', y='max_memory_mb', ax=ax1)
        ax1.set_title('Memory Usage Distribution by Classifier')
        ax1.set_xlabel('Classifier')
        ax1.set_ylabel('Peak Memory Usage (MB)')
        ax1.tick_params(axis='x', rotation=45)
        
        # 2. Memory vs sequence count (line plot)
        ax2 = axes[0, 1]
        memory_by_size = self.df_success.groupby(['classifier', 'sequence_count'])['max_memory_mb'].mean().reset_index()
        for classifier in memory_by_size['classifier'].unique():
            data = memory_by_size[memory_by_size['classifier'] == classifier]
            ax2.plot(data['sequence_count'], data['max_memory_mb'], marker='o', label=classifier)
        
        ax2.set_title('Memory Usage vs Sequence Count')
        ax2.set_xlabel('Number of Sequences')
        ax2.set_ylabel('Average Peak Memory (MB)')
        ax2.legend()
        ax2.set_xscale('log')
        
        # 3. Memory by V-region (bar plot)
        ax3 = axes[1, 0]
        memory_by_region = self.df_success.groupby(['v_region', 'classifier'])['max_memory_mb'].mean().reset_index()
        memory_pivot = memory_by_region.pivot(index='v_region', columns='classifier', values='max_memory_mb')
        memory_pivot.plot(kind='bar', ax=ax3)
        ax3.set_title('Average Memory Usage by V-Region')
        ax3.set_xlabel('V-Region')
        ax3.set_ylabel('Average Peak Memory (MB)')
        ax3.legend(title='Classifier', bbox_to_anchor=(1.05, 1), loc='upper left')
        ax3.tick_params(axis='x', rotation=45)
        
        # 4. Memory efficiency (MB per sequence)
        ax4 = axes[1, 1]
        efficiency_data = self.df_success.copy()
        efficiency_data['memory_per_sequence'] = efficiency_data['max_memory_mb'] / efficiency_data['sequence_count']
        
        sns.boxplot(data=efficiency_data, x='classifier', y='memory_per_sequence', ax=ax4)
        ax4.set_title('Memory Efficiency (MB per Sequence)')
        ax4.set_xlabel('Classifier')
        ax4.set_ylabel('Memory per Sequence (MB)')
        ax4.tick_params(axis='x', rotation=45)
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'memory_analysis.png'), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(self.output_dir, 'memory_analysis.pdf'), bbox_inches='tight')
        plt.close()
        
    def plot_combined_performance(self):
        """Generate combined performance plots"""
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Combined Performance Analysis', fontsize=16, fontweight='bold')
        
        # 1. Runtime vs Memory scatter plot
        ax1 = axes[0, 0]
        for classifier in self.df_success['classifier'].unique():
            data = self.df_success[self.df_success['classifier'] == classifier]
            ax1.scatter(data['runtime_seconds'], data['max_memory_mb'], 
                       label=classifier, alpha=0.7, s=50)
        
        ax1.set_title('Runtime vs Memory Usage')
        ax1.set_xlabel('Runtime (seconds)')
        ax1.set_ylabel('Peak Memory (MB)')
        ax1.legend()
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        
        # 2. Performance heatmap
        ax2 = axes[0, 1]
        # Calculate average runtime and memory for each classifier-region combination
        heatmap_data = self.df_success.groupby(['classifier', 'v_region']).agg({
            'runtime_seconds': 'mean',
            'max_memory_mb': 'mean'
        }).round(2)
        
        # Create efficiency score (inverse of runtime + memory, normalized)
        runtime_norm = (heatmap_data['runtime_seconds'] - heatmap_data['runtime_seconds'].min()) / \
                       (heatmap_data['runtime_seconds'].max() - heatmap_data['runtime_seconds'].min())
        memory_norm = (heatmap_data['max_memory_mb'] - heatmap_data['max_memory_mb'].min()) / \
                      (heatmap_data['max_memory_mb'].max() - heatmap_data['max_memory_mb'].min())
        
        heatmap_data['efficiency_score'] = 1 - (runtime_norm + memory_norm) / 2
        efficiency_pivot = heatmap_data['efficiency_score'].reset_index().pivot(
            index='classifier', columns='v_region', values='efficiency_score'
        )
        
        sns.heatmap(efficiency_pivot, annot=True, cmap='RdYlGn', center=0.5, ax=ax2)
        ax2.set_title('Performance Efficiency Score\n(Higher = Better)')
        
        # 3. Scalability analysis
        ax3 = axes[1, 0]
        scalability_data = []
        for classifier in self.df_success['classifier'].unique():
            classifier_data = self.df_success[self.df_success['classifier'] == classifier]
            if len(classifier_data['sequence_count'].unique()) > 1:
                # Calculate slope of runtime vs sequence count (log-log)
                sizes = classifier_data['sequence_count'].values
                runtimes = classifier_data['runtime_seconds'].values
                
                if len(sizes) > 1:
                    log_sizes = np.log10(sizes)
                    log_runtimes = np.log10(runtimes)
                    slope = np.polyfit(log_sizes, log_runtimes, 1)[0]
                    
                    scalability_data.append({
                        'classifier': classifier,
                        'scalability_slope': slope
                    })
        
        if scalability_data:
            scalability_df = pd.DataFrame(scalability_data)
            bars = ax3.bar(scalability_df['classifier'], scalability_df['scalability_slope'])
            ax3.set_title('Runtime Scalability\n(Lower slope = Better scaling)')
            ax3.set_xlabel('Classifier')
            ax3.set_ylabel('Log-Log Slope (Runtime vs Size)')
            ax3.tick_params(axis='x', rotation=45)
            
            # Color bars based on slope (green = better scaling)
            slopes = scalability_df['scalability_slope']
            colors = plt.cm.RdYlGn_r((slopes - slopes.min()) / (slopes.max() - slopes.min()))
            for bar, color in zip(bars, colors):
                bar.set_color(color)
        
        # 4. Resource efficiency by dataset size
        ax4 = axes[1, 1]
        efficiency_by_size = self.df_success.groupby(['classifier', 'sequence_count']).agg({
            'runtime_seconds': 'mean',
            'max_memory_mb': 'mean'
        }).reset_index()
        
        # Calculate throughput (sequences per second)
        efficiency_by_size['throughput'] = efficiency_by_size['sequence_count'] / efficiency_by_size['runtime_seconds']
        
        for classifier in efficiency_by_size['classifier'].unique():
            data = efficiency_by_size[efficiency_by_size['classifier'] == classifier]
            ax4.plot(data['sequence_count'], data['throughput'], marker='o', label=classifier)
        
        ax4.set_title('Throughput (Sequences per Second)')
        ax4.set_xlabel('Dataset Size (Number of Sequences)')
        ax4.set_ylabel('Throughput (sequences/second)')
        ax4.legend()
        ax4.set_xscale('log')
        ax4.set_yscale('log')
        
        plt.tight_layout()
        plt.savefig(os.path.join(self.output_dir, 'combined_performance.png'), dpi=300, bbox_inches='tight')
        plt.savefig(os.path.join(self.output_dir, 'combined_performance.pdf'), bbox_inches='tight')
        plt.close()
        
    def generate_performance_tables(self):
        """Generate detailed performance tables"""
        
        # 1. Summary table by classifier
        summary_table = self.df_success.groupby('classifier').agg({
            'runtime_seconds': ['mean', 'std', 'min', 'max'],
            'max_memory_mb': ['mean', 'std', 'min', 'max'],
            'sequence_count': ['count']
        }).round(3)
        
        summary_table.columns = ['_'.join(col).strip() for col in summary_table.columns.values]
        summary_table.columns = [col.replace('sequence_count_count', 'total_runs') for col in summary_table.columns]
        
        # 2. Performance by dataset size
        size_table = self.df_success.groupby(['classifier', 'sequence_count']).agg({
            'runtime_seconds': ['mean', 'std'],
            'max_memory_mb': ['mean', 'std']
        }).round(3)
        
        size_table.columns = ['_'.join(col).strip() for col in size_table.columns.values]
        
        # 3. Performance by V-region
        region_table = self.df_success.groupby(['classifier', 'v_region']).agg({
            'runtime_seconds': ['mean', 'std'],
            'max_memory_mb': ['mean', 'std']
        }).round(3)
        
        region_table.columns = ['_'.join(col).strip() for col in region_table.columns.values]
        
        # 4. Ranking table
        ranking_data = []
        for metric in ['runtime_seconds', 'max_memory_mb']:
            metric_avg = self.df_success.groupby('classifier')[metric].mean().sort_values()
            for rank, (classifier, value) in enumerate(metric_avg.items(), 1):
                ranking_data.append({
                    'classifier': classifier,
                    'metric': metric,
                    'rank': rank,
                    'value': round(value, 3)
                })
        
        ranking_table = pd.DataFrame(ranking_data)
        ranking_pivot = ranking_table.pivot(index='classifier', columns='metric', values='rank')
        ranking_pivot['average_rank'] = ranking_pivot.mean(axis=1)
        ranking_pivot = ranking_pivot.sort_values('average_rank')
        
        # Save all tables
        with pd.ExcelWriter(os.path.join(self.output_dir, 'performance_tables.xlsx')) as writer:
            summary_table.to_excel(writer, sheet_name='Summary_by_Classifier')
            size_table.to_excel(writer, sheet_name='Performance_by_Size')
            region_table.to_excel(writer, sheet_name='Performance_by_Region')
            ranking_pivot.to_excel(writer, sheet_name='Classifier_Rankings')
        
        # Also save as CSV files
        summary_table.to_csv(os.path.join(self.output_dir, 'summary_by_classifier.csv'))
        size_table.to_csv(os.path.join(self.output_dir, 'performance_by_size.csv'))
        region_table.to_csv(os.path.join(self.output_dir, 'performance_by_region.csv'))
        ranking_pivot.to_csv(os.path.join(self.output_dir, 'classifier_rankings.csv'))
        
        return summary_table, size_table, region_table, ranking_pivot
    
    def generate_report(self):
        """Generate a comprehensive HTML report"""
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Taxonomic Classifier Performance Benchmark Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 40px; }}
                h1, h2 {{ color: #333; }}
                .summary {{ background-color: #f5f5f5; padding: 20px; margin: 20px 0; }}
                .plot {{ text-align: center; margin: 20px 0; }}
                table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
                th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
                th {{ background-color: #f2f2f2; }}
                .metric {{ font-weight: bold; color: #007acc; }}
            </style>
        </head>
        <body>
            <h1>Taxonomic Classifier Performance Benchmark Report</h1>
            <p><strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            
            <div class="summary">
                <h2>Executive Summary</h2>
                <p>This report presents the performance benchmarking results for {len(self.df_success['classifier'].unique())} 
                taxonomic classifiers across {len(self.df_success['v_region'].unique())} amplicon regions 
                and {len(self.df_success['sequence_count'].unique())} different dataset sizes.</p>
                
                <ul>
                    <li><strong>Total successful runs:</strong> {len(self.df_success)}</li>
                    <li><strong>Classifiers tested:</strong> {', '.join(sorted(self.df_success['classifier'].unique()))}</li>
                    <li><strong>V-regions tested:</strong> {', '.join(sorted(self.df_success['v_region'].unique()))}</li>
                    <li><strong>Dataset sizes:</strong> {', '.join(map(str, sorted(self.df_success['sequence_count'].unique())))} sequences</li>
                </ul>
            </div>
            
            <h2>Performance Analysis</h2>
            
            <div class="plot">
                <h3>Runtime Analysis</h3>
                <img src="runtime_analysis.png" alt="Runtime Analysis" style="max-width: 100%;">
            </div>
            
            <div class="plot">
                <h3>Memory Usage Analysis</h3>
                <img src="memory_analysis.png" alt="Memory Analysis" style="max-width: 100%;">
            </div>
            
            <div class="plot">
                <h3>Combined Performance Analysis</h3>
                <img src="combined_performance.png" alt="Combined Performance" style="max-width: 100%;">
            </div>
            
            <h2>Key Findings</h2>
            <div class="summary">
        """
        
        # Add key findings
        fastest_classifier = self.df_success.groupby('classifier')['runtime_seconds'].mean().idxmin()
        most_efficient_memory = self.df_success.groupby('classifier')['max_memory_mb'].mean().idxmin()
        
        html_content += f"""
                <ul>
                    <li><span class="metric">Fastest Classifier:</span> {fastest_classifier}</li>
                    <li><span class="metric">Most Memory Efficient:</span> {most_efficient_memory}</li>
                    <li><span class="metric">Runtime Range:</span> 
                        {self.df_success['runtime_seconds'].min():.2f} - {self.df_success['runtime_seconds'].max():.2f} seconds</li>
                    <li><span class="metric">Memory Range:</span> 
                        {self.df_success['max_memory_mb'].min():.1f} - {self.df_success['max_memory_mb'].max():.1f} MB</li>
                </ul>
            </div>
            
            <h2>Detailed Tables</h2>
            <p>For detailed performance tables, please refer to the Excel file: 
            <a href="performance_tables.xlsx">performance_tables.xlsx</a></p>
            
            <h2>Data Files</h2>
            <ul>
                <li><a href="summary_by_classifier.csv">Summary by Classifier (CSV)</a></li>
                <li><a href="performance_by_size.csv">Performance by Dataset Size (CSV)</a></li>
                <li><a href="performance_by_region.csv">Performance by V-Region (CSV)</a></li>
                <li><a href="classifier_rankings.csv">Classifier Rankings (CSV)</a></li>
            </ul>
            
        </body>
        </html>
        """
        
        # Save HTML report
        with open(os.path.join(self.output_dir, 'benchmark_report.html'), 'w') as f:
            f.write(html_content)
    
    def run_complete_analysis(self, results_file=None):
        """Run complete analysis pipeline"""
        print("Loading benchmark results...")
        self.load_results(results_file)
        
        print("Generating summary statistics...")
        self.generate_summary_statistics()
        
        print("Creating runtime analysis plots...")
        self.plot_runtime_comparison()
        
        print("Creating memory analysis plots...")
        self.plot_memory_comparison()
        
        print("Creating combined performance plots...")
        self.plot_combined_performance()
        
        print("Generating performance tables...")
        self.generate_performance_tables()
        
        print("Generating HTML report...")
        self.generate_report()
        
        print(f"Analysis complete! Results saved to: {self.output_dir}")

def main():
    """Main function"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze benchmark results')
    parser.add_argument('--results-dir', required=True, 
                       help='Directory containing benchmark results')
    parser.add_argument('--results-file', 
                       help='Specific results CSV file to analyze')
    
    args = parser.parse_args()
    
    analyzer = BenchmarkAnalyzer(args.results_dir)
    analyzer.run_complete_analysis(args.results_file)

if __name__ == "__main__":
    main()
