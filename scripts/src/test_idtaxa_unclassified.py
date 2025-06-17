#!/usr/bin/env python3

import pandas as pd
import tempfile
import os
import sys

# Add the scripts directory to the path to import the processing function
scripts_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scripts', 'src')
sys.path.append(scripts_dir)

from process_idtaxa_output import process_idtaxa_output
from evaluate_results import evaluate_classification

def test_idtaxa_unclassified_behavior():
    """Test how IDTAXA processing handles the 'unclassified_' format"""
    
    # Test data with the exact format from your examples
    test_data = [
        # Case 1: unclassified at species level (should be treated as correctly unclassified)
        ["3bac14d2805bf995afe9d316a5777b1e", "Root; d__Bacteria; p__Bacillota_I; c__Bacilli_A; o__Lactobacillales; f__Lactobacillaceae; g__Lactobacillus; unclassified_g__Lactobacillus", "1.0000; 1.0000; 1.0000; 1.0000; 0.8942; 0.8942; 0.8942; 0.8942"],
        
        # Case 2: Normal species classification
        ["f7447c8e079023f4f2579d8580575dd3", "Root; d__Bacteria; p__Bacillota_I; c__Bacilli_A; o__Lactobacillales; f__Lactobacillaceae; g__Lactobacillus; s__Lactobacillus iners", "1.0000; 1.0000; 1.0000; 1.0000; 1.0000; 1.0000; 0.9918; 0.9688"],
        
        # Case 3: unclassified at genus level
        ["5ab13a4d1b95b491938622224ce5016e", "Root; d__Bacteria; p__Pseudomonadota; c__Gammaproteobacteria; o__Pseudomonadales_A_641035; f__Nitrincolaceae; g__Marinobacterium_D; unclassified_g__Marinobacterium_D", "1.0000; 1.0000; 1.0000; 1.0000; 1.0000; 0.9925; 0.9722; 0.9722"],
        
        # Case 4: unclassified at genus level from another example
        ["e51f253cbc31450a3cab5eb0d5d787cb", "Root; d__Bacteria; p__Bacillota_I; c__Bacilli_A; o__Lactobacillales; f__Carnobacteriaceae; g__Carnobacterium_A_320743; unclassified_g__Carnobacterium_A_320743", "1.0000; 1.0000; 1.0000; 1.0000; 1.0000; 0.9388; 0.9241; 0.9241"],
    ]
    
    # Create a temporary input file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
        input_file = f.name
        f.write("Feature.ID\tTaxon\tConfidence\n")
        for row in test_data:
            f.write(f"{row[0]}\t{row[1]}\t{row[2]}\n")
    
    # Create a temporary output file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
        output_file = f.name
    
    try:
        print("=== Testing IDTAXA 'unclassified_' format processing ===\n")
        
        # Process the file
        process_idtaxa_output(input_file, output_file)
        
        # Read and display the results
        result_df = pd.read_csv(output_file, sep='\t')
        
        print("Raw IDTAXA output vs Processed output:")
        print("=" * 100)
        
        for i, (_, row) in enumerate(result_df.iterrows()):
            print(f"\nExample {i+1}:")
            print(f"Feature ID: {row['Feature ID']}")
            print(f"Raw input:  {test_data[i][1]}")
            print(f"Processed:  {row['Taxon']}")
            print(f"Confidence: {row['Confidence']}")
            
            # Check if this contains unclassified format
            if 'unclassified_' in row['Taxon']:
                print(f"⚠️  Contains 'unclassified_' format: This may be incorrectly evaluated!")
            
        print("\n" + "=" * 100)
        print("\n=== Testing evaluation behavior with different expected taxonomies ===\n")
        
        # Test evaluation with different expected taxonomies
        test_cases = [
            {
                'name': 'Case 1: Observed unclassified vs Expected unclassified',
                'observed': 'd__Bacteria; p__Bacillota_I; c__Bacilli_A; o__Lactobacillales; f__Lactobacillaceae; g__Lactobacillus; s__unclassified_g__Lactobacillus',
                'expected': 'd__Bacteria; p__Bacillota_I; c__Bacilli_A; o__Lactobacillales; f__Lactobacillaceae; g__Lactobacillus; s__'
            },
            {
                'name': 'Case 2: Observed unclassified vs Expected specific species',
                'observed': 'd__Bacteria; p__Bacillota_I; c__Bacilli_A; o__Lactobacillales; f__Lactobacillaceae; g__Lactobacillus; s__unclassified_g__Lactobacillus',
                'expected': 'd__Bacteria; p__Bacillota_I; c__Bacilli_A; o__Lactobacillales; f__Lactobacillaceae; g__Lactobacillus; s__Lactobacillus iners'
            },
            {
                'name': 'Case 3: Observed specific vs Expected unclassified',
                'observed': 'd__Bacteria; p__Bacillota_I; c__Bacilli_A; o__Lactobacillales; f__Lactobacillaceae; g__Lactobacillus; s__Lactobacillus iners',
                'expected': 'd__Bacteria; p__Bacillota_I; c__Bacilli_A; o__Lactobacillales; f__Lactobacillaceae; g__Lactobacillus; s__'
            }
        ]
        
        for case in test_cases:
            result = evaluate_classification(case['observed'], case['expected'])
            print(f"{case['name']}:")
            print(f"  Observed: {case['observed']}")
            print(f"  Expected: {case['expected']}")
            print(f"  Evaluation result: {result}")
            print()
            
    finally:
        # Clean up temporary files
        if os.path.exists(input_file):
            os.unlink(input_file)
        if os.path.exists(output_file):
            os.unlink(output_file)

if __name__ == "__main__":
    test_idtaxa_unclassified_behavior()
