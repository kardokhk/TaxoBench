#!/usr/bin/env python3

import pandas as pd
import tempfile
import os
import sys

# Add the scripts directory to the path to import the processing function
scripts_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(scripts_dir)

from process_idtaxa_output import process_idtaxa_output

def test_idtaxa_processing():
    # Create test data that matches the problematic idtaxa output format shown by the user
    test_data = [
        ["3bac14d2805bf995afe9d316a5777b1e", "d__Root; p__d__Bacteria; c__p__Bacillota_I; o__c__Bacilli_A; f__o__Lactobacillales; g__f__Lactobacillaceae; s__g__Lactobacillus", "1.0000; 1.0000; 1.0000; 1.0000; 0.8473; 0.8473; 0.8277"],
        ["f7447c8e079023f4f2579d8580575dd3", "d__Root; p__d__Bacteria; c__p__Bacillota_I; o__c__Bacilli_A; f__o__Lactobacillales; g__f__Lactobacillaceae; s__g__Lactobacillus", "1.0000; 1.0000; 1.0000; 1.0000; 1.0000; 1.0000; 1.0000"],
        ["5ab13a4d1b95b491938622224ce5016e", "d__Root; p__d__Bacteria; c__p__Pseudomonadota; o__c__Gammaproteobacteria; f__o__Pseudomonadales_A_641035; g__f__Nitrincolaceae; s__g__Marinobacterium_D", "1.0000; 1.0000; 1.0000; 1.0000; 1.0000; 1.0000; 0.9654"]
    ]
    
    # Create a temporary input file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
        input_file = f.name
        f.write("Feature ID\tTaxon\tConfidence\n")
        for row in test_data:
            f.write(f"{row[0]}\t{row[1]}\t{row[2]}\n")
    
    # Create a temporary output file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
        output_file = f.name
    
    try:
        # Process the file
        process_idtaxa_output(input_file, output_file)
        
        # Read and display the results
        result_df = pd.read_csv(output_file, sep='\t')
        print("\nProcessed results:")
        print("=" * 80)
        for _, row in result_df.iterrows():
            print(f"Feature ID: {row['Feature ID']}")
            print(f"Taxon: {row['Taxon']}")
            print(f"Confidence: {row['Confidence']}")
            print("-" * 40)
              # Verify the results look correct
        first_row = result_df.iloc[0]
        expected_taxonomy = "d__Bacteria; p__Bacillota_I; c__Bacilli_A; o__Lactobacillales; f__Lactobacillaceae; g__Lactobacillus; s__"
        
        print(f"\nFirst row taxonomy: {first_row['Taxon']}")
        print(f"Expected taxonomy:  {expected_taxonomy}")
        
        if first_row['Taxon'] == expected_taxonomy:
            print("✅ Test PASSED: Taxonomy formatting is correct!")
        else:
            print("❌ Test FAILED: Taxonomy formatting is incorrect!")
            
    finally:
        # Clean up temporary files
        os.unlink(input_file)
        os.unlink(output_file)

if __name__ == "__main__":
    test_idtaxa_processing()
