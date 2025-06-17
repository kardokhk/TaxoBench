#!/usr/bin/env python3

import tempfile
import os
from process_spingo_output import process_spingo_output

def test_spingo_processing():
    """Test SPINGO output processing with sample data"""
    
    print("Starting SPINGO processing test...")
      # Sample SPINGO output (new single-line format)
    sample_input = """e09095b66760ea8feab73c7d5b843a19	1.00	Persephonella	0.98	Persephonella marina	0.95
b15193fce14759d1c06728933e044af6	1.00	Bacteroides_H_857956	1.00	Bacteroides_H_857956 uniformis	0.95
31cd3a00123ecbf1b4fbc3e1def0ad	1.00	Chlorobium	1.00	Chlorobium phaeoclathratiforme	0.47
e3055f1b3a2ef5ffe239567f02e0e758	1.00	Streptococcus	1.00	Streptococcus agalactiae	0.77
c06ccd1f7d57566ef669942328b1a946	1.00	AMBIGUOUS	0.00	AMBIGUOUS	0.00"""
    
    print("Creating temporary files...")
    
    # Create temporary input file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as temp_input:
        temp_input.write(sample_input)
        temp_input_path = temp_input.name
    
    print(f"Input file created: {temp_input_path}")
    
    # Create temporary output file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as temp_output:
        temp_output_path = temp_output.name
    
    print(f"Output file created: {temp_output_path}")
    
    try:
        print("About to call process_spingo_output...")
        # Process the file
        process_spingo_output(temp_input_path, temp_output_path)
        print("process_spingo_output completed!")
        
        # Read and display results
        print("Reading results...")
        with open(temp_output_path, 'r') as f:
            result = f.read()
            print("SPINGO Processing Test Results:")
            print("=" * 50)
            print(result)
            print("=" * 50)
            print("Test completed successfully!")
            
    except Exception as e:
        print(f"Error occurred: {e}")
        import traceback
        traceback.print_exc()
            
    finally:
        print("Cleaning up...")
        # Clean up
        if os.path.exists(temp_input_path):
            os.unlink(temp_input_path)
        if os.path.exists(temp_output_path):
            os.unlink(temp_output_path)
        print("Cleanup completed.")

if __name__ == "__main__":
    test_spingo_processing()
