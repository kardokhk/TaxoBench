#!/usr/bin/env python3

import os
import sys
import pandas as pd

def process_spingo_output(input_tsv, output_tsv):
    """
    Process SPINGO output format and convert to standardized format
    
    SPINGO output format (single-line with 6 columns):
    sequence_id    ignored_score    genus    genus_confidence    species    species_confidence
    
    Example:
    e09095b66760ea8feab73c7d5b843a19	1.00	Persephonella	0.98	Persephonella marina	0.95
    b15193fce14759d1c06728933e044af6	1.00	Bacteroides_H_857956	1.00	Bacteroides_H_857956 uniformis	0.95
    
    Where:
    - Column 1: sequence_id
    - Column 2: 1.00 (ignored)
    - Column 3: genus
    - Column 4: genus_confidence  
    - Column 5: species
    - Column 6: species_confidence
      Expected output format:
    Feature ID    Taxon    Confidence
    """
    print(f"Processing input file: {input_tsv}")
    
    # Read SPINGO output file
    with open(input_tsv, 'r') as f:
        lines = f.readlines()
    
    print(f"Read {len(lines)} lines from input file")
    
    output_rows = []
    
    print("Starting to process lines...")
    
    for i, line in enumerate(lines):
        line = line.strip()
        
        # Skip empty lines
        if not line:
            print(f"Skipping empty line {i+1}")
            continue
            
        # Split the line into parts
        parts = line.split('\t')
        
        # Check if this is the new single-line format (6 columns)
        if len(parts) == 6:
            sequence_id = parts[0]
            # parts[1] is ignored (1.00 in the example)
            genus = parts[2]
            genus_confidence = float(parts[3])
            species = parts[4] 
            species_confidence = float(parts[5])
            
            print(f"Processing sequence: {sequence_id}, genus: {genus} ({genus_confidence}), species: {species} ({species_confidence})")
            
        else:
            print(f"Skipping line {i+1} with unexpected format: {len(parts)} columns")
            continue
        
        # Build taxonomy string
        taxonomy_parts = []
        confidence_parts = []
        
        # Add empty levels for higher taxonomy (domain through family)
        prefixes = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
        for j, prefix in enumerate(prefixes):
            if j < 5:  # domain, phylum, class, order, family - empty
                taxonomy_parts.append(prefix)
                confidence_parts.append("0.0000")
            elif j == 5:  # genus
                if genus == "AMBIGUOUS":
                    taxonomy_parts.append(prefix)
                    confidence_parts.append("0.0000")
                else:
                    taxonomy_parts.append(f"{prefix}{genus}")
                    confidence_parts.append(f"{genus_confidence:.4f}")
            elif j == 6:  # species
                if species and species != "AMBIGUOUS" and genus != "AMBIGUOUS":
                    taxonomy_parts.append(f"{prefix}{species}")
                    confidence_parts.append(f"{species_confidence:.4f}")
                else:
                    taxonomy_parts.append(prefix)
                    confidence_parts.append("0.0000")
        
        taxonomy_str = "; ".join(taxonomy_parts)
        confidence_str = "; ".join(confidence_parts)
        
        output_rows.append((sequence_id, taxonomy_str, confidence_str))
    
    # Create output DataFrame
    out_df = pd.DataFrame(output_rows, columns=["Feature ID", "Taxon", "Confidence"])
    
    # Save to output file
    out_df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Processed {len(output_rows)} sequences from SPINGO output")

def main():
    if len(sys.argv) != 3:
        print("Usage: python process_spingo_output.py <input_tsv> <output_tsv>")
        sys.exit(1)
    
    input_tsv = sys.argv[1]
    output_tsv = sys.argv[2]
    
    if not os.path.exists(input_tsv):
        print(f"Error: Input file {input_tsv} does not exist")
        sys.exit(1)
    
    process_spingo_output(input_tsv, output_tsv)

if __name__ == "__main__":
    main()
