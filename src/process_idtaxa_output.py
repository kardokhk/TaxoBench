#!/usr/bin/env python3

import os
import sys
import pandas as pd
import re

def process_idtaxa_output(input_tsv, output_tsv):
    """
    Process idtaxa raw output to create properly formatted final taxonomy file.
    
    idtaxa specific processing:
    1. Remove the "Root" domain from taxonomy
    2. Remove duplicate prefixes that occur when idtaxa taxonomy already has prefixes
    3. Ensure consistent format with other classifiers
    """
    print(f"Processing idTaxa output from {input_tsv} to {output_tsv}")
      # Read input file
    try:
        df = pd.read_csv(input_tsv, sep='\t')
        print(f"Read {len(df)} sequences from input file")
    except Exception as e:
        print(f"Error reading input file: {e}")
        sys.exit(1)
    if df.empty:
        print(f"Error: {input_tsv} holds no classified sequence")
        sys.exit(1)
    
    output_rows = []
    
    # Handle potential column name variations (R converts spaces to dots)
    feature_id_col = 'Feature ID' if 'Feature ID' in df.columns else 'Feature.ID'
    
    for _, row in df.iterrows():
        feature_id = row[feature_id_col]
        taxonomy = row['Taxon']
        confidence = row['Confidence']
        
        # Split taxonomy and confidence
        tax_parts = [part.strip() for part in taxonomy.split(';')]
        conf_parts = [part.strip() for part in confidence.split(';')]
        
        # Process taxonomy parts
        processed_tax_parts = []
        processed_conf_parts = []
          # Standard prefixes for the 7 taxonomic levels
        prefixes = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']        # Check if input starts with Root (idtaxa specific)
        starts_with_root = len(tax_parts) > 0 and tax_parts[0].endswith('Root')
        
        # Process taxonomy levels
        if starts_with_root:
            # Skip Root and process remaining levels as domain, phylum, class, etc.
            start_index = 1
        else:
            start_index = 0
        
        # Process each taxonomic level
        for level_idx in range(7):  # Always process 7 levels
            if start_index + level_idx < len(tax_parts):
                tax_part = tax_parts[start_index + level_idx]
                conf_part = conf_parts[start_index + level_idx] if start_index + level_idx < len(conf_parts) else '0.0000'
                
                # Remove existing prefixes if present to avoid duplication
                # Handle cases like "p__d__Bacteria" by removing all prefixes
                clean_taxon = tax_part
                # Keep removing prefixes until no more are found
                while re.match(r'^[a-z]__', clean_taxon):
                    clean_taxon = re.sub(r'^[a-z]__', '', clean_taxon)
                
                # If the cleaned taxon is empty or just whitespace, treat as unclassified
                if not clean_taxon.strip():
                    processed_tax_parts.append(prefixes[level_idx])
                    processed_conf_parts.append('0.0000')
                else:
                    processed_tax_parts.append(f"{prefixes[level_idx]}{clean_taxon}")
                    processed_conf_parts.append(conf_part)
            else:
                # No more taxonomy data, fill with empty levels
                processed_tax_parts.append(prefixes[level_idx])
                processed_conf_parts.append('0.0000')
        
        # Pad to 7 levels if needed
        while len(processed_tax_parts) < 7:
            processed_tax_parts.append(prefixes[len(processed_tax_parts)])
            processed_conf_parts.append('0.0000')
        
        # Join back
        final_taxonomy = '; '.join(processed_tax_parts)
        final_confidence = '; '.join(processed_conf_parts)
        
        output_rows.append((feature_id, final_taxonomy, final_confidence))
    
    # Create output DataFrame
    out_df = pd.DataFrame(output_rows, columns=["Feature ID", "Taxon", "Confidence"])
    
    # Save to output file
    out_df.to_csv(output_tsv, sep='\t', index=False)
    print(f"Processed {len(output_rows)} sequences and saved to {output_tsv}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python process_idtaxa_output.py <input_tsv> <output_tsv>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
    
    process_idtaxa_output(input_file, output_file)
