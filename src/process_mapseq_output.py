#!/usr/bin/env python3
"""
Process MAPseq output format and convert to standardized format

MAPseq output format (.mseq):
Tab-delimited file with taxonomic classifications and confidence scores
Columns include: query, dbhit, bitscore, identity, matches, mismatches, gaps, 
query_start, query_end, dbhit_start, dbhit_end, strand, and then pairs of 
(taxon_name, combined_cf, score_cf) for each taxonomic level.

Example line (simplified):
query_1    dbhit_1    1492    0.995    ...    Bacteria    1    1    Pseudomonadota    0.927    1    ...

Expected output format:
Feature ID    Taxon    Confidence

Cutoff functionality:
When a cutoff threshold is specified (0.0-1.0), classifications below this threshold
will be considered unclassified. This allows for testing different confidence levels
and their impact on classification accuracy, similar to other classifiers like idTaxa
and taxonova.
"""

import os
import sys
import pandas as pd

def parse_mapseq_taxonomy_line(line, cutoff_threshold=0.0):
    """
    Parse a single line from MAPseq output and extract taxonomic information.
    
    Args:
        line: A line from the .mseq file
        cutoff_threshold: Minimum confidence threshold (0.0-1.0)
        
    Returns:
        tuple: (sequence_id, taxonomy_string, confidence_string)
    """
    parts = line.strip().split('\t')
    
    if len(parts) < 13:  # Minimum expected columns
        return None, None, None
    
    sequence_id = parts[0]
    
    # Taxonomic information starts after the alignment columns
    # Expected format: ..., Domain, combined_cf, score_cf, Phylum, combined_cf, score_cf, ...
    taxonomy_start_idx = 13  # After strand column
    
    # Extract taxonomic levels and confidence scores
    taxonomic_levels = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    rank_prefixes = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    taxonomy_parts = []
    confidence_parts = []
    
    # Process each taxonomic level
    for i, level in enumerate(taxonomic_levels):
        try:
            # Calculate the index for this taxonomic level
            # Each level has 3 columns: taxon_name, combined_cf, score_cf
            base_idx = taxonomy_start_idx + (i * 3)
            
            if base_idx + 2 < len(parts):
                taxon_name = parts[base_idx].strip()
                combined_cf = parts[base_idx + 1].strip()
                score_cf = parts[base_idx + 2].strip()
                
                # Use score_cf as the confidence (this is the actual confidence score)
                if taxon_name and taxon_name != '' and score_cf and score_cf != '':
                    try:
                        confidence_value = float(score_cf)
                        
                        # Apply cutoff threshold
                        if confidence_value >= cutoff_threshold:
                            # Add rank prefix to taxon name
                            prefixed_taxon = f"{rank_prefixes[i]}{taxon_name}"
                            taxonomy_parts.append(prefixed_taxon)
                            confidence_parts.append(f"{confidence_value:.4f}")
                        else:
                            # Below threshold - add empty level
                            if len(taxonomy_parts) > 0:  # Only add if we have some classification
                                empty_taxon = f"{rank_prefixes[i]}"
                                taxonomy_parts.append(empty_taxon)
                                confidence_parts.append("0.0000")
                            else:
                                break  # Stop if we hit an unconfident classification at the beginning
                    except ValueError:
                        # If score_cf is not a valid float, use combined_cf or default
                        try:
                            confidence_value = float(combined_cf) if combined_cf else 0.0
                            
                            # Apply cutoff threshold
                            if confidence_value >= cutoff_threshold:
                                prefixed_taxon = f"{rank_prefixes[i]}{taxon_name}"
                                taxonomy_parts.append(prefixed_taxon)
                                confidence_parts.append(f"{confidence_value:.4f}")
                            else:
                                # Below threshold - add empty level
                                if len(taxonomy_parts) > 0:  # Only add if we have some classification
                                    empty_taxon = f"{rank_prefixes[i]}"
                                    taxonomy_parts.append(empty_taxon)
                                    confidence_parts.append("0.0000")
                                else:
                                    break  # Stop if we hit an unconfident classification at the beginning
                        except ValueError:
                            # Could not parse confidence - treat as below threshold
                            if len(taxonomy_parts) > 0:  # Only add if we have some classification
                                empty_taxon = f"{rank_prefixes[i]}"
                                taxonomy_parts.append(empty_taxon)
                                confidence_parts.append("0.0000")
                            else:
                                break  # Stop if we hit an invalid classification at the beginning
                else:
                    # Add empty level with 0 confidence for incomplete classifications
                    if len(taxonomy_parts) > 0:  # Only add if we have some classification
                        empty_taxon = f"{rank_prefixes[i]}"
                        taxonomy_parts.append(empty_taxon)
                        confidence_parts.append("0.0000")
                    else:
                        break  # Stop if we hit an empty classification at the beginning
            else:
                break  # Not enough columns for this level
                
        except (IndexError, ValueError) as e:
            # Stop processing if we encounter an error
            break
    
    # Build final taxonomy and confidence strings
    if taxonomy_parts:
        # Create hierarchical taxonomy string
        taxonomy_str = "; ".join(taxonomy_parts)
        confidence_str = "; ".join(confidence_parts)
        return sequence_id, taxonomy_str, confidence_str
    else:
        # Unclassified sequence
        return sequence_id, "Unclassified", "0.0000"


def process_mapseq_output(input_mseq, output_tsv, cutoff_threshold=0.0):
    """
    Process MAPseq output format and convert to standardized format
    
    Args:
        input_mseq: Path to input .mseq file from MAPseq
        output_tsv: Path to output standardized TSV file
        cutoff_threshold: Minimum confidence threshold (0.0-1.0)
    """
    print(f"Processing MAPseq output file: {input_mseq}")
    if cutoff_threshold > 0.0:
        print(f"Using confidence cutoff threshold: {cutoff_threshold}")
    
    # Read MAPseq output file
    with open(input_mseq, 'r') as f:
        lines = f.readlines()
    
    print(f"Read {len(lines)} lines from input file")
    
    output_rows = []
    processed_count = 0
    skipped_count = 0
    
    print("Starting to process lines...")
    
    for i, line in enumerate(lines):
        line = line.strip()
        
        # Skip empty lines and comment lines
        if not line or line.startswith('#'):
            skipped_count += 1
            continue
        
        try:
            sequence_id, taxonomy_str, confidence_str = parse_mapseq_taxonomy_line(line, cutoff_threshold)
            
            if sequence_id and taxonomy_str and confidence_str:
                output_rows.append((sequence_id, taxonomy_str, confidence_str))
                processed_count += 1
            else:
                print(f"Warning: Could not parse line {i+1}: {line[:100]}...")
                skipped_count += 1
                
        except Exception as e:
            print(f"Error processing line {i+1}: {e}")
            print(f"Line content: {line[:100]}...")
            skipped_count += 1
    
    print(f"Processed {processed_count} sequences successfully")
    print(f"Skipped {skipped_count} lines")
    
    # Create output DataFrame
    if output_rows:
        out_df = pd.DataFrame(output_rows, columns=["Feature ID", "Taxon", "Confidence"])
        
        # Save to output file
        out_df.to_csv(output_tsv, sep='\t', index=False)
        print(f"Saved results to: {output_tsv}")
        print(f"Output contains {len(output_rows)} classified sequences")
    else:
        print("Warning: No sequences were successfully processed")
        # Create empty output file with headers
        out_df = pd.DataFrame(columns=["Feature ID", "Taxon", "Confidence"])
        out_df.to_csv(output_tsv, sep='\t', index=False)


def main():
    if len(sys.argv) < 3 or len(sys.argv) > 4:
        print("Usage: python process_mapseq_output.py <input_mseq> <output_tsv> [cutoff_threshold]")
        print("  cutoff_threshold: Optional minimum confidence threshold (0.0-1.0), default: 0.0")
        sys.exit(1)
    
    input_mseq = sys.argv[1]
    output_tsv = sys.argv[2]
    cutoff_threshold = float(sys.argv[3]) if len(sys.argv) == 4 else 0.0
    
    if not os.path.exists(input_mseq):
        print(f"Error: Input file {input_mseq} does not exist")
        sys.exit(1)
    
    if cutoff_threshold < 0.0 or cutoff_threshold > 1.0:
        print(f"Error: cutoff_threshold must be between 0.0 and 1.0, got {cutoff_threshold}")
        sys.exit(1)
    
    process_mapseq_output(input_mseq, output_tsv, cutoff_threshold)


if __name__ == "__main__":
    main()
