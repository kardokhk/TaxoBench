#!/usr/bin/env python3
"""
Process Kraken2 output format and convert to standardized format

Kraken2 output format (--output):
C/U    sequence_id    taxonomy_id    sequence_length    lca_mapping

Example:
C	23334	Bacteria (taxid 1)	253	1:223
U	66653	0	252	0:252

Where:
- Column 1: C (classified) or U (unclassified)
- Column 2: sequence_id
- Column 3: taxonomy assignment (with taxid) or 0 for unclassified
- Column 4: sequence_length
- Column 5: LCA mapping information

Kraken2 report format (--report):
percentage  reads_clade reads_direct rank taxid taxonomy_name

Expected output format:
Feature ID    Taxon    Confidence
"""

import os
import sys
import pandas as pd
import re

def parse_kraken2_taxonomy(taxonomy_str, taxid_to_info=None):
    """
    Parse Kraken2 taxonomy string and extract taxonomic information.
    
    Args:
        taxonomy_str: String like "Bacteria (taxid 1)" when --use-names=True,
                     or just "1" when --use-names=False, or "0" for unclassified
        taxid_to_info: Optional dict mapping taxid to taxonomic info from report file
        
    Returns:
        tuple: (taxonomy_name, taxid)
    """
    if taxonomy_str == "0":
        return "Unclassified", 0
    
    # Case 1: --use-names=True format: "Bacteria (taxid 1)"
    match = re.match(r'^(.+?)\s+\(taxid\s+(\d+)\)$', taxonomy_str)
    if match:
        taxonomy_name = match.group(1).strip()
        taxid = int(match.group(2))
        return taxonomy_name, taxid
    
    # Case 2: --use-names=False format: just the taxid as a string (e.g., "1", "594")
    if taxonomy_str.isdigit():
        taxid = int(taxonomy_str)
        # Look up the taxonomy name from the report file if available
        if taxid_to_info and taxid in taxid_to_info:
            taxonomy_name = taxid_to_info[taxid]['name']
            return taxonomy_name, taxid
        else:
            # If no report info available, return taxid as name (fallback)
            return f"taxid_{taxid}", taxid
    
    # Fallback if format is different
    return taxonomy_str, 0

def parse_kraken2_report(report_file):
    """
    Parse Kraken2 report file to build taxonomic hierarchy information.
    
    Args:
        report_file: Path to Kraken2 report file
        
    Returns:
        dict: Mapping of taxid to taxonomic level information
        dict: Mapping of taxonomy name to taxonomic level information
        dict: Mapping of taxid to complete lineage  
    """
    taxid_to_info = {}
    name_to_info = {}
    taxid_to_lineage = {}
    
    if not os.path.exists(report_file):
        print(f"Warning: Report file {report_file} not found. Using basic processing.")
        return taxid_to_info, name_to_info, taxid_to_lineage
    
    rank_mapping = {
        'D': 0,  # Domain
        'P': 1,  # Phylum
        'C': 2,  # Class
        'O': 3,  # Order
        'F': 4,  # Family
        'G': 5,  # Genus
        'S': 6   # Species
    }
    
    # Stack to track hierarchy based on indentation
    hierarchy_stack = []
    
    with open(report_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                percentage = float(parts[0])
                rank = parts[3]
                taxid = int(parts[4])
                taxonomy_name = parts[5]
                
                # Extract indentation level and clean name
                indent_level = len(taxonomy_name) - len(taxonomy_name.lstrip())
                clean_name = taxonomy_name.strip()
                
                info = {
                    'percentage': percentage,
                    'rank': rank,
                    'taxid': taxid,
                    'name': clean_name,
                    'indent_level': indent_level,
                    'taxonomic_level': rank_mapping.get(rank, -1)
                }
                
                # Maintain hierarchy stack based on indentation
                # Remove items from stack that are at same or deeper level
                while hierarchy_stack and hierarchy_stack[-1]['indent_level'] >= indent_level:
                    hierarchy_stack.pop()
                
                # Add current item to stack
                hierarchy_stack.append(info)
                
                # Build complete lineage for this taxid
                lineage = [''] * 7  # Initialize with empty strings for d__, p__, c__, o__, f__, g__, s__
                
                for ancestor in hierarchy_stack:
                    level = ancestor['taxonomic_level']
                    if 0 <= level <= 6:  # Valid taxonomic level
                        lineage[level] = ancestor['name']
                
                taxid_to_lineage[taxid] = lineage
                taxid_to_info[taxid] = info
                name_to_info[clean_name] = info
    
    return taxid_to_info, name_to_info, taxid_to_lineage

def build_hierarchical_taxonomy(taxonomy_name, taxid, name_to_info, taxid_to_info, taxid_to_lineage):
    """
    Build hierarchical taxonomy string using report information.
    
    Args:
        taxonomy_name: The classified taxonomy name
        taxid: The taxonomy ID
        name_to_info: Mapping from names to taxonomic info
        taxid_to_info: Mapping from taxids to taxonomic info
        taxid_to_lineage: Mapping from taxids to complete lineage
        
    Returns:
        tuple: (taxonomy_string, confidence_string, taxonomic_level)
    """
    prefixes = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    taxonomy_parts = [''] * 7  # Initialize with empty strings
    confidence_parts = ['0.0000'] * 7
    
    if taxonomy_name == "Unclassified":
        # All unclassified
        taxonomy_str = "; ".join([prefix for prefix in prefixes])
        confidence_str = "; ".join(confidence_parts)
        return taxonomy_str, confidence_str, -1
    
    # First, try to get the complete lineage from the report file
    if taxid in taxid_to_lineage:
        lineage = taxid_to_lineage[taxid]
        info = taxid_to_info.get(taxid, {})
        
        # Fill in the lineage with appropriate confidence scores
        for i, (prefix, taxon_name) in enumerate(zip(prefixes, lineage)):
            if taxon_name:  # If we have a taxon name at this level
                taxonomy_parts[i] = taxon_name
                # Use the percentage from the classified taxon for the specific level
                # and a default confidence for higher levels in the lineage
                if i == info.get('taxonomic_level', -1):
                    # This is the level where the classification was made
                    confidence_parts[i] = f"{info.get('percentage', 0.0)/100.0:.4f}"
                else:
                    # Higher level in the lineage - use default confidence
                    confidence_parts[i] = "1.0000"  # Assume high confidence for lineage
        
        actual_level = info.get('taxonomic_level', -1)
    else:
        # Fallback: try to find the taxonomic information for this classification
        info = None
        if taxonomy_name in name_to_info:
            info = name_to_info[taxonomy_name]
        elif taxid in taxid_to_info:
            info = taxid_to_info[taxid]
        
        if info and info['taxonomic_level'] >= 0:
            # Place the classification at its correct taxonomic level
            level = info['taxonomic_level']
            taxonomy_parts[level] = taxonomy_name
            confidence_parts[level] = f"{info['percentage']/100.0:.4f}"
            actual_level = level
        else:
            # Fallback: try to guess the taxonomic level based on common patterns
            level = guess_taxonomic_level(taxonomy_name)
            if level >= 0:
                taxonomy_parts[level] = taxonomy_name
                confidence_parts[level] = "1.0000"  # Default confidence
                actual_level = level
            else:
                # Place at genus level as a fallback (like SPINGO does)
                taxonomy_parts[5] = taxonomy_name  # genus
                confidence_parts[5] = "1.0000"
                actual_level = 5
    
    # Build the final taxonomy string
    final_taxonomy_parts = []
    final_confidence_parts = []
    
    for i, (prefix, tax_part, conf_part) in enumerate(zip(prefixes, taxonomy_parts, confidence_parts)):
        if tax_part:
            final_taxonomy_parts.append(f"{prefix}{tax_part}")
        else:
            final_taxonomy_parts.append(prefix)
        final_confidence_parts.append(conf_part)
    
    taxonomy_str = "; ".join(final_taxonomy_parts)
    confidence_str = "; ".join(final_confidence_parts)
    
    return taxonomy_str, confidence_str, actual_level

def guess_taxonomic_level(taxonomy_name):
    """
    Guess taxonomic level based on common naming patterns.
    
    Args:
        taxonomy_name: The taxonomy name to analyze
        
    Returns:
        int: Taxonomic level (0-6) or -1 if unknown
    """
    name_lower = taxonomy_name.lower()
    
    # Common domain names
    if name_lower in ['bacteria', 'archaea', 'eukaryota']:
        return 0  # Domain
    
    # Common phylum endings
    if (name_lower.endswith('ota') or name_lower.endswith('phyta') or 
        name_lower.endswith('mycota') or 'proteobacter' in name_lower or
        'bacillota' in name_lower or 'pseudomonadota' in name_lower):
        return 1  # Phylum
    
    # Common class endings  
    if (name_lower.endswith('ia') or name_lower.endswith('i') or
        'bacilli' in name_lower or 'clostridia' in name_lower or
        'proteobacteria' in name_lower):
        return 2  # Class
    
    # Common order endings
    if (name_lower.endswith('ales') or name_lower.endswith('iales')):
        return 3  # Order
    
    # Common family endings
    if (name_lower.endswith('aceae') or name_lower.endswith('idae')):
        return 4  # Family
    
    # Common genus patterns (single word, capitalized)
    if len(taxonomy_name.split()) == 1 and taxonomy_name[0].isupper():
        return 5  # Genus
    
    # Species (usually two words)
    if len(taxonomy_name.split()) >= 2:
        return 6  # Species
    
    return -1  # Unknown

def process_kraken2_output(input_file, output_file, confidence_threshold=0.0, report_file=None):
    """
    Process Kraken2 output and convert to standardized format.
    
    Args:
        input_file: Path to Kraken2 output file
        output_file: Path to output TSV file
        confidence_threshold: Minimum confidence threshold (not used in basic implementation)
        report_file: Path to Kraken2 report file for enhanced processing
    """
    print(f"Processing Kraken2 output: {input_file}")
    
    # Try to find report file automatically if not provided
    if report_file is None:
        # Look for report file in same directory
        base_dir = os.path.dirname(input_file)
        possible_names = [
            os.path.join(base_dir, "kraken2_report.txt"),
            input_file.replace("_classification.txt", "_report.txt"),
            input_file.replace("_classified.txt", "_report.txt"),
            input_file.replace("_output.txt", "_report.txt")
        ]
        for possible_report in possible_names:
            if os.path.exists(possible_report):
                report_file = possible_report
                print(f"Found report file: {report_file}")
                break
    
    # Parse report file for taxonomic hierarchy information
    taxid_to_info, name_to_info, taxid_to_lineage = parse_kraken2_report(report_file)
    
    output_rows = []
    
    with open(input_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            
            if not line:
                continue
                
            parts = line.split('\t')
            
            if len(parts) != 5:
                print(f"Warning: Line {line_num} has unexpected format: {len(parts)} columns")
                continue
            
            classification_status = parts[0]  # C or U
            sequence_id = parts[1]
            taxonomy_str = parts[2]
            sequence_length = parts[3]
            lca_mapping = parts[4]
            
            # Parse taxonomy
            taxonomy_name, taxid = parse_kraken2_taxonomy(taxonomy_str, taxid_to_info)
            
            # Build standardized taxonomy string using hierarchical information
            if classification_status == 'C':  # Classified
                taxonomy_string, confidence_string, level = build_hierarchical_taxonomy(
                    taxonomy_name, taxid, name_to_info, taxid_to_info, taxid_to_lineage
                )
            else:  # Unclassified
                taxonomy_string, confidence_string, level = build_hierarchical_taxonomy(
                    "Unclassified", 0, name_to_info, taxid_to_info, taxid_to_lineage
                )
            
            output_rows.append((sequence_id, taxonomy_string, confidence_string))
    
    # Create output DataFrame
    out_df = pd.DataFrame(output_rows, columns=["Feature ID", "Taxon", "Confidence"])
    
    # Save to output file
    out_df.to_csv(output_file, sep='\t', index=False)
    print(f"Processed {len(output_rows)} sequences from Kraken2 output")
    print(f"Output saved to: {output_file}")
    
    # Print summary statistics
    if name_to_info:
        print(f"Used taxonomic hierarchy information from report file with {len(name_to_info)} entries")

def process_kraken2_report_for_confidence(report_file):
    """
    Optional: Process Kraken2 report file to extract confidence information.
    This could be used to enhance confidence scores.
    
    Args:
        report_file: Path to Kraken2 report file
        
    Returns:
        dict: Mapping of taxonomy names to confidence scores
    """
    confidence_map = {}
    
    if not os.path.exists(report_file):
        return confidence_map
    
    with open(report_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                percentage = float(parts[0])
                taxonomy_name = parts[5].strip()
                confidence_map[taxonomy_name] = percentage / 100.0
    
    return confidence_map

def main():
    if len(sys.argv) < 3:
        print("Usage: python process_kraken2_output.py <input_file> <output_file> [confidence_threshold] [report_file]")
        print("  input_file: Kraken2 classification output file")
        print("  output_file: Output TSV file in standardized format")
        print("  confidence_threshold: Optional minimum confidence (default: 0.0)")
        print("  report_file: Optional Kraken2 report file for enhanced processing")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    confidence_threshold = float(sys.argv[3]) if len(sys.argv) > 3 else 0.0
    report_file = sys.argv[4] if len(sys.argv) > 4 else None
    
    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} does not exist")
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    process_kraken2_output(input_file, output_file, confidence_threshold, report_file)

if __name__ == "__main__":
    main()
