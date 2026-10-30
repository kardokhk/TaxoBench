#!/usr/bin/env python3

import os
import sys
import pandas as pd

def convert_taxonomy_with_all_confidences(input_tsv, output_tsv, ncbi_style_taxonomy=False):
    # Read input file
    df = pd.read_csv(input_tsv, sep='\t')

    # Set prefixes
    domain_prefix = "k__" if ncbi_style_taxonomy else "d__"
    levels = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    prefixes = [domain_prefix, 'p__', 'c__', 'o__', 'f__', 'g__', 's__']

    output_rows = []

    for _, row in df.iterrows():
        sample_id = row['sequence_id']
        taxonomy_parts = []
        confidence_parts = []

        unclassified_parent = False
        for level, prefix in zip(levels, prefixes):
            name = str(row[level])
            conf = float(row[f"{level}_confidence"])

            # if parent was unclassified, or current level is unclassified, then force unclassified
            # Handle both "Unclassified" and "Unknown" designations for unclassified sequences
            if unclassified_parent or name.lower() in ['unclassified', 'unknown']:
                taxonomy_parts.append(f"{prefix}")
                confidence_parts.append("0.0000")
                unclassified_parent = True
            else:
                if ncbi_style_taxonomy and prefix == 's__':
                    genus = str(row['genus'])
                    species_parts = name.split()
                    if species_parts[0].lower() == genus.lower():
                        species_clean = " ".join(species_parts[1:])
                    else:
                        species_clean = name
                    taxonomy_parts.append(f"{prefix}{species_clean}")
                else:
                    taxonomy_parts.append(f"{prefix}{name}")
                confidence_parts.append(f"{conf:.4f}")

        taxonomy_str = "; ".join(taxonomy_parts)
        confidence_str = "; ".join(confidence_parts)

        output_rows.append((sample_id, taxonomy_str, confidence_str))

    out_df = pd.DataFrame(output_rows, columns=["Feature ID", "Taxon", "Confidence"])
    out_df.to_csv(output_tsv, sep='\t', index=False)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python process_taxonova_output.py <input_tsv> <output_tsv> <ncbi_style_taxonomy>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    ncbi_style_taxonomy = sys.argv[3].lower() == "true"

    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)

    convert_taxonomy_with_all_confidences(input_file, output_file, ncbi_style_taxonomy)
