#!/usr/bin/env python3
"""
Helper script for QIIME2 operations in benchmarking
"""

import sys
import pandas as pd
from qiime2 import Artifact

def convert_qza_to_tsv(qza_file, output_tsv):
    """Convert QIIME2 QZA file to TSV format"""
    try:
        df = Artifact.load(qza_file).view(pd.DataFrame)
        df.to_csv(output_tsv, sep='\t', index=True)
        print(f"Successfully converted {qza_file} to {output_tsv}")
        return True
    except Exception as e:
        print(f"Error converting {qza_file}: {e}")
        return False

def import_fasta_to_qza(fasta_file, qza_file):
    """Import FASTA file to QIIME2 QZA format"""
    try:
        artifact = Artifact.import_data('FeatureData[Sequence]', fasta_file)
        artifact.save(qza_file)
        print(f"Successfully imported {fasta_file} to {qza_file}")
        return True
    except Exception as e:
        print(f"Error importing {fasta_file}: {e}")
        return False

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python qiime_helper.py <command> <args...>")
        print("Commands:")
        print("  convert_qza_to_tsv <qza_file> <output_tsv>")
        print("  import_fasta_to_qza <fasta_file> <qza_file>")
        sys.exit(1)
    
    command = sys.argv[1]
    
    if command == "convert_qza_to_tsv":
        if len(sys.argv) != 4:
            print("Usage: python qiime_helper.py convert_qza_to_tsv <qza_file> <output_tsv>")
            sys.exit(1)
        success = convert_qza_to_tsv(sys.argv[2], sys.argv[3])
        sys.exit(0 if success else 1)
    
    elif command == "import_fasta_to_qza":
        if len(sys.argv) != 4:
            print("Usage: python qiime_helper.py import_fasta_to_qza <fasta_file> <qza_file>")
            sys.exit(1)
        success = import_fasta_to_qza(sys.argv[2], sys.argv[3])
        sys.exit(0 if success else 1)
    
    else:
        print(f"Unknown command: {command}")
        sys.exit(1)
