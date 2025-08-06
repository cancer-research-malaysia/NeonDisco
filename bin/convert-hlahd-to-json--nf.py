#!/usr/bin/env python3
"""
Convert HLA-HD TSV output to arcasHLA-compatible JSON format.

Usage: convert_hlahd_to_json.py <hlahd_result_file> <sample_name>

Input format (HLA-HD TSV):
A	HLA-A*11:02:01	HLA-A*74:02:01
B	HLA-B*40:01:02	HLA-B*54:01:01
C	HLA-C*07:02:01	HLA-C*01:217:01

Output: JSON compatible with arcasHLA format
"""

import json
import sys
import os

def convert_hlahd_to_json(result_file, sample_name):
    """Convert HLA-HD TSV output to arcasHLA JSON format."""
    
    output = {}
    
    try:
        if not os.path.exists(result_file):
            raise FileNotFoundError(f"HLA-HD result file not found: {result_file}")
            
        with open(result_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    parts = line.split('\t')
                    if len(parts) >= 3:  # Locus + at least 2 alleles
                        locus = parts[0]  # A, B, or C
                        if locus in ['A', 'B', 'C']:
                            alleles = []
                            # Add all alleles for this locus (typically 2)
                            for i in range(1, len(parts)):
                                if parts[i].startswith('HLA-'):
                                    # Remove HLA- prefix to match expected format
                                    allele = parts[i].replace('HLA-', '')
                                    alleles.append(allele)
                            
                            # Only add locus if we found alleles
                            if alleles:
                                output[locus] = alleles
        
        return json.dumps(output, indent=2)
        
    except Exception as e:
        # If conversion fails, create empty output
        error_output = {"error": str(e)}
        return json.dumps(error_output, indent=2)

def main():
    if len(sys.argv) != 3:
        print("Usage: convert-hlahd-to-json--nf.py <hlahd_result_file> <sample_name>", file=sys.stderr)
        sys.exit(1)
    
    result_file = sys.argv[1]
    sample_name = sys.argv[2]
    
    json_output = convert_hlahd_to_json(result_file, sample_name)
    print(json_output)

if __name__ == "__main__":
    main()
